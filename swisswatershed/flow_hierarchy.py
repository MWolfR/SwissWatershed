import pandas, numpy, os
import geopandas as gpd


def recursive_build(curr_tezgnr, subtree_data, dict_children, dict_parent):
    lst_children = []
    while len(subtree_data) > 0:
        child_idx = (subtree_data.H2 - subtree_data.H1).argmax()
        child_tezgnr = subtree_data.iloc[child_idx].name
        lst_children.append(child_tezgnr)
        dict_parent[child_tezgnr] = curr_tezgnr

        child_h1 = subtree_data.iloc[child_idx].H1
        child_h2 = subtree_data.iloc[child_idx].H2
        child_data = subtree_data.loc[(subtree_data.H1 >= child_h1)
                                      & (subtree_data.H2 < child_h2)]
        recursive_build(child_tezgnr, child_data, dict_children, dict_parent)
        subtree_data = subtree_data.drop(child_data.index).drop(child_tezgnr)
    dict_children[curr_tezgnr] = lst_children


def build_parent_children_dict(data):
    dict_c = {}
    dict_p = {}
    recursive_build(-1, data, dict_c, dict_p)
    parents = pandas.DataFrame.from_dict(dict_p, orient="index")[0]
    return parents, dict_c


def expand_generations(parents):
    lst_gen =[parents.reset_index().rename(columns={"index": 0, 0: 1}).set_index(parents.index)]
    i = 2
    while (parents > -1).any():
        parents = pandas.Series(lst_gen[0][1][parents].values, index=parents.index, name=i)
        i += 1
        lst_gen.append(parents)
    return pandas.concat(lst_gen, axis=1)


def find_neighbors(flow_integrator, max_dist=20000):
    from scipy.spatial import KDTree

    c_xy = numpy.hstack(flow_integrator._centroids.apply(lambda _x: _x.xy).values).transpose()
    tree = KDTree(c_xy)

    pairs = numpy.vstack(list(tree.query_pairs(max_dist)))
    A = flow_integrator._geometry.iloc[pairs[:, 0]]
    B = flow_integrator._geometry.iloc[pairs[:, 1]]
    passes = A.reset_index().geometry.combine(B.reset_index().geometry, lambda _a, _b: _a.intersects(_b))
    pair_idx = pandas.DataFrame({"a": A.index[passes], "b": B.index[passes]})
    return pair_idx


class SourceData(object):
    LOC_BASINS = "EZG_Gewaesser.gpkg"
    LOC_OUTFLOW = os.path.join("EZG_Gewaesser.gdb", "a00000032.gdbtable")
    LOC_CITY_NAMES = "a00000050.gdbtable"
    LOC_CITY_LOCATION = "a00000055.gdbtable"
    col_id = "TEILEZGNR"
    col_name = "LANGTEXT"
    col_uuid = "PLZO_OS_UUID"

    def __init__(self, root, root_cities=None):
        self._basins = gpd.read_file(os.path.join(root, SourceData.LOC_BASINS)).set_index(SourceData.col_id)
        self._outflow = gpd.read_file(os.path.join(root, SourceData.LOC_OUTFLOW)).set_index(SourceData.col_id)
        if root_cities is not None:
            self._city_lo1 = gpd.read_file(os.path.join(root_cities, SourceData.LOC_CITY_NAMES)).set_index(SourceData.col_name)
            self._city_lo2 = gpd.read_file(os.path.join(root_cities, SourceData.LOC_CITY_LOCATION)).set_index(SourceData.col_uuid)
        else:
            self._city_lo1 = None
            self._city_lo2 = None


class FlowIntegrator(object):
    def __init__(self, tbl_basins, tbl_flowout):
        self._basins = tbl_basins
        self._flowout = tbl_flowout
        self._centroids = tbl_basins.geometry.apply(lambda _x: _x.centroid)
        self._geometry = tbl_basins.geometry
        self._flowout = tbl_flowout.geometry
        self._flowout = pandas.concat([
            self._flowout,
            self._centroids[self._centroids.index.difference(self._flowout.index)]
        ], axis=0)
        
        print("Building hierarchy based on outflow...")
        parents, children = build_parent_children_dict(self._basins)
        print("Done!")
        parents[-1] = -1
        self._gens = expand_generations(parents)
    
    def integrate_length(self, row, stop_idx=-1, count_first=True):
        lens = []
        if count_first: pt_fr = self._centroids[row.name]
        else: pt_fr = self._flowout[row.name]
        for _idx in row.values:
            if _idx == -1: break
            pt_to = self._flowout[_idx]
            if _idx == stop_idx:
                break
            lens.append(pt_fr.distance(pt_to))
            pt_fr = pt_to
        return numpy.sum(lens), pt_to
    
    def cumulative_length_for_pair(self, idx1, idx2):
        lca_idx = numpy.nonzero(numpy.in1d(self._gens.loc[idx1], self._gens.loc[idx2]))[0][0]
        lca = self._gens.loc[idx1][lca_idx]
        not_ancestor = (lca != idx1) and (lca != idx2)
        a, pt_a = self.integrate_length(self._gens.loc[idx1], lca, count_first=not_ancestor)
        b, pt_b = self.integrate_length(self._gens.loc[idx2], lca, count_first=not_ancestor)
        if lca == -1:
            return a + b + pt_a.distance(pt_b)
        assert pt_a.distance(pt_b) == 0.0
        return a + b


class WatershedSeparation(object):
    def __init__(self, flow_integrator, pair_max_dist=20000):
        self._flow = flow_integrator
        print("Finding neighboring areas...")
        self.pairs = find_neighbors(flow_integrator, max_dist=pair_max_dist)
        print("Done!")
        print("Calculating separation of {0} pairs...".format(len(self.pairs)))
        self._payload = self.calculate_separation(self.pairs)
        print("Done!")
        self.payload = self._payload

    @staticmethod
    def facultative_conversion(obj):
        import shapely
        from shapely import ops
        if isinstance(obj, shapely.MultiLineString):
            return ops.linemerge(obj)
        return obj

    def calculate_separation(self, pairs):
        pair_sep = pairs.apply(lambda _x: self._flow.cumulative_length_for_pair(*_x), axis=1)
        pair_geom = pairs.applymap(lambda _x: self._flow._geometry[_x])
        pair_inters = pair_geom.apply(lambda _x: _x[0].intersection(_x[1]), axis=1)
        pair_inters = pair_inters.apply(self.facultative_conversion)
        pair_sep.name = "separation"
        pair_inters.name = "border"
        return pandas.concat([pairs, pair_sep, pair_inters], axis=1)

    def filter_lower_res(self, str_prop):
        tbl = self.payload
        series = self._flow._basins[str_prop]
        valid = (series[tbl["a"]].values != series[tbl["b"]].values)
        self.payload = tbl.loc[valid]
        return self
    
    def filter_min_val(self, min_val):
        tbl = self.payload
        valid = tbl["separation"] > min_val
        self.payload = tbl.loc[valid]
        return self
    
    def filter_reset(self):
        self.payload = self._payload
        return self

