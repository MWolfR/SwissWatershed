import numpy
import pandas
import shapely

from matplotlib import pyplot as plt
from matplotlib import cm


DEFAULT_COL = {
    "scale": 0.5,
    "slope": 1.0/400.,
    "offset": -0.75,
    "clip": [0.0, 1.0]
}
DEFAULT_LW = {
    "scale": 0.5,
    "slope": 1.0/175,
    "offset": -0.1,
    "clip": [0.1, 3.0]
}

def make_processor(**kwargs):
    def process(val):
        scale = kwargs.get("scale", None)
        clip = kwargs.get("clip", None)
        if scale == "log":
            val = numpy.log10(val)
        elif isinstance(scale, float):
            val = val ** scale
        val = val * kwargs.get("slope", 1.0)
        val = val + kwargs.get("offset", 0.0)
        if clip is not None:
            val = numpy.maximum(numpy.minimum(val, clip[1]), clip[0])
            
        return val
    return process

def _plot_hlp(plt_specs):
    col = plt_specs["color"]
    lw = plt_specs["lw"]
    coords = plt_specs["border"]
    if isinstance(coords, shapely.LineString):
        plt.plot(*coords.xy, color=col, lw=lw)

def apply_filters(watershed, **kwargs):
    watershed.filter_reset()
    if "min_val" in kwargs:
        watershed.filter_min_val(kwargs["min_val"])
    if "lower_res" in kwargs:
        watershed.filter_lower_res(kwargs["lower_res"])

def plot_according_to_specs(watershed, config):
    print("Plotting the data...")
    proc_col = make_processor(**config.get("color", DEFAULT_COL))
    proc_lw = make_processor(**config.get("linewidth", DEFAULT_LW))
    cmap = cm.__dict__[config.get("colormap", "copper")]
    apply_filters(watershed, **config.get("filters", {"min_val": 0}))

    specs = pandas.concat([watershed.payload["separation"].apply(proc_col).apply(cmap),
                           watershed.payload["separation"].apply(proc_lw),
                           watershed.payload["border"]],
                          keys=["color", "lw", "border"], axis=1)
    fig = plt.figure(figsize=config.get("figsize", (12, 8)))
    ax = fig.add_axes([0.01, 0.01, 0.85, 0.98])
    _ = specs.apply(_plot_hlp, axis=1)
    ax.set_frame_on(False); ax.set_xticks([]); ax.set_yticks([])
    plt.axis("equal")

    ax = fig.add_axes([0.9, 0.25, 0.05, 0.5])
    min_val = 100
    max_val = numpy.percentile(watershed.payload["separation"], 99)
    vals = numpy.linspace(min_val, max_val, 100)
    val_centers = 0.5 * (vals[:-1] + vals[1:])
    cols = proc_col(val_centers)
    lw = proc_lw(val_centers)
    vals = vals / 1000

    for a, b, c, l in zip(vals[:-1], vals[1:], cols, lw):
        plt.plot([0, 0], [a, b], color=cm.copper(c), lw=l)
    ax.set_frame_on(False); ax.set_xticks([])
    ax.yaxis.tick_right()
    ax.set_ylabel("Km")
    print("Done!")
    return fig

def main():
    import sys, json
    from flow_hierarchy import FlowIntegrator, WatershedSeparation, SourceData
    fn_cfg = sys.argv[1]
    with open(fn_cfg, "r") as fid:
        config = json.load(fid)
    data = SourceData(config["data"]["root"])
    flow = FlowIntegrator(data._basins, data._outflow)
    watershed = WatershedSeparation(flow)
    fig = plot_according_to_specs(watershed, config.get("plotting", {}))
    fig.savefig(config.get("output", "watershed.svg"))

if __name__ == "__main__":
    main()
    