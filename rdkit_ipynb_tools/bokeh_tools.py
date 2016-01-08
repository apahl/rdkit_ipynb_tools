#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#####
Tools
#####

*Created on 2015-12-12 by A. Pahl*

Bokeh plotting functionality for Mol_lists.
"""

import colorsys

import numpy as np

from bokeh.plotting import figure, ColumnDataSource
from bokeh.io import output_notebook, show
from bokeh.models import HoverTool

AVAIL_COLORS = ["cornflowerblue", "firebrick", "goldenrod", "aqua", "brown", "chartreuse", "darkmagenta"
                "aquamarine", "blue", "red", "blueviolet", "darkorange", "forestgreen", "lime"]

output_notebook()


class ColorScale():
    """Used for continuous coloring."""

    def __init__(self, num_values, val_min, val_max):
        self.num_values = num_values
        self.num_val_1 = num_values - 1
        self.value_min = val_min
        self.value_max = val_max
        self.value_range = self.value_max - self.value_min
        self.color_scale = []
        hsv_tuples = [(0.35 + ((x * 0.65) / (self.num_val_1)), 0.9, 0.9) for x in range(self.num_values)]
        rgb_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), hsv_tuples)
        for rgb in rgb_tuples:
            rgb_int = [int(255 * x) for x in rgb]
            self.color_scale.append('#{:02x}{:02x}{:02x}'.format(*rgb_int))

    def __call__(self, value, reverse=False):
        """return the color from the scale corresponding to the place in the value_min ..  value_max range"""
        pos = int(((value - self.value_min) / self.value_range) * self.num_val_1)

        if reverse:
            pos = self.num_val_1 - pos

        return self.color_scale[pos]


class Chart():
    """A Bokeh Plot."""

    def __init__(self, kind="scatter", **kwargs):
        self.kind = kind
        self.height = kwargs.get("height", 450)
        self.title = kwargs.get("title", "Bokeh Plot")

        self.plot = figure(plot_height=self.height, title=self.title, tools="pan,wheel_zoom,box_zoom,reset,resize")


    def add_data(self, d, x, y, **kwargs):
        colors = "cornflowerblue"
        series_by = kwargs.get("series_by", None)
        color_by = kwargs.get("color_by", None)
        tooltip = get_tooltip(x, y,
                              kwargs.get("pid", None),
                              series_by,
                              color_by,
                              kwargs.get("tooltip", None))

        self.plot.add_tools(tooltip)
        self.plot.xaxis.axis_label = x
        self.plot.yaxis.axis_label = y

        radius = kwargs.get("radius", kwargs.get("r", 10))

        if series_by:
            avail_colors = AVAIL_COLORS
            series_colors = {item: 0 for item in d[series_by] if item is not None and item != np.nan}
            if len(series_colors) > len(avail_colors):
                raise LookupError("Too many series values (more than available colors {}).".format(len(avail_colors)))

            for idx, series_item in enumerate(series_colors):
                series_colors[series_item] = avail_colors[idx]

            colors = []
            for val in d[series_by]:
                colors.append(series_colors.get(val, "black"))

            d["colors"] = colors

        elif color_by:
            color_by_min = min(d[color_by])
            color_by_max = max(d[color_by])
            color_scale = ColorScale(20, color_by_min, color_by_max)
            colors = []
            for val in d[color_by]:
                if val is not None and val != np.nan:
                    colors.append(color_scale(val))
                else:
                    colors.append("black")

            d["colors"] = colors

        d["x"] = d[x]
        d["y"] = d[y]
        self.plot.circle(x, y, size=radius, color=colors, source=ColumnDataSource(d))

    def show(self):
        show(self.plot)


def get_tooltip(x, y, pid=None, series_by=None, color_by=None, tooltip=None):
    if pid:
        pid_tag = '<span style="font-size: 13px; color: #000000;">{pid}: @{pid}</span><br>'.format(pid=pid)
    else:
        pid_tag = ""

    if series_by:
        series_tag = '<span style="font-size: 13px; color: #000000;"><b>{series_by}: @{series_by}</b>&nbsp;&nbsp;</span><span style="font-size: 14px; color: @colors;">&#9899</span><br>'.format(series_by=series_by)
        color_tag = ""
    elif color_by:
        series_tag = ""
        color_tag = '<span style="font-size: 13px; color: #000000;"><br>{color_by}: @{color_by}&nbsp;&nbsp;</span><span style="font-size: 14px; color: @colors;">&#9899</span>'.format(color_by=color_by)
    else:
        color_tag = ""
        series_tag = ""

    if tooltip == "struct":
        templ = HoverTool(
            tooltips="""
            <div>
                <div style="width: 200px; height: 200px;">
                    <img
                        src="data:image/png;base64,@mol"

                        border="2" alt="Mol"
                    ></img>
                </div>
                <div>
                    {series_tag}{pid_tag}
                    <span style="font-size: 13px; color: #000000;">{x}: @x<br>
                    {y}: @y</span>{color_tag}
                </div>
            </div>
            """.format(pid_tag=pid_tag, series_tag=series_tag, x=x, y=y)
        )
    else:
        templ = HoverTool(
            tooltips="""
            <div>
                {series_tag}{pid_tag}
                <span style="font-size: 13px; color: #000000;">{x}: @x<br>
                {y}: @y</span>{color_tag}
            </div>
            """.format(pid_tag=pid_tag, series_tag=series_tag, color_tag=color_tag, x=x, y=y)
        )
        # templ = HoverTool(tooltips=[(x, "@x"), (y, "@y")])

    return templ


def guess_id_prop(prop_list):  # try to guess an id_prop
    for prop in prop_list:
        if prop.lower().endswith("id"):
            return prop
    return None


def cpd_scatter(df, x, y, r=7, pid=None, **kwargs):
    """Predefined Plot #1.
    Quickly plot an RDKit Pandas dataframe or a molecule dictionary with structure tooltips."""

    if not pid:
        if isinstance(df, dict):
            prop_list = df.keys()
        else:
            prop_list = [df.index.name]
            prop_list.extend(df.columns.values)

        pid = guess_id_prop(prop_list)

    title = kwargs.get("title", "Compound Scatter Plot")
    scatter = Chart(title=title, r=r)
    scatter.add_data(df, x, y, pid=pid, **kwargs)
    return scatter.show()
