#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#####
Bokeh Tools
#####

*Created on 2015-12-12 by A. Pahl*

Bokeh plotting functionality for Mol_lists.
"""

import colorsys

import numpy as np

from bokeh.plotting import figure, ColumnDataSource
from bokeh.io import output_notebook, show
from bokeh.models import HoverTool

AVAIL_COLORS = ["#1F77B4", "firebrick", "goldenrod", "aqua", "brown", "chartreuse", "darkmagenta"
                "aquamarine", "blue", "red", "blueviolet", "darkorange", "forestgreen", "lime"]
# AVAIL_MARKERS: circle, diamond, triangle, square, inverted_triangle, asterisk,
#                circle_cross, circle_x, cross, diamond_cross, square_cross, square_x, asterisk, diamond

output_notebook()


class ColorScale():

    def __init__(self, num_values, val_min, val_max, middle_color="yellow", reverse=False):
        self.num_values = num_values
        self.num_val_1 = num_values - 1
        self.value_min = val_min
        self.value_max = val_max
        self.reverse = reverse
        self.value_range = self.value_max - self.value_min
        self.color_scale = []
        if middle_color.startswith("y"):  # middle color yellow
            hsv_tuples = [(0.0 + ((x * 0.35) / (self.num_val_1)), 0.99, 0.9) for x in range(self.num_values)]
            self.reverse = not self.reverse
        else:  # middle color blue
            hsv_tuples = [(0.35 + ((x * 0.65) / (self.num_val_1)), 0.9, 0.9) for x in range(self.num_values)]
        rgb_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), hsv_tuples)
        for rgb in rgb_tuples:
            rgb_int = [int(255 * x) for x in rgb]
            self.color_scale.append('#{:02x}{:02x}{:02x}'.format(*rgb_int))

        if self.reverse:
            self.color_scale.reverse()

    def __call__(self, value):
        """return the color from the scale corresponding to the place in the value_min .. value_max range"""
        pos = int(((value - self.value_min) / self.value_range) * self.num_val_1)

        return self.color_scale[pos]


    def legend(self):
        """Return the value_range and a list of tuples (value, color) to be used in a legend."""
        legend = []
        for idx, color in enumerate(self.color_scale):
            val = self.value_min + idx / self.num_val_1 * self.value_range
            legend.append((val, color))

        return legend


class Chart():
    """A Bokeh Plot."""

    def __init__(self, kind="scatter", **kwargs):
        self.data = {}
        self.kind = kind
        self.height = kwargs.get("height", 450)
        self.title = kwargs.get("title", "Bokeh Plot")
        self.position = kwargs.get("position", kwargs.get("pos", "top_left"))

        self.series_counter = 0
        self.tools_added = False

        self.plot = figure(plot_height=self.height, title=self.title, tools="pan,wheel_zoom,box_zoom,reset,resize,save")


    def _add_series(self, x, y, series, size, source=None):
        color = AVAIL_COLORS[self.series_counter]

        if self.series_counter == 0:
            self.plot_type = self.plot.circle

        elif self.series_counter == 1:
            self.plot_type = self.plot.diamond
            size += 3  # diamonds appear smaller than circles of the same size
        elif self.series_counter == 2:
            self.plot_type = self.plot.triangle
        elif self.series_counter == 4:
            self.plot_type = self.plot.inverted_triangle
        elif self.series_counter == 5:
            self.plot_type = self.plot.asterisk
        elif self.series_counter == 6:
            self.plot_type = self.plot.circle_cross
        elif self.series_counter == 7:
            self.plot_type = self.plot.circle_x
        elif self.series_counter == 8:
            self.plot_type = self.plot.cross
        elif self.series_counter == 9:
            self.plot_type = self.plot.diamond_cross
        elif self.series_counter == 10:
            self.plot_type = self.plot.square_cross
        elif self.series_counter == 11:
            self.plot_type = self.plot.square_x
        else:
            self.plot_type = self.plot.asterisk

        self.plot_type(x, y, legend=series, size=size, color=color, source=source)
        self.plot.legend.orientation = self.position

        self.series_counter += 1
        if self.series_counter >= len(AVAIL_COLORS):
            print("* series overflow, starting again.")
            self._add_series = 0


    def add_data(self, d, x, y, **kwargs):
        colors = "#1F77B4"
        series = kwargs.get("series", None)
        if series is not None:
            series_by = "Series"
        else:
            series_by = kwargs.get("series_by", None)

        color_by = kwargs.get("color_by", None)

        self.plot.xaxis.axis_label = x
        self.plot.yaxis.axis_label = y

        tooltip = get_tooltip(x, y,
                              kwargs.get("pid", None),
                              series,
                              series_by,
                              color_by,
                              kwargs.get("tooltip", None))

        self.plot.add_tools(tooltip)


        size = kwargs.get("radius", kwargs.get("r", kwargs.get("size", kwargs.get("s", 10))))
        reverse = kwargs.get("invert", False)

        if series:
            d["x"] = d[x]
            d["y"] = d[y]
            d["series"] = [series] * len(d[x])

            self._add_series(x, y, series, size=size, source=ColumnDataSource(d))

        elif series_by:
            series_keys = set()
            for idx, item in enumerate(d[series_by]):
                if item is None:
                    d[series_by][idx] = "None"
                elif item is np.nan:
                    d[series_by][idx] = "NaN"

                series_keys.add(d[series_by][idx])

            for series in series_keys:
                d_series = {x: [], y: [], "series": []}
                for idx, el in enumerate(d[x]):
                    if d[series_by][idx] == series:
                        d_series[x].append(d[x][idx])
                        d_series[y].append(d[y][idx])
                        d_series["series"].append(d[series_by][idx])

                d_series["x"] = d_series[x]
                d_series["y"] = d_series[y]

                self._add_series(x, y, series, size=size, source=ColumnDataSource(d_series))


        elif color_by:
            color_by_min = min(d[color_by])
            color_by_max = max(d[color_by])
            color_scale = ColorScale(20, color_by_min, color_by_max, reverse)
            colors = []
            for val in d[color_by]:
                if val is not None and val != np.nan:
                    colors.append(color_scale(val))
                else:
                    colors.append("black")

            d["colors"] = colors
            d["x"] = d[x]
            d["y"] = d[y]
            self.plot.circle(x, y, size=size, color=colors, source=ColumnDataSource(d))

        else:
            d["x"] = d[x]
            d["y"] = d[y]
            self.plot.circle(x, y, size=size, source=ColumnDataSource(d))


    def show(self):
        show(self.plot)


def get_tooltip(x, y, pid=None, series=None, series_by=None, color_by=None, tooltip=None):
    if pid:
        pid_tag = '<span style="font-size: 13px; color: #000000;">{pid}: @{pid}</span><br>'.format(pid=pid)
    else:
        pid_tag = ""

    if series_by:
        series_tag = '<span style="font-size: 13px; color: #000000;"><b>{series_by}: @series</b>&nbsp;&nbsp;</span><br>'.format(series_by=series_by)
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
            """.format(pid_tag=pid_tag, series_tag=series_tag, color_tag=color_tag, x=x, y=y)
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
