#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#####
Tools
#####

*Created on 2015-12-12 by A. Pahl*

Bokeh plotting functionality for Mol_lists.
"""

from bokeh.plotting import figure, ColumnDataSource
from bokeh.io import output_notebook, show
from bokeh.models import HoverTool

output_notebook()


class Chart():
    """A Bokeh Plot."""

    def __init__(self, kind="scatter", **kwargs):
        self.kind = kind
        self.height = kwargs.get("height", 450)
        self.title = kwargs.get("title", "Bokeh Plot")

        self.plot = figure(plot_height=self.height, title=self.title, tools="pan,wheel_zoom,box_zoom,reset,resize")


    def add_data(self, d, x, y, **kwargs):
        self.tooltip = get_tooltip(x, y,
                                   kwargs.get("pid", None),
                                   kwargs.get("series", None),
                                   kwargs.get("tooltip", None))

        self.plot.add_tools(self.tooltip)
        self.plot.xaxis.axis_label = x
        self.plot.yaxis.axis_label = y

        self.radius = kwargs.get("radius", kwargs.get("r", 10))

        d["x"] = d[x]
        d["y"] = d[y]
        self.plot.circle(x, y, size=self.radius, source=ColumnDataSource(d))

    def show(self):
        show(self.plot)


def get_tooltip(x, y, pid=None, series=None, tooltip=None):
    if pid:
        pid_tag = '<span style="font-size: 12px; color: #000000;">{pid}: @{pid}<br>'.format(pid=pid)
    else:
        pid_tag = ""

    if series:
        series_tag = '<span style="font-size: 12px; color: #000000;"><b>series: {series}</b><br>'.format(series)
    else:
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
                    <span style="font-size: 12px; color: #000000;">{x}: @x<br>
                    {y}: @y</span>
                </div>
            </div>
            """.format(pid_tag=pid_tag, series_tag=series_tag, x=x, y=y)
        )
    else:
        templ = HoverTool(
            tooltips="""
            <div>
                {series_tag}{pid_tag}
                <span style="font-size: 12px; color: #000000;">{x}: @x<br>
                {y}: @y</span>
            </div>
            """.format(pid_tag=pid_tag, series_tag=series_tag, x=x, y=y)
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
