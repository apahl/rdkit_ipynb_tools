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
        self.tooltip = get_tooltip(kwargs.get("tooltip", None))

        self.plot = figure(plot_height=self.height, tools=[self.tooltip])


    def add_data(self, d, x, y, **kwargs):
        self.plot.circle('x', 'y', size=20, source=ColumnDataSource(d))

    def show(self):
        show(self.plot)


def get_tooltip(x, y, tooltip):
    if tooltip == "struct":
        templ = HoverTool(
            tooltips="""
            <div>
                <div>
                    <img
                        src="@mol" height="200" width="200"
                        style="float: left; margin: 0px 15px 15px 0px;"
                        border="2"
                    ></img>
                </div>
                <div>
                    <span style="font-size: 12px; color: #000000;">{x}: @x<br>
                    {y}: @y</span>
                </div>
            </div>
            """.format(x=x, y=y)
        )
    else:
        templ = HoverTool(tooltips=[(x, "@x"), (y, "@y")])

    return templ
