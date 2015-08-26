# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 08:39:28 2015

@author: Axel Pahl
@title: hc_tools.py
"""

# 1. stdlib imports
import time
import string
import json

# 2. third-party imports

# 3. internal project imports
from . import tools

# 4. IPython imports
# from IPython.html import widgets
from IPython.display import HTML, display

HIGHCHARTS = """
<script src="lib/highcharts.js"></script>
<script src="lib/highcharts-more.js"></script>
<script src="lib/modules/exporting.js"></script>
"""

CHART_TEMPL = """<div id="container_${id}" style="height: ${height}px"></div>
<script>
$$(function () {
    $$('#container_${id}').highcharts(

$chart

    );
});
</script>
"""

CHART_KINDS = ["scatter", "column"]
TOOLTIP_OPTIONS = "struct"

print("- loading highcharts...")
display(HTML(HIGHCHARTS))


class Chart():
    def __init__(self, kind="scatter", **kwargs):
        if not kind in CHART_KINDS:
            raise ValueError("{} is not a supported chart kind ({})".format(kind, CHART_KINDS))

        self.kind = kind
        self.height= kwargs.get("height", 400)
        radius = kwargs.get("r", kwargs.get("radius", 5)) # accept "r" or "radius" for this option
        self.legend = kwargs.get("legend", None)
        self.chart_id = time.strftime("%y%m%d%H%M%S")
        self.chart = {}
        self.chart["title"] = {"text": kwargs.get("title", "{} plot".format(self.kind))}
        self.chart["subtitle"] = {"text": kwargs.get("subtitle")}
        self.chart["series"] = []
        self.chart["plotOptions"] = {"scatter": {"marker": {"radius": radius}}}
        
    
    def _structure_tooltip(self, i):
        return str(self.dpid[i])+"<br>"+str(self.dmol[i])


    def _data_columns(self):
        """Generate the data for the Column plot"""
        data = []
        cats = []

        for i in range(self.dlen):
            cats.append(str(self.dx[i]))
            data.append(float(self.dy[i]))
        
        self.chart["series"].append({"name": self.arg_y, "data": data})
        self.chart["xAxis"]["categories"] = cats
        

    def _data_tuples(self, d):
        """Generate the data tuples required for Highcharts scatter plot."""
        data = []
        dx = d["x"]
        dy = d["y"]
        if self.arg_z:
            dz = d["z"]
        if self.arg_pid or self.arg_struct:
            dpid = d["id"]

        for i in range(len(dx)):
            tmp_d = {"x": float(dx[i]), "y": float(dy[i])}
            if self.arg_z:
                tmp_d["z"] = float(dz[i])
            if self.arg_pid:
                tmp_d["id"] = str(dpid[i])

            data.append(tmp_d)
        
        return data
    

    def _data_series(self):
        # [{"name": "A", "data": [{"x": 1, "y": 2}, {"x": 2, "y": 3}]},
        #  {"name": "B", "data": [{"x": 2, "y": 3}, {"x": 3, "y": 4}]}]
        self.arg_z = None  # not implemented yet
        series = []

        names = set(str(c) for c in self.dcolor_by)
        color_series_x = {name: [] for name in names}
        color_series_y = {name: [] for name in names}
        if self.arg_pid:
            color_series_id = {name: [] for name in names}
        if self.arg_struct:
            color_series_mol = {name: [] for name in names}
        
        for i in range(self.dlen):
            col_by_str = str(self.dcolor_by[i])
            color_series_x[col_by_str].append(float(self.dx[i]))
            color_series_y[col_by_str].append(float(self.dy[i]))
            if self.arg_struct:
                color_series_mol[col_by_str].append(self._structure_tooltip(i))
            elif self.arg_pid:
                color_series_id[col_by_str].append(str(self.dpid[i]))

        for name in names:
            tmp_d = {"x": color_series_x[name], "y": color_series_y[name]}
            if self.arg_struct:
                tmp_d["id"] = color_series_mol[name]
            elif self.arg_pid:
                tmp_d["id"] = color_series_id[name]

            series_dict = {"name": name}
            series_dict["data"] = self._data_tuples(tmp_d)
            series.append(series_dict)
        
        return series
        


    def add_data(self, d, x="x", y="y", pid=None, z=None, color_by=None, **kwargs):
        """Add the data to the chart.
        d is the input dictionary, x, y [, and z] are the keys for the properties to plot.
        pid is the optional key to a (compound) id to be displayed in the tooltip."""
        if not x in d or not y in d:
            raise KeyError("'{x}' and '{y}' are required parameters for scatter plot, but could not all be found in dict.".format(x=x, y=y))
        
        if len(d[x]) != len(d[y]):
            raise ValueError("'{x}' and '{y}' must have the same length.".format(x=self.arg_x, y=self.arg_y))
        
        self.arg_x = x
        self.arg_y = y
        self.arg_z = z
        self.arg_color_by = color_by
        self.arg_pid = pid
        
        self.dx = list(d[x])
        self.dy = list(d[y])
        self.dlen = len(self.dx)
        
        if self.arg_pid:
            # pandas data series and pid == index
            if not isinstance(d[x], list) and self.arg_pid == d.index.name:
                self.dpid = list(d.index)
            else:
                self.dpid = list(d[pid])

            if self.dlen != len(self.dpid):
                raise ValueError("'{x}' and '{pid}' must have the same length.".format(x=self.arg_x, pid=self.arg_pid))
        else:
            self.arg_pid = None

        # plot-specific options
        if self.kind in ["scatter"]:
            self.arg_tooltip = kwargs.get("tooltip", "")
            if self.arg_tooltip not in TOOLTIP_OPTIONS:
                print("- unknown tooltip option {}, setting to empty.".format(self.arg_tooltip))
                self.arg_tooltip = ""
            self.arg_struct = "struct" in self.arg_tooltip
            self.arg_mol_col = kwargs.get("mol_col", "mol")
            if self.arg_struct:
                self.dmol = list(d[self.arg_mol_col])
                print("- try to display structure tooltips")

        self.chart["credits"] = {'enabled': False}
        self.chart["xAxis"] = {"title": {"enabled": True, "text": self.arg_x}}
        self.chart["yAxis"] = {"title": {"enabled": True, "text": self.arg_y}}
        
        if self.kind == "scatter":
            # defining the tooltip
            self.chart["tooltip"] = {"useHTML": True}
            # self.chart["tooltip"]["headerFormat"] = "{y} vs. {x}<br>".format(x=x, y=y)
            self.chart["tooltip"]["headerFormat"] = ""
            
            point_format = ["<b>{x}:</b> {{point.x}}<br><b>{y}:</b> {{point.y}}".format(x=self.arg_x, y=self.arg_y)]
            if self.arg_pid or self.arg_struct:
                point_format.append("<b>{pid}:</b> {{point.id}}".format(pid=self.arg_pid))
            self.chart["tooltip"]["pointFormat"] = "<br>".join(point_format)

            if self.arg_color_by:
                if self.dlen != len(d[color_by]):
                    raise ValueError("'{x}' and '{color_by}' must have the same length.".format(x=self.arg_x, color_by=self.arg_color_by))
                self.dcolor_by = list(d[color_by])
            
            if self.arg_z:
                if self.dlen != len(d[z]):
                    raise ValueError("'{x}' and '{z}' must have the same length.".format(x=self.arg_x, pid=self.arg_pid))
                self.dz = list(d[z])

            if not self.legend:
                self.chart["legend"] = {'enabled': False}
            else:
                self.chart["legend"] = {'enabled': True}
            self.chart["chart"] = {"type": "scatter", "zoomType": "xy"}
            
            if self.arg_color_by:
                if self.legend != False:
                    self.chart["legend"] = {'enabled': True}
                self.chart["tooltip"]["headerFormat"] = '<b>{color_by}: {{series.name}}</b><br>'.format(color_by=self.arg_color_by)
                series = self._data_series()
                self.chart["series"].extend(series)
            else:
                tmp_d = {"x": self.dx, "y": self.dy}
                if self.arg_struct:
                    tmp_d["id"] = [self._structure_tooltip(i) for i in range(self.dlen)]
                elif self.arg_pid:
                    tmp_d["id"] = self.dpid
                data = self._data_tuples(tmp_d)
                self.chart["series"].append({"name": "series", "data": data})
        
        if self.kind == "column":
            self.chart["chart"] = {"type": "column", "zoomType": "xy"}
            self._data_columns()

    
    def show(self, debug=False):
        formatter = string.Template(CHART_TEMPL)
        # if debug:
        #     print(self.chart)
        chart_json = json.dumps(self.chart)
        html = formatter.substitute({"id": self.chart_id, "chart": chart_json,
                                     "height": self.height})
        if debug:
            print(html)
        return HTML(html)
