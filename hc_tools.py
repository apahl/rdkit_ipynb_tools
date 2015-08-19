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

CHART_KINDS = ["scatter"]

print("- loading highcharts...")
display(HTML(HIGHCHARTS))


class Chart():
    def __init__(self, kind="scatter", **kwargs):
        if not kind in CHART_KINDS:
            raise ValueError("{} is not a supported chart kind ({})".format(kind, CHART_KINDS))

        self.kind = kind
        self.height= kwargs.get("height", 400)
        radius = kwargs.get("r", kwargs.get("radius", 5)) # accept "r" or "radius" for this option
        self.chart_id = time.strftime("%y%m%d%H%M%S")
        self.chart = {}
        self.chart["title"] = {"text": kwargs.get("title", "{} plot".format(self.kind))}
        self.chart["subtitle"] = {"text": kwargs.get("subtitle")}
        self.chart["series"] = []
        self.chart["plotOptions"] = {"scatter": {"marker": {"radius": radius}}}
        
    
    def _data_tuples(self, d, x, y, z, pid):
        """Generate the data tuples required for Highcharts scatter plot."""
        data = []
        dlen = len(d[x])
        for i in range(dlen):
            tmp_d = {"x": d[x][i], "y": d[y][i]}
            if z:
                tmp_d["z"] = d[z][i]
            if pid:
                tmp_d["id"] = d[pid][i]

            data.append(tmp_d)
        
        return data
    

    def _data_series(self, d, x, y, color_by, pid):
        # [{"name": "A", "data": [{"x": 1, "y": 2}, {"x": 2, "y": 3}]},
        #  {"name": "B", "data": [{"x": 2, "y": 3}, {"x": 3, "y": 4}]}]
        z = None  # not implemented yet
        series = []
        names = set(d[color_by])
        color_series_x = {name: [] for name in names}
        color_series_y = {name: [] for name in names}
        if pid:
            color_series_id = {name: [] for name in names}
        
        for i in range(len(d[x])):
            color_series_x[d[color_by][i]].append(d[x][i])
            color_series_y[d[color_by][i]].append(d[y][i])
            if pid:
                color_series_id[d[color_by][i]].append(d[pid][i])

        for name in names:
            tmp_d = {x: color_series_x[name], y: color_series_y[name]}
            if pid:
                tmp_d[pid] = color_series_id[name]
            series_dict = {"name": name}
            series_dict["data"] = self._data_tuples(tmp_d, x, y, z, pid)
            series.append(series_dict)
        
        return series
        


    def add_data(self, d, x="x", y="y", z="z", pid="id", **kwargs):
        """Add the data to the chart.
        d is the input dictionary, x, y [, and z] are the keys for the properties to plot.
        pid is the optional key to a (compound) id to be displayed in the tooltip."""
        if not x in d or not y in d:
            raise KeyError("'{x}' and '{y}' are required parameters for scatter plot, but could not all be found in dict.".format(x=x, y=y))
        
        if len(d[x]) != len(d[y]):
            raise ValueError("'{x}' and '{y}' must have the same length.".format(x=x, y=y))
            
        if pid in d:
            if len(d[x]) != len(d[pid]):
                raise ValueError("'{x}' and '{pid}' must have the same length.".format(x=x, pid=pid))
        else:
            pid = None

            
        color_by = kwargs.get("color_by")
        if color_by:
            if color_by in d:
                if len(d[x]) != len(d[pid]):
                    raise ValueError("'{x}' and '{color_by}' must have the same length.".format(x=x, color_by=color_by))
            else:
                raise KeyError("'{}' was not found in d".format(color_by))
        
        self.chart["credits"] = {'enabled': False}
        self.chart["xAxis"] = {"title": {"enabled": True, "text": x}}
        self.chart["yAxis"] = {"title": {"enabled": True, "text": y}}
        
        # defining the tooltip
        self.chart["tooltip"] = {}
        # self.chart["tooltip"]["headerFormat"] = "{y} vs. {x}<br>".format(x=x, y=y)
        self.chart["tooltip"]["headerFormat"] = ""
        
        point_format = ["<b>{x}:</b> {{point.x}}<br><b>{y}:</b> {{point.y}}".format(x=x, y=y)]
        if pid:
            point_format.append("<b>{pid}:</b>  {{point.id}}".format(pid=pid))
        self.chart["tooltip"]["pointFormat"] = "<br>".join(point_format)

        if self.kind == "scatter":
        
            if z in d:
                if len(d[x]) != len(d[z]):
                    raise ValueError("'{x}' and '{z}' must have the same length.".format(x=x, pid=pid))
            else:
                z = None

            self.chart["legend"] = {'enabled': False}
            self.chart["chart"] = {"type": "scatter", "zoomType": "xy"}
            
            if color_by:
                self.chart["tooltip"]["headerFormat"] = '<b>{color_by}: {{series.name}}</b><br>'.format(color_by=color_by)
                series = self._data_series(d, x, y, color_by, pid)
                self.chart["series"].extend(series)
            else:
                data = self._data_tuples(d, x, y, z, pid)
                self.chart["series"].append({"name": "series", "data": data})

    
    def show(self, debug=False):
        formatter = string.Template(CHART_TEMPL)
        chart_json = json.dumps(self.chart)
        html = formatter.substitute({"id": self.chart_id, "chart": chart_json,
                                     "height": self.height})
        if debug:
            print(html)
        return HTML(html)
