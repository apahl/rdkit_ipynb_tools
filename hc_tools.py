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
from IPython.html import widgets
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
    #TODO: fill this stub
    def __init__(self, kind="scatter", **kwargs):
        if not kind in CHART_KINDS:
            raise ValueError("{} is not a supported chart kind ({})".format(kind, CHART_KINDS))

        self.kind = kind
        self.height= kwargs.get("height", 400)
        self.chart = {}
        self.chart_id = time.strftime("%y%m%d%H%M%S")
        self.chart["title"] = {"text": kwargs.get("title", "{} plot".format(self.kind))}
        self.chart["subtitle"] = {"text": kwargs.get("subtitle")}
    
    
    def _data_tuples(self, d, x, y, z, pid):
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
    

    def add_data(self, d, x="x", y="y", z="z", pid="id", **kwargs):
        if not x in d or not y in d:
            raise ValueError("'{x}' and '{y}' are required parameters for scatter plot, but could not be found in dict.".format(x=x, y=y))
        
        if len(d[x]) != len(d[y]):
            raise ValueError("'{x}' and '{y}' must have the same length.".format(x=x, y=y))
            
        if pid in d:
            if len(d[x]) != len(d[pid]):
                raise ValueError("'{x}' and '{pid}' must have the same length.".format(x=x, pid=pid))
        else:
            pid = None
        
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
            
            data = self._data_tuples(d, x, y, z, pid)
            self.chart["series"] = []
            self.chart["series"].append({"data": data})

    
    def show(self, debug=False):
        formatter = string.Template(CHART_TEMPL)
        chart_json = json.dumps(self.chart)
        html = formatter.substitute({"id": self.chart_id, "chart": chart_json,
                                     "height": self.height})
        if debug:
            print(html)
        return HTML(html)
