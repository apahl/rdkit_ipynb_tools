#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
##############
File Templates
##############

*Created on Sun Apr 3 19:30 2016 by A. Pahl*

File templates for Reporting.
"""


STYLE = """<style>
  body{
  background-color: #FFFFFF;
  font-family: freesans, arial, verdana, sans-serif;
}
th {
  border-collapse: collapse;
  border-width: thin;
  border-style: solid;
  border-color: black;
  text-align: left;
  font-weight: bold;
}
td {
  border-collapse:collapse;
  border-width:thin;
  border-style:solid;
  border-color:black;
  padding: 5px;
}
table {
  border-collapse:collapse;
  border-width:thin;
  border-style:solid;
  border-color:black;
  background-color: #FFFFFF;
  text-align: left;
}
</style>"""

CLUSTER_JS_PY = """# list method index is not (yet) available
def index(l, elmnt):
    for idx, el in enumerate(l):
        if el == elmnt:
            return idx
    return -1  # not found

class Cluster:
    def __init__(self):
        self.current_cluster = 0
        self.list_idx = 0
        self.cluster_list = [${cluster_list}]

    def show_cluster(self):
        document.getElementById("display_no").innerHTML = self.current_cluster
        document.getElementById("cluster_frame").src = "clusters/cluster_{}.html".format(self.current_cluster)

    def get_cluster(self):
        cluster = int(document.getElementById("inp_cluster_no").value)
        idx = index(self.cluster_list, cluster)
        if idx >= 0:
            self.current_cluster = cluster
            self.list_idx = idx
            self.show_cluster()

    def next_cluster(self):
        if self.list_idx < len(self.cluster_list) - 1:
            self.list_idx += 1
            self.current_cluster = self.cluster_list[self.list_idx]
            self.show_cluster()

    def prev_cluster(self):
        if self.list_idx > 0:
            self.list_idx -= 1
            self.current_cluster = self.cluster_list[self.list_idx]
            self.show_cluster()

cluster = Cluster()
"""

CLUSTER_HTML = """<!DOCTYPE html>
<html>
<head>
${style}
</head>
<body>
<div style="position:fixed; left:0; width:100%; top:0; height:100%;">
    <iframe name="cluster" id="cluster_frame" src="empty.html" width=80% height=100% align=right>
        no iframes spported.
    </iframe>
    <h3>Cluster Report</h3>
    <p>
        Cluster No.:<br>
        <input type="text" id="inp_cluster_no"><button onclick="report_clusters.cluster.get_cluster()">show</button>
    </p><p>
        <button onclick="report_clusters.cluster.prev_cluster()">Previous</button>
        <span id="display_no" style="margin: 30px;" align="center">&nbsp;</span>
        <button onclick="report_clusters.cluster.next_cluster()">Next</button>
    </p>
</div>
<script src="__javascript__/cluster_js.min.js"></script>
</body>
</html>
"""
