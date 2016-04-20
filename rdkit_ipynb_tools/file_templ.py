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

CLUSTER_REPORT_INTRO = """<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <link rel="stylesheet" type="text/css" href="css/style.css" />
  <link rel="stylesheet" type="text/css" href="css/collapsible_list.css" />
  <link rel="shortcut icon" href="icons/chart_bar.png" />

  <script src="lib/jquery.min.js"></script>
  <script src="lib/folding.js"></script>

  <title>Clustering</title>
</head>
<body>
<p><button onclick="expand_all()">expand all clusters</button>
<button onclick="collapse_all()">collapse all clusters</button></p>
"""

CLUSTER_REPORT_EXTRO = """
<p><button onclick="expand_all()">expand all clusters</button>
<button onclick="collapse_all()">collapse all clusters</button></p>

  <script>
    folding();
  </script>

</body>
</html>"""

# string template ${cluster_list}:
CLUSTER_PY = """# list method index is not (yet) available
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

# string template ${cluster_list}:
CLUSTER_JS = """'use strict';function report_clusters(){function D(a,b,c){"undefined"==typeof b&&(b=a,a=0);"undefined"==typeof c&&(c=1);if(0<c&&a>=b||0>c&&a<=b)return[];for(var e=[];0<c?a<b:a>b;a+=c)e.push(a);return e}function E(a){return F(D(n(a)),a)}function z(a){if(null==a||"object"==typeof a)return a;var b={},c;for(c in obj)a.hasOwnProperty(c)&&(b[c]=a[c]);return b}function G(a){if(null==a||"object"==typeof a)return a;var b={},c;for(c in obj)a.hasOwnProperty(c)&&(b[c]=G(a[c]));return b}function m(a){return a?
[].slice.apply(a):[]}function h(a){a=a?[].slice.apply(a):[];a.__class__=h;return a}function k(a){var b=[];if(a)for(var c=0;c<a.length;c++)b.add(a[c]);b.__class__=k;return b}function H(){var a=[],b;for(b in this)r(b)||a.push(b);return a}function I(){var a=[],b;for(b in this)r(b)||a.push([b,this[b]]);return a}function J(a){delete this[a]}function t(a){if(!a||a instanceof Array){var b={};if(a)for(var c=0;c<a.length;c++){var e=a[c];b[e[0]]=e[1]}}else b=a;Object.defineProperty(b,"__class__",{value:t,enumerable:!1,
writable:!0});Object.defineProperty(b,"py_keys",{value:H,enumerable:!1});Object.defineProperty(b,"py_items",{value:I,enumerable:!1});Object.defineProperty(b,"py_del",{value:J,enumerable:!1});return b}function p(a){try{return a.__str__()}catch(b){return new String(a)}}var d={},u=function(a,b,c){if(""!=b){b=b.split(".");for(var e=b.length,g=0;g<b.length;g++){if(!a.hasOwnProperty(b[g])){e=g;break}a=a[b[g]]}for(g=e;g<b.length;g++)a[b[g]]={},a=a[b[g]]}for(var d in c)a[d]=c[d]};d.__nest__=u;var A=function(a){a.__inited__||
a.__all__.__init__(a.__all__);return a.__all__};d.__init__=A;var f=function(a,b,c){return a&&(a.hasOwnProperty("__class__")||"string"==typeof a||a instanceof String)?(c&&Object.defineProperty(a,c,{value:function(){var c=[].slice.apply(arguments);return b.apply(null,[a].concat(c))},writable:!0,enumerable:!0,configurable:!0}),function(){var c=[].slice.apply(arguments);return b.apply(null,[a].concat(c))}):b};d.__get__=f;var v=function(a,b,c){for(var e=function(){var a=[].slice.apply(arguments);return e.__new__(a)},
g=b.length-1;0<=g;g--){var d=b[g],q;for(q in d){var B=Object.getOwnPropertyDescriptor(d,q);Object.defineProperty(e,q,B)}}e.__name__=a;e.__bases__=b;for(q in c)B=Object.getOwnPropertyDescriptor(c,q),Object.defineProperty(e,q,B);return e};d.__class__=v;var w=d.__class__("object",[],{__init__:function(a){},__name__:"object",__bases__:[],__new__:function(a){var b=Object.create(this,{__class__:{value:this,enumerable:!0}});this.__init__.apply(null,[b].concat(a));return b}});d.object=w;d.__pragma__=function(){};
u(d,"org.transcrypt.__base__",{__all__:{__inited__:!1,__init__:function(a){var b=v("__Envir__",[w],{get __init__(){return f(this,function(a){a.transpiler_name="transcrypt";a.transpiler_version="3.5.138";a.target_subdir="__javascript__"})}}),c=b();a.__Envir__=b;a.__envir__=c}}});u(d,"org.transcrypt.__standard__",{__all__:{__inited__:!1,__init__:function(a){var b=v("Exception",[w],{get __init__(){return f(this,function(a){var b=h([].slice.apply(arguments).slice(1));a.args=b})},get __repr__(){return f(this,
function(a){return n(a.args)?"{}{}".format(a.__class__.__name__,K(h(a.args))):"???"})},get __str__(){return f(this,function(a){return 1<n(a.args)?p(h(a.args)):n(a.args)?p(a.args[0]):"???"})}}),c=v("ValueError",[b],{}),e=function(a,b,c){if("undefined"==typeof b||null!=b&&b.__class__==l)b=null;if("undefined"==typeof c||null!=c&&c.__class__==l)c=!1;if(arguments.length){var e=arguments.length-1;if(arguments[e]&&arguments[e].__class__==l){var e=arguments[e--],d;for(d in e)switch(d){case "iterable":a=e[d];
break;case "key":b=e[d];break;case "reverse":c=e[d]}}}b?a.sort(function(a,c){if(arguments.length){var e=arguments.length-1;if(arguments[e]&&arguments[e].__class__==l){var e=arguments[e--],d;for(d in e)switch(d){case "a":a=e[d];break;case "b":c=e[d]}}}return b(a)>b(c)}):a.sort();c&&a.reverse()};a.Exception=b;a.ValueError=c;a.__sort__=e;a.sorted=function(a,b,c){if("undefined"==typeof b||null!=b&&b.__class__==l)b=null;if("undefined"==typeof c||null!=c&&c.__class__==l)c=!1;if(arguments.length){var d=
arguments.length-1;if(arguments[d]&&arguments[d].__class__==l){var d=arguments[d--],f;for(f in d)switch(f){case "iterable":a=d[f];break;case "key":b=d[f];break;case "reverse":c=d[f]}}}f=C(a)==t?z(a.py_keys()):z(a);e(f,b,c);return f}}}});u(d,"",A(d.org.transcrypt.__base__));var L=d.__envir__;u(d,"",A(d.org.transcrypt.__standard__));var O=d.__sort__;L.executor_name=L.transpiler_name;d.main={__file__:""};d.__except__=null;var l=function(a){a.__class__=l;a.constructor=Object;return a};d.___kwargdict__=
l;d.property=function(a,b){b||(b=function(){});return{get:function(){return a(this)},set:function(a){b(this,a)},enumerable:!0}};d.__merge__=function(a,b){var c={},e;for(e in a)c[e]=a[e];for(e in b)c[e]=b[e];return c};var M=function(){for(var a=[].slice.apply(arguments),b="",c=0;c<a.length;c++)b+=p(a[c])+" ";console.log(b)};d.print=M;console.log.apply=function(){M([].slice.apply(arguments).slice(1))};d.__in__=function(a,b){return C(b)==t?-1<b.py_keys().indexOf(a):-1<b.indexOf(a)};var r=function(a){return a.startswith("__")&&
a.endswith("__")||"constructor"==a||a.startswith("py_")};d.__specialattrib__=r;var n=function(a){try{return a.length}catch(c){var b=0;for(attrib in a)r(attrib)||b++;return b}};d.len=n;var N={__name__:"bool"};d.bool=N;var x=function(a){if(isNaN(a))throw"ValueError";return+a};x.__name__="float";d.float=x;var y=function(a){return x(a)|0};y.__name__="int";d.int=y;var C=function(a){try{return a.__class__}catch(c){var b=typeof a;return"boolean"==b?N:"number"==b?0==a%1?y:x:b}};d.type=C;d.isinstance=function(a,
b){function c(a){if(a==b)return!0;for(var d=0;d<a.__bases__.length;d++)if(c(a.__bases__[d],b))return!0;return!1}return c(a.__class__)};var K=function(a){try{return a.__repr__()}catch(d){try{return a.__str__()}catch(f){try{if(a.constructor==Object){var b="{",c=!1,e;for(e in a)if(!r(e)){var g=e.isnumeric()?e:"'"+e+"'";c?b+=", ":c=!0;try{b+=g+": "+a[e].__repr__()}catch(h){b+=g+": "+a[e].toString()}}return b+"}"}return"boolean"==typeof a?a.toString().capitalize():a.toString()}catch(h){return console.log("ERROR: Could not evaluate repr (<object of type "+
typeof a+">)"),"???"}}}};d.repr=K;d.chr=function(a){return String.fromCharCode(a)};d.org=function(a){return a.charCodeAt(0)};var F=function(){var a=[].slice.call(arguments);return(0==a.length?[]:a.reduce(function(a,c){return a.length<c.length?a:c})).map(function(b,c){return a.map(function(a){return a[c]})})};d.zip=F;d.range=D;d.enumerate=E;d.copy=z;d.deepcopy=G;d.list=m;Array.prototype.__class__=m;m.__name__="list";Array.prototype.__getslice__=function(a,b,c){0>a&&(a=this.length+a);null==b?b=this.length:
0>b&&(b=this.length+b);for(var e=m([]);a<b;a+=c)e.push(this[a]);return e};Array.prototype.__setslice__=function(a,b,c,e){0>a&&(a=this.length+a);null==b?b=this.length:0>b&&(b=this.length+b);if(null==c)Array.prototype.splice.apply(this,[a,b-a].concat(e));else for(var d=0;a<b;a+=c)this[a]=e[d++]};Array.prototype.__repr__=function(){if(this.__class__==k&&!this.length)return"set()";for(var a=this.__class__&&this.__class__!=m?this.__class__==h?"(":"{":"[",b=0;b<this.length;b++){b&&(a+=", ");try{a+=this[b].__repr__()}catch(c){a+=
this[b].toString()}}this.__class__==h&&1==this.length&&(a+=",");return a+=this.__class__&&this.__class__!=m?this.__class__==h?")":"}":"]"};Array.prototype.__str__=Array.prototype.__repr__;Array.prototype.append=function(a){this.push(a)};Array.prototype.clear=function(){this.length=0};Array.prototype.extend=function(a){this.push.apply(this,a)};Array.prototype.insert=function(a,b){this.splice(a,0,b)};Array.prototype.remove=function(a){a=this.indexOf(a);if(-1==a)throw"KeyError";this.splice(a,1)};Array.prototype.py_pop=
function(a){return void 0==a?this.pop():this.splice(a,1)[0]};Array.prototype.py_sort=function(){O.apply(null,[this].concat([].slice.apply(arguments)))};d.tuple=h;h.__name__="tuple";d.set=k;k.__name__="set";Array.prototype.__bindexOf__=function(a){a+="";for(var b=0,c=this.length-1;b<=c;){var e=(b+c)/2|0,d=this[e]+"";if(d<a)b=e+1;else if(d>a)c=e-1;else return e}return-1};Array.prototype.add=function(a){-1==this.indexOf(a)&&this.push(a)};Array.prototype.discard=function(a){a=this.indexOf(a);-1!=a&&this.splice(a,
1)};Array.prototype.isdisjoint=function(a){this.sort();for(var b=0;b<a.length;b++)if(-1!=this.__bindexOf__(a[b]))return!1;return!0};Array.prototype.issuperset=function(a){this.sort();for(var b=0;b<a.length;b++)if(-1==this.__bindexOf__(a[b]))return!1;return!0};Array.prototype.issubset=function(a){return k(a.slice()).issuperset(this)};Array.prototype.union=function(a){for(var b=k(this.slice().sort()),c=0;c<a.length;c++)-1==b.__bindexOf__(a[c])&&b.push(a[c]);return b};Array.prototype.intersection=function(a){this.sort();
for(var b=k(),c=0;c<a.length;c++)-1!=this.__bindexOf__(a[c])&&b.push(a[c]);return b};Array.prototype.difference=function(a){a=k(a.slice().sort());for(var b=k(),c=0;c<this.length;c++)-1==a.__bindexOf__(this[c])&&b.push(this[c]);return b};Array.prototype.symmetric_difference=function(a){return this.union(a).difference(this.intersection(a))};Array.prototype.update=function(){var a=[].concat.apply(this.slice(),arguments).sort();this.clear();for(var b=0;b<a.length;b++)a[b]!=a[b-1]&&this.push(a[b])};d.__keys__=
H;d.__items__=I;d.__del__=J;d.dict=t;t.__name__="dict";d.str=p;String.prototype.__class__=p;p.__name__="str";String.prototype.__repr__=function(){return(-1==this.indexOf("'")?"'"+this+"'":'"'+this+'"').replace("\n","\\n")};String.prototype.__str__=function(){return this};String.prototype.capitalize=function(){return this.charAt(0).toUpperCase()+this.slice(1)};String.prototype.endswith=function(a){return""==a||this.slice(-a.length)==a};String.prototype.find=function(a,b){return this.indexOf(a,b)};
Object.defineProperty(String.prototype,"format",{get:function(){return f(this,function(a){var b=h([].slice.apply(arguments).slice(1)),c=0;return a.replace(/\{(\w*)\}/g,function(a,d){""==d&&(d=c++);if(d==+d)return"undefined"==b[d]?a:b[d];for(var f=0;f<b.length;f++)if("object"==typeof b[f]&&"undefined"!=typeof b[f][d])return b[f][d];return a})})},enumerable:!0});String.prototype.isnumeric=function(){return!isNaN(parseFloat(this))&&isFinite(this)};String.prototype.join=function(a){return a.join(this)};
String.prototype.lower=function(){return this.toLowerCase()};String.prototype.py_replace=function(a,b,c){return this.split(a,c).join(b)};String.prototype.lstrip=function(){return this.replace(/^\s*/g,"")};String.prototype.rfind=function(a,b){return this.lastIndexOf(a,b)};String.prototype.rsplit=function(a,b){var c=this.split(a||/s+/);return b?[c.slice(0,-b).join(a)].concat(c.slice(-b)):c};String.prototype.rstrip=function(){return this.replace(/\s*$/g,"")};String.prototype.py_split=function(a,b){a||
(a=" ");return this.split(a||/s+/,b)};String.prototype.startswith=function(a){return 0==this.indexOf(a)};String.prototype.strip=function(){return this.trim()};String.prototype.upper=function(){return this.toUpperCase()};d.__neg__=function(a){return"object"==typeof a&&"__neg__"in a?a.__neg__():-a};d.__matmul__=function(a,b){return a.__matmul__(b)};d.__mul__=function(a,b){return"object"==typeof a&&"__mul__"in a?a.__mul__(b):"object"==typeof b&&"__rmul__"in b?b.__rmul__(a):a*b};d.__div__=function(a,
b){return"object"==typeof a&&"__div__"in a?a.__div__(b):"object"==typeof b&&"__rdiv__"in b?b.__rdiv__(a):a/b};d.__add__=function(a,b){return"object"==typeof a&&"__add__"in a?a.__add__(b):"object"==typeof b&&"__radd__"in b?b.__radd__(a):a+b};d.__sub__=function(a,b){return"object"==typeof a&&"__sub__"in a?a.__sub__(b):"object"==typeof b&&"__rsub__"in b?b.__rsub__(a):a-b};d.__getitem__=function(a,b){return"object"==typeof a&&"__getitem__"in a?a.__getitem__(b):a[b]};d.__setitem__=function(a,b,c){"object"==
typeof a&&"__setitem__"in a?a.__setitem__(b,c):a[b]=c};d.__getslice__=function(a,b,c,d){return"object"==typeof a&&"__getitem__"in a?a.__getitem__([b,c,d]):a.__getslice__(b,c,d)};d.__setslice__=function(a,b,c,d,f){"object"==typeof a&&"__setitem__"in a?a.__setitem__([b,c,d],f):a.__setslice__(b,c,d,f)};d.__call__=function(){var a=[].slice.apply(arguments);return"object"==typeof a[0]&&"__call__"in a[0]?a[0].__call__.apply(null,a.slice(1)):a[0].apply(null,a.slice(1))};(function(){var a=function(a,b){for(var c=
E(a),d=0;d<c.length;d++){var f=c[d],h=f[0];if(f[1]==b)return h}return-1},b=v("Cluster",[w],{get __init__(){return f(this,function(a){a.current_cluster=0;a.list_idx=0;a.cluster_list=m([${cluster_list}])})},get show_cluster(){return f(this,function(a){document.getElementById("display_no").innerHTML=a.current_cluster;document.getElementById("cluster_frame").src="clusters/cluster_{}.html".format(a.current_cluster)})},get get_cluster(){return f(this,function(b){var c=y(document.getElementById("inp_cluster_no").value),
d=a(b.cluster_list,c);0<=d&&(b.current_cluster=c,b.list_idx=d,b.show_cluster())})},get next_cluster(){return f(this,function(a){a.list_idx<n(a.cluster_list)-1&&(a.list_idx++,a.current_cluster=a.cluster_list[a.list_idx],a.show_cluster())})},get prev_cluster(){return f(this,function(a){0<a.list_idx&&(a.list_idx--,a.current_cluster=a.cluster_list[a.list_idx],a.show_cluster())})}}),c=b();d.Cluster=b;d.cluster=c;d.index=a})();return d}window.report_clusters=report_clusters();
"""
