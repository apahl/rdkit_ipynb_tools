#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# html_templates.py
"""
Created on Wed Jun 01 2015

@author: pahl
"""

import time
import os.path as op

TABLE_OPTIONS = {"cellspacing": "1", "cellpadding": "1", "border": "1", 
                 "align": "", "height": "60px", "summary": "Table", } # "width": "800px",

# PAGE_OPTIONS = {"icon": "icons/chart_bar.png", "css": ["css/style.css", "css/collapsible_list.css"], 
#                 "scripts": ["lib/jquery.js", "lib/highcharts.js", "script/folding.js"]}
PAGE_OPTIONS = {"icon": "icons/benzene.png", "css": ["css/style.css"]}
JSME = "lib/jsme/jsme.nocache.js"

HTML_FILE_NAME = "html/mol_table/index.html"


def tag(name, content, options=None, lf_open=False, lf_close=False):
    """creates a HTML stub with closed tags of type <name> around <content> 
    with additional <options::dict> in the opening tag
    when lf_(open|close)==True, the respective tag will be appended with a line feed.
    returns: html stub as list"""
    
    if lf_open:
        lf_open_str = "\n"
    else:
        lf_open_str = ""
    if lf_close:
        lf_close_str = "\n"
    else:
        lf_close_str = ""
    
    option_str = ""
    if options:
        option_list = [" "]
        for option in options:
            option_list.extend([option, '="', str(options[option]), '" '])
    
        option_str = "".join(option_list)
    
    stub = ["<{}{}>{}".format(name, option_str, lf_open_str)]
    if type(content) == list:
        stub.extend(content)
    else:
        stub.append(content)
    
    stub.append("</{}>{}".format(name, lf_close_str))
    
    return stub


def page(content, title="Results", options=PAGE_OPTIONS):
    """create a full HTML page from a list of stubs below
    options dict:
      css:     list of CSS style file paths to include.
      scripts: list of javascript library file paths to include.
      icon: path to icon image
    returns HTML page as STRING !!!"""
    
    # override the title if there is a title in <options>
    if "title" in options and len(options["title"]) > 2:
        title = options["title"]

    if "icon" in options and len(options["icon"]) > 2:
        icon_str = '<link rel="shortcut icon" href="{}" />'.format(options["icon"])
    else:
        icon_str = ""

    if "css" in options and options["css"]:
        css = options["css"]
        if type(css) != list:
            css = [css]
        
        css_str = "".join(['  <link rel="stylesheet" type="text/css" href="{}">\n'.format(file_name) for file_name in css])

    else:
        css_str = ""

    if "scripts" in options and options["scripts"]:
        scripts = options["scripts"]
        if type(scripts) != list:
            scripts = [scripts]

        js_str  = "".join(['  <script src="{}"></script>\n'.format(file_name) for file_name in scripts])
    
    else:
        js_str = ""
    
    if type(content) == list:
        content_str = "".join(content)
    else:
        content_str = content
    
    html_page = """<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <title>{title}</title>
  {icon_str}
{css_str}
{js_str}
</head>
<body>
{content_str}
</body>
</html>
""".format(title=title, icon_str=icon_str, css_str=css_str, js_str=js_str, content_str=content_str)
    
    return html_page


def write(text, fn=HTML_FILE_NAME):
    with open(fn, "w") as f:
        f.write(text)


def script(content):
    return tag("script", content, lf_open=True, lf_close=True)


def img(src, options=None):
    """takes a src, returns an img tag"""

    option_str = ""
    if options:
        option_list = [" "]
        for option in options:
            option_list.extend([option, '="', str(options[option]), '" '])
    
        option_str = "".join(option_list)
    
    stub = ['<img {}src="{}" alt="icon" />'.format(option_str, src)]
    
    return stub


def table(content, options=TABLE_OPTIONS):
    tbody = tag("tbody", content, lf_open=True, lf_close=True)
    return tag("table", tbody, options, lf_open=True, lf_close=True)


def tr(content, options=None):
    return tag("tr", content, options, lf_close=True)    


def td(content, options=None):
    return tag("td", content, options, lf_close=False)


def p(content):
    return tag("p", content, lf_open=True, lf_close=True)


def div(content, options=None):
    return tag("div", content, options, lf_close=False)


def ul(content):
    return tag("ul", content, lf_open=True, lf_close=True)


def li(content):
    return tag("li", content, lf_open=False, lf_close=True)


def li_lf(content): # list item with opening line feed
    return tag("li", content, lf_open=True, lf_close=True)
    

def b(content, options=None):
    return tag("b", content, options, lf_close=False)


def a(content, options):
    """the anchor tag requires an "href" in options,
    therefore the "options" parameter is not optional in this case (klingt komisch, ist aber so)"""
    return tag("a", content, options, lf_close=False)
