#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# nb_tools.py
"""
##############
Notebook Tools
##############

*Created on Wed Apr 22 16:37:35 2015 by A. Pahl*

A set of tools to use in the IPython (JuPyTer) Notebook
"""

# from rdkit.Chem import AllChem as Chem
# from rdkit.Chem import Draw
# Draw.DrawingOptions.atomLabelFontSize = 18
# import rdkit.Chem.Descriptors as Desc

# import sys
import time


def is_interactive_ipython():
    try:
        get_ipython()
        ipy = True
        print("> interactive IPython session.")
    except NameError:
        ipy = False
    return ipy


IPYTHON = is_interactive_ipython()

if IPYTHON:
    from IPython.core.display import HTML, Javascript, display
    import uuid


class ProgressbarJS():
    """A class to display a Javascript progressbar in the IPython notebook."""
    def __init__(self, color="#43ace8"):
        if IPYTHON:
            self.bar_id = str(uuid.uuid4())
            self.eta_id = str(uuid.uuid4())
            # possible colours: #94CAEF (blue from HTML reports),
            #                   #d6d2d0 (grey from window decorations)
            #
            self.pb = HTML(
                """
                <table style="border: none;"><tbody><tr style="border: none;">
                <td style="border: none;"><div style="border: 1px solid black; height:6px; width:500px">
                  <div id="{}" style="background-color:{}; height:4px; width:0%">&nbsp;</div>
                </div></td><td style="border: none;">&nbsp;ETA:&nbsp;</td><td style="border: none;"><div id="{}" width=100px></div></td>
                </tr></tbody></table>
                """.format(self.bar_id, color, self.eta_id))
            self.prev_time = 0.0
            self.start_time = time.time()
            display(self.pb)


    def update(self, perc, force=False):
        """update the progressbar
        in: progress in percent"""
        if IPYTHON:
            # make sure that the update function is not called too often:
            self.cur_time = time.time()
            if force or (self.cur_time - self.prev_time >= 0.25):
                if perc > 100: perc = 100
                if perc >= 25:
                    eta = (100 - perc) * (self.cur_time - self.start_time) / perc
                    eta_str = format_seconds(eta)
                else:
                    eta_str = "..."
                self.prev_time = self.cur_time
                display(Javascript("""
                                     $('div#{}').width('{}%');
                                     $('div#{}').text('{}');
                                   """.format(self.bar_id, perc, self.eta_id, eta_str)))

    def done(self):
        """finalize with a full progressbar for aesthetics"""
        if IPYTHON:
            display(Javascript("""
                                 $('div#{}').width('{}%');
                                 $('div#{}').text('{}');
                               """.format(self.bar_id, 100, self.eta_id, "done")))


def show_progress(iterable, iter_len=0):
    """A convenience wrapper for the ProgressBar class around iterables.

    Parameters:
        iterable (list, generator): The iterable object over which to loop.
        iter_len (int): Optional length of the object. This can be given if the iter object is a generator.

    Returns:
        A generator from iterable and displays the Javascript toolbar in the notebook."""

    if iter_len == 0:
        iter_len = len(iterable)

    steps = iter_len // 100
    if steps < 1:
        steps = 1
    pb = ProgressbarJS()

    for x, item in enumerate(iterable):
        if x % steps == 0:
            pb.update(100 * x / iter_len)

        yield item

    pb.done()


def format_seconds(seconds):
    seconds = int(seconds)
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    t_str = "{:02.0f}h {:02d}m {:02d}s".format(h, m, s)
    return t_str
