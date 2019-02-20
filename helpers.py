#!/usr/bin/env python2
# -*- encoding: utf-8 -*-

from __future__ import print_function

import numpy
from IPython.display import Image, HTML, display
from cgi import escape

try:
    import ete3
    def draw(tree):
        return ete3.Tree(tree.newick+";", format=3).render("%%inline")
except ImportError:
    def draw(tree):
        print(tree.ascii_art())

from model import invent_random_word, naive_monte_carlo, random_observed_data

die_faces = "⚀⚁⚂⚃⚄⚅"

def roll_die():
    return die_faces[numpy.random.randint(len(die_faces))]

def big(text):
    display(HTML('&nbsp<p style="font-size: 400%">{:}</p>&nbsp'.format(escape(text))))

