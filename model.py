#!/usr/bin/env python2
# -*- encoding: utf-8 -*-
from __future__ import print_function

import itertools
import collections

import json
import numpy

import newick

all_cards = """🂱 🂲 🂳 🂴 🂵 🂶 🂷 🂸 🂹 🂺
🂡 🂢 🂣 🂤 🂥 🂦 🂧 🂨 🂩 🂪
🃁 🃂 🃃 🃄 🃅 🃆 🃇 🃈 🃉 🃊
🃑 🃒 🃓 🃔 🃕 🃖 🃗 🃘 🃙 🃚""".split()

def draw_card(cards=all_cards):
    return cards[numpy.random.randint(len(cards))]


with open("real_forms.json") as lexicon:
    real_words = json.load(lexicon)

real_words = [w for w in real_words if w]

with open("real_forms.json", "w") as lexicon:
    json.dump(real_words, lexicon, indent=2)

def invent_random_word():
    return ''.join(draw_card(real_words))


def random_tree(depth=15, split_on="🂱🂡🃁🃑", root=None):
    if root is None:
        root = newick.Node('0')
    if depth > 0:
        root.length += 1
        depth -= 1
        if draw_card() in split_on:
            left = newick.Node(root.name+"l")
            root.add_descendant(left)
            random_tree(depth, split_on, root=left)
            right = newick.Node(root.name+"r")
            root.add_descendant(right)
            random_tree(depth, split_on, root=right)
        else:
            random_tree(depth, split_on, root=root)
    return root

def general_binary_model_on_tree(tree, gain_on="🂡🂢🂣🂤🂥🂦🂧🂨🂩🂪", loss_on="🂡🂢🂣🂤🂥🂦🂧🂨", at_root="🂩🂪"):
    data = {}
    if draw_card() in at_root:
        data[None] = True
    else:
        data[None] = False
    for node in tree.walk(mode="preorder"):
        state = data[node.ancestor.name]
        for i in range(int(node.length)):
            if state:
                if draw_card() in loss_on:
                    state = False
            else:
                if draw_card() in gain_on:
                    state = True
        data[node.name] = state
    return data


def random_observed_data(tree, probabilities=[("🂡🂢🂣🂤🂥", "🂡🂢🂣🂤", "🂡🂢🂣🂤🂥🂦🂧🂨🂩🂪")]*4):
    leaves = [n.name for n in tree.get_leaves()]
    data = {l: [] for l in leaves}
    for prob in probabilities:
            for language, value in general_binary_model_on_tree(tree, *prob).items():
                if language in leaves:
                    data[language].append(value)
    return data


def naive_monte_carlo(data, depth=15, log_bad_trees=True):
    recent_languages = len(data)
    data_as_counter = collections.Counter(
        tuple(x) for x in data.values())
    tree = None
    while True:
        suggested_tree = random_tree(depth, root=newick.Node('s'))
        if len(suggested_tree.get_leaves()) != len(data):
            if log_bad_trees:
                print("Wrong number of leaves: {:}".format(suggested_tree.newick))
                yield None
                continue
        else:
            tree = suggested_tree
        if tree is None:
            continue
        s_data = random_observed_data(tree)
        s_data_as_counter = collections.Counter(
            tuple(x) for x in s_data.values())
        if data_as_counter != s_data_as_counter:
            print("Not generating the right data: {:}".format(tree.newick))
            yield None
            continue
        used = set()
        for s_l in tree.get_leaves():
            for l, v in data.items():
                if l not in used and v == s_data[s_l.name]:
                    s_l.name = l
                    used.add(l)
                    break
            else:
                raise RuntimeError
        yield tree


def dollo_model_on_tree(tree, new_form="🂡🂢🂣🂤🂥", existing_forms=None):
    if existing_forms is None:
        existing_forms = set()
    data = {}
    form = max(existing_forms, default=-1) + 1
    existing_forms.add(form)
    for node in tree.walk(mode="preorder"):
        for i in range(int(node.length)):
            if draw_card() in new_form:
                form += 1
                existing_forms.add(form)
        data[node.name] = form
    return data


def random_observed_data(tree, models=[dollo_model_on_tree]*3):
    leaves = [n.name for n in tree.get_leaves()]
    data = {l: [] for l in leaves}
    for model in models:
        for language, value in model(tree).items():
            if language in leaves:
                data[language].append(value)
    return data


def data_pattern(data):
    """

    """
    data = list(data.items())
    data.sort(key=lambda x: x[1])
    data.sort(key=lambda row:
              [[r[i] for n, r in data].count(x)
               for i, x in enumerate(row[1])])
    for c in range(len(data[0][1])):
        values = {row[c]: i for i, (name, row) in enumerate(data)}
        for name, row in data:
            row[c] = values[row[c]]
    return tuple(zip(*data))

assert data_pattern({1: [1]})[1] == data_pattern({1: [2]})[1]
assert data_pattern({1: [1], 2: [2]})[1] != data_pattern({1: [2], 2: [2]})[1]
assert data_pattern({1: [3], 2: [2]})[1] == data_pattern({2: [2], 1: [1]})[1]
assert data_pattern({1: [1, 3], 2: [1, 2]})[1] == data_pattern({2: [4, 2], 1: [4, 1]})[1]
assert data_pattern({1: [1, 3], 2: [1, 2], 3: [2, 3], 4:[2, 2]})[1] == data_pattern({1: [1, 2], 2: [1, 1], 3: [2, 2], 4:[2, 1]})[1]
assert data_pattern({1: [1, 3], 2: [1, 2], 3: [2, 3], 4:[2, 2]})[1] != data_pattern({1: [1, 1], 2: [1, 1], 3: [2, 2], 4:[2, 2]})[1]
assert data_pattern({1: [1, 1], 2: [1, 2], 3: [2, 1], 4:[2, 2]})[1] == data_pattern({1: [1, 1], 2: [1, 2], 3: [2, 2], 4:[2, 1]})[1]

def naive_monte_carlo(data, tree_generator=random_tree, models=[dollo_model_on_tree]*3, verbose=True):
    recent_languages = len(data)
    names, classes = data_pattern(data)
    tree = None
    while True:
        tree = tree_generator()
        if tree is None:
            continue
        if len(tree.get_leaves()) != len(data):
            if verbose:
                print("Wrong number of leaves: {:}".format(tree.newick))
            yield None
            continue
        s_data = random_observed_data(tree, models)
        s_names, s_classes = data_pattern(s_data)
        if classes != s_classes:
            if verbose:
                print("Not generating the right data: {:}".format(tree.newick))
            yield None
            continue
        for old, new in zip(s_names, names):
            tree.get_node(old).name = new
        yield tree

try:
    with open("true.tree") as treefile:
        true_tree = newick.load(treefile)[0]
except (ValueError, FileNotFoundError):
    true_tree = random_tree()
    with open("true.tree", "w") as treefile:
            print(true_tree.newick, file=treefile, end=";\n", flush=True)

try:
    with open("data.json") as datafile:
        true_data = json.load(datafile)
except (json.JSONDecodeError, FileNotFoundError):
    true_data = random_observed_data(true_tree)
    with open("data.json", "w") as datafile:
        json.dump(true_data, datafile, indent=2)

if __name__ == "__main__":
    print(true_tree.newick)
    print(true_data)
    input()

    with open("fine_trees.nex", "w") as treeset:
        for i, tree in enumerate(naive_monte_carlo(true_data)):
            if tree:
                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                print(tree.newick)
                for node in tree.walk():
                    if node.name in true_data:
                        continue
                    node.name = None
                print(tree.newick, file=treeset, end=";\n", flush=True)
                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            if i > 8000:
                break
        for i, tree in enumerate(naive_monte_carlo(true_data, visit_bad_trees=False)):
            if tree:
                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                print(tree.newick)
                for node in tree.walk():
                    if node.name in true_data:
                        continue
                    node.name = None
                print(tree.newick, file=treeset, end=";\n", flush=True)
                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            if i > 8000:
                break
        for i, tree in enumerate(naive_monte_carlo(true_data, visit_bad_trees=False, verbose=False)):
            if tree:
                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                print(tree.newick)
                for node in tree.walk():
                    if node.name in true_data:
                        continue
                    node.name = None
                print(tree.newick, file=treeset, end=";\n", flush=True)
                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
