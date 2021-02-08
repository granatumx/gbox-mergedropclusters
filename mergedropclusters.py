#!/usr/bin/env python

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from colour import Color
from matplotlib.patches import Polygon
import statistics as st
import time
from random import randrange
import re
import sys

from granatum_sdk import Granatum

def parse(st):
    return list(map(lambda s: s.strip(), list(filter(lambda s: s != "", st.split(',')))))

def main():
    tic = time.perf_counter()

    gn = Granatum()
    assay = gn.pandas_from_assay(gn.get_import('assay'))
    groups = gn.get_import('groups')

    inv_map = {}
    for k, v in groups.items():
        inv_map[v] = inv_map.get(v, []) + [k]

    drop_set = parse(gn.get_arg('drop_set'))
    merge_set_1 = parse(gn.get_arg('merge_set_1'))
    merge_set_2 = parse(gn.get_arg('merge_set_2'))
    merge_set_3 = parse(gn.get_arg('merge_set_3'))
    relabel_set_1 = gn.get_arg('relabel_set_1')
    relabel_set_2 = gn.get_arg('relabel_set_2')
    relabel_set_3 = gn.get_arg('relabel_set_3')

    if len(merge_set_1) > 0:
        if relabel_set_1 == "":
            relabel_set_1 = " + ".join(merge_set_1)

    if len(merge_set_2) > 0:
        if relabel_set_2 == "":
            relabel_set_2 = " + ".join(merge_set_2)

    if len(merge_set_3) > 0:
        if relabel_set_3 == "":
            relabel_set_3 = " + ".join(merge_set_3)

    try:
        for ds in drop_set:
            cells = inv_map[ds]
            gn.add_result("Dropping {} cells that match {}".format(len(cells), ds), "markdown")
            assay = assay.drop(cells, axis=1)
            groups = {key:val for key, val in groups.items() if val != ds}
    except Exception as e:
        gn.add_result("Error found in drop set, remember it should be comma separated: {}".format(e), "markdown")
    
    try:
        if len(merge_set_1) > 0:
            merge_set_1_cells = []
            for ms1 in merge_set_1:
                merge_set_1_cells = merge_set_1_cells + inv_map[ms1]
            groups[merge_set_1_cells] = relabel_set_1
            for cell in merge_set_1_cells:
                groups[cell] = relabel_set_1

        if len(merge_set_2) > 0:
            merge_set_2_cells = []
            for ms2 in merge_set_2:
                merge_set_2_cells = merge_set_2_cells + inv_map[ms2]
            for cell in merge_set_2_cells:
                groups[cell] = relabel_set_2

        if len(merge_set_3) > 0:
            merge_set_3_cells = []
            for ms3 in merge_set_3:
                merge_set_3_cells = merge_set_3_cells + inv_map[ms3]
            for cell in merge_set_3_cells:
                groups[cell] = relabel_set_3
    except Exception as e:
        gn.add_result("Error found in merge sets, remember it should be comma separated: {}".format(e), "markdown")

    toc = time.perf_counter()
    time_passed = round(toc - tic, 2)

    gn.export_statically(gn.assay_from_pandas(assay), "Label adjusted assay")
    gn.export_statically(groups, "Adjusted labels")

    timing = "* Finished sample coloring step in {} seconds*".format(time_passed)
    gn.add_result(timing, "markdown")

    gn.commit()


if __name__ == "__main__":
    main()
