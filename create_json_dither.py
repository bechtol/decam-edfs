#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
import os
from collections import OrderedDict as odict
import copy
import datetime
import json

import numpy as np
import pylab as plt
import pandas as pd

# EDFS footprint
EDFS = np.rec.fromrecords([
    [63.25, -45.67],
    [65.35, -46.10],
    [66.40, -47.25],
    [65.99, -48.72],
    [59.25, -51.19],
    [56.95, -50.82],
    [55.90, -49.40],
    [56.80, -47.99],
    [63.25, -45.67],
],names=['ra','dec'])

PROGRAM = 'EDFS'
EXPTIME = 300
OVERHEAD = 30
BANDS = odict([
    ('u',12), # 10 is sufficent, but 12 is multiple of 2
    ('g',24),
    ('r',24),
    ('i',24),
    ('z',24),
])

# Baseline
#HEX = odict([
#    [1, (60.88, -47.65)],
#    [2, (61.70, -49.32)],
#    [3, (63.38, -46.75)],
#    [4, (64.20, -48.42)],
#    [5, (58.38, -48.55)],
#    [6, (59.20, -50.22)],
#])
# Compact
HEX = odict([
    [1, (60.88, -47.75)],
    [2, (61.70, -49.22)],
    [3, (62.88, -47.05)],
    [4, (63.70, -48.52)],
    [5, (58.88, -48.45)],
    [6, (59.70, -49.92)],
    [7, (64.88, -46.35)],
    [8, (65.70, -47.82)],
    [9, (56.88, -49.15)],
    [10, (57.70, -50.62)],
])

# http://www.ctio.noao.edu/noao/content/DECam-What
RA_GAP = (153 * 0.2637) / 3600. # deg
DEC_GAP = (201 * 0.2637) / 3600. # deg
#DITHER = [[0., 0.],
#          [1.2 * RA_GAP, 1.2 * DEC_GAP],
#          [-2.2 * RA_GAP, -2.2 * DEC_GAP]]
DITHER = [[0., 0.],
          [1.2 * RA_GAP, 1.2 * DEC_GAP],
          [-1.2 * RA_GAP, 2.2 * DEC_GAP]]

# Default configuration ('None' values must be filled)
SISPI_DICT = odict([
    ("object",  None),
    ("expTime", EXPTIME),
    ("RA",      None),
    ("dec",     None),
    ("filter",  None),
    #("count",   1),
    ("expType", "object"),
    ("program", PROGRAM),
    #("wait",    "False"),
    ("propid",  "2020B-0278"),
    ("proposer","Bechtol"),
    ("comment", ""),
])

def write_json(outfile,data,**kwargs):
    kwargs.setdefault('indent',4)
    json.encoder.FLOAT_REPR = lambda o: format(o, '.4f')

    with open(outfile,'w') as out:
        # It'd be nice to have a header
        #out.write(header())
        out.write(json.dumps(data,**kwargs))

def create_json():
    for ii,(ra,dec) in HEX.items():
        for band,nexp in BANDS.items():
            fields = []
            for jj in range(0, int(nexp//2)):
                dither_index = jj % len(DITHER)
                name = 'hex_%i'%(ii)
                d = dict(copy.deepcopy(SISPI_DICT))
                d['object'] = "EDFS field: "+name+" dither_%i"%(dither_index)
                d['RA'] = ra + (DITHER[dither_index][0] / np.cos(dec))
                d['dec'] = dec + DITHER[dither_index][1]
                d['filter'] = band
                #d['count'] = int(nexp//2)
                #d['count'] = 1
                d['seqtot'] = int(nexp//2)
                d['seqnum'] = jj + 1
                fields.append(d)
                
            filename = 'json_dither/edfs_%s_%s.json'%(name,band)
            print("  Writing %s..."%filename)
            write_json(filename,fields)

    #filename = 'json/edfs_all.json'
    #print("  Writing %s..."%filename)
    #write_json(filename,fields)

    #return pd.DataFrame(fields)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()

    create_json()
