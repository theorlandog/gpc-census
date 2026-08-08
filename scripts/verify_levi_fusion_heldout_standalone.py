#!/usr/bin/env python3
"""Independent standard-library verifier for the held-out Levi audit.

The exact primitive tau rows and unique wider-system signatures are embedded
below. The verifier recomputes every scored row-level quantity and compares the
result with ``results/data/levi_fusion_heldout_results.json``.

This file imports neither gpc_census nor third-party packages.
"""

from __future__ import annotations

import json
import math
import pathlib
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product

ROOT = pathlib.Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "data" / "levi_fusion_heldout_results.json"
INPUTS = json.loads(r'''
{
 "(3,10)": {
  "d": 10,
  "label_counts": {
   "PADDING": 52,
   "PRIMITIVE": 41
  },
  "n": 3,
  "rows": [
   {
    "label": "PADDING",
    "tau": [
     2,
     2,
     -4,
     2,
     -4,
     -4,
     -1,
     -1,
     -1,
     -1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     2,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     2
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     1,
     -1,
     -1,
     -1,
     0,
     0,
     0,
     0,
     -1,
     1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     1,
     0,
     -1,
     -1,
     -1,
     -1,
     0,
     0,
     0,
     1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     2,
     1,
     -3,
     -2,
     -1,
     -1,
     -1,
     -1,
     0,
     1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     2,
     -2,
     -1,
     -1,
     -1,
     -1,
     0,
     1,
     -3,
     1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     -2,
     1,
     1,
     -2,
     1,
     -2,
     -2,
     -2,
     -2
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     1,
     -2,
     1,
     -2,
     -2,
     1,
     -2,
     -2,
     -2
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     1,
     -2,
     -2,
     1,
     1,
     -2,
     -2,
     -2,
     -2
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     -2,
     1,
     1,
     1,
     1,
     -2,
     -2,
     -2,
     -2,
     -2
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     1,
     -1,
     0,
     1,
     -2,
     0,
     -1,
     -1,
     -1,
     0
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     1,
     -1,
     0,
     0,
     -1,
     1,
     -2,
     -1,
     -1,
     0
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     0,
     -1,
     1,
     1,
     -1,
     0,
     -2,
     -1,
     -1,
     0
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     0,
     -1,
     -1,
     -1,
     0,
     0,
     1,
     -2,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     -1,
     -1,
     0,
     0,
     1,
     -2,
     0,
     -1,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     -1,
     -1,
     0,
     0,
     0,
     -1,
     1,
     -2,
     -1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     1,
     0,
     -1,
     -1,
     0,
     1,
     -2,
     -1,
     -1,
     0
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     0,
     1,
     -1,
     -1,
     1,
     0,
     -2,
     -1,
     -1,
     0
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     -1,
     1,
     0,
     0,
     1,
     -1,
     -2,
     -1,
     -1,
     0
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     -1,
     0,
     1,
     1,
     0,
     -1,
     -2,
     -1,
     -1,
     0
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     -1,
     1,
     0,
     1,
     0,
     -2,
     -1,
     -1,
     -1,
     0
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     1,
     0,
     -1,
     1,
     -2,
     -1,
     0,
     -1,
     -1,
     0
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     0,
     1,
     -1,
     1,
     -1,
     -2,
     0,
     -1,
     -1,
     0
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     0,
     -1,
     1,
     -2,
     -1,
     -1,
     0,
     0,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     0,
     1,
     -1,
     1,
     -1,
     -2,
     -1,
     0,
     0,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     0,
     1,
     -1,
     1,
     -2,
     -1,
     0,
     -1,
     0,
     -1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     1,
     1,
     -2,
     0,
     -1,
     -1,
     0,
     -1,
     -1,
     0
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     1,
     -2,
     0,
     -1,
     -1,
     -1,
     0,
     0,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     -1,
     1,
     0,
     0,
     1,
     -2,
     -1,
     -1,
     0,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     -1,
     0,
     1,
     1,
     0,
     -2,
     -1,
     -1,
     0,
     -1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     1,
     1,
     -2,
     -1,
     0,
     0,
     -1,
     -1,
     -1,
     0
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     1,
     -2,
     -1,
     -1,
     0,
     0,
     0,
     -1,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     4,
     -5,
     -2,
     1,
     -5,
     -2,
     1,
     -5,
     -2
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     -2,
     1,
     4,
     -5,
     -5,
     -2,
     1,
     -5,
     -2
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     -2,
     1,
     1,
     4,
     -5,
     -2,
     -5,
     1,
     -5,
     -2
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     4,
     -5,
     1,
     -5,
     -2,
     -2,
     1,
     -5,
     -2
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     -1,
     2,
     -1,
     2,
     -1,
     -4,
     -1,
     -1,
     -1,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     -1,
     2,
     -1,
     -1,
     2,
     -1,
     -4,
     -1,
     -1,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     2,
     -1,
     -1,
     2,
     -4,
     -1,
     -1,
     -1,
     -1,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     2,
     -1,
     -1,
     -1,
     -1,
     2,
     -4,
     -1,
     -1,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     -1,
     -1,
     2,
     2,
     -1,
     -1,
     -4,
     -1,
     -1,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     2,
     2,
     -4,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     2,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     2,
     -4,
     -1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     1,
     -1,
     -1,
     -1,
     0,
     0,
     0,
     1,
     -2,
     0
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     1,
     1,
     -2,
     -1,
     -1,
     -1,
     0,
     0,
     0,
     0
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     0,
     -1,
     0,
     -1,
     0,
     -1,
     0,
     -1,
     0
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     0,
     1,
     -1,
     0,
     0,
     -1,
     -1,
     0,
     -1,
     0
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     0,
     0,
     0,
     1,
     -1,
     -1,
     -1,
     0,
     -1,
     0
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     1,
     4,
     -5,
     -2,
     1,
     -5,
     -2,
     -2,
     -2,
     1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     1,
     -2,
     1,
     4,
     -5,
     -5,
     -2,
     -2,
     -2,
     1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     -2,
     1,
     1,
     4,
     -5,
     -2,
     -5,
     -2,
     -2,
     1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     2,
     -3,
     -1,
     -2,
     -1,
     0,
     1,
     -3,
     -2
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     -1,
     0,
     1,
     2,
     -3,
     -2,
     -1,
     1,
     -3,
     -2
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     2,
     -2,
     0,
     1,
     -3,
     1,
     -3,
     -2,
     -1,
     -1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     1,
     -2,
     1,
     2,
     -3,
     0,
     -3,
     -2,
     -1,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     2,
     1,
     -3,
     -2,
     -1,
     -1,
     0,
     1,
     -3,
     -2
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     2,
     -2,
     -1,
     -1,
     0,
     1,
     -3,
     1,
     -3,
     -2
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     2,
     1,
     -3,
     -2,
     0,
     1,
     -3,
     -2,
     -1,
     -1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     1,
     2,
     -3,
     -2,
     1,
     0,
     -3,
     -2,
     -1,
     -1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     -2,
     2,
     0,
     1,
     1,
     -3,
     -3,
     -2,
     -1,
     -1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     -2,
     1,
     1,
     2,
     0,
     -3,
     -3,
     -2,
     -1,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     2,
     -3,
     1,
     -3,
     -2,
     -1,
     0,
     -1,
     -2
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     1,
     2,
     -3,
     1,
     -3,
     -2,
     0,
     -1,
     -2,
     -1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     2,
     1,
     -3,
     1,
     -3,
     -2,
     0,
     -2,
     -1,
     -1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     1,
     1,
     -2,
     2,
     -3,
     -3,
     0,
     -2,
     -1,
     -1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     1,
     2,
     -3,
     1,
     -2,
     -3,
     0,
     -2,
     -1,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     2,
     1,
     -3,
     1,
     -3,
     -2,
     -1,
     -1,
     0,
     -2
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     1,
     -2,
     2,
     -3,
     -3,
     -1,
     -1,
     0,
     -2
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     2,
     -3,
     1,
     -2,
     -3,
     -1,
     -1,
     0,
     -2
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     3,
     2,
     -5,
     -3,
     -2,
     -1,
     -2,
     -1,
     0,
     1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     3,
     -3,
     -2,
     -1,
     -2,
     -1,
     0,
     2,
     -5,
     1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     5,
     -4,
     -1,
     -4,
     -1,
     2,
     -7,
     -4,
     -1,
     2
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     5,
     -4,
     -1,
     2,
     -7,
     -4,
     -1,
     -4,
     -1,
     2
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     5,
     2,
     -7,
     -4,
     -1,
     -4,
     -1,
     -4,
     -1,
     2
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     -1,
     -1,
     2,
     5,
     -7,
     -4,
     -4,
     -4,
     -1,
     2
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     5,
     -4,
     -1,
     -4,
     -1,
     -4,
     -1,
     2,
     -7,
     2
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     2,
     5,
     -7,
     -1,
     -4,
     -4,
     -1,
     2,
     -7,
     -4
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     -1,
     -1,
     2,
     5,
     -7,
     -4,
     -4,
     2,
     -7,
     -4
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     5,
     -4,
     -1,
     2,
     -7,
     2,
     -7,
     -4,
     -1,
     -4
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     5,
     -4,
     -1,
     2,
     -7,
     -4,
     -1,
     2,
     -7,
     -4
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     2,
     -4,
     2,
     5,
     -7,
     -1,
     -7,
     -4,
     -1,
     -4
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     5,
     2,
     -7,
     -4,
     -1,
     -4,
     -1,
     2,
     -7,
     -4
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     5,
     -4,
     -1,
     -4,
     -1,
     2,
     -7,
     2,
     -7,
     -4
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     5,
     2,
     -7,
     -4,
     -1,
     2,
     -7,
     -4,
     -1,
     -4
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     2,
     5,
     -7,
     -4,
     2,
     -1,
     -7,
     -4,
     -1,
     -4
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     -4,
     5,
     -1,
     2,
     2,
     -7,
     -7,
     -4,
     -1,
     -4
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     -4,
     2,
     2,
     5,
     -1,
     -7,
     -7,
     -4,
     -1,
     -4
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     2,
     5,
     -7,
     2,
     -7,
     -4,
     -1,
     -1,
     -4,
     -4
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     5,
     2,
     -7,
     2,
     -7,
     -4,
     -1,
     -4,
     -1,
     -4
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     2,
     2,
     -4,
     5,
     -7,
     -7,
     -1,
     -4,
     -1,
     -4
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     2,
     5,
     -7,
     2,
     -4,
     -7,
     -1,
     -4,
     -1,
     -4
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     3,
     -3,
     -2,
     -1,
     -2,
     -1,
     0,
     1,
     -4,
     2
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     3,
     1,
     -4,
     -3,
     -2,
     -1,
     -2,
     -1,
     0,
     2
    ]
   }
  ]
 },
 "(3,9)": {
  "d": 9,
  "label_counts": {
   "PADDING": 31,
   "PRIMITIVE": 21
  },
  "n": 3,
  "rows": [
   {
    "label": "PRIMITIVE",
    "tau": [
     2,
     2,
     -4,
     2,
     -4,
     -4,
     -1,
     -1,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     0,
     -1,
     0,
     -1,
     0,
     -1,
     0,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     0,
     1,
     -1,
     0,
     0,
     -1,
     -1,
     0,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     0,
     0,
     0,
     1,
     -1,
     -1,
     -1,
     0,
     -1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     1,
     1,
     -2,
     0,
     -1,
     -1,
     -1,
     0,
     0
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     1,
     -2,
     -1,
     -1,
     0,
     0,
     0,
     -1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     1,
     0,
     -1,
     1,
     -2,
     -1,
     -1,
     0,
     0
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     0,
     -1,
     -1,
     -1,
     0,
     0,
     1,
     -2
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     -1,
     -1,
     0,
     0,
     0,
     -1,
     1,
     -2
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     -1,
     -1,
     0,
     0,
     1,
     -2,
     0,
     -1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     0,
     1,
     -1,
     1,
     -1,
     -2,
     -1,
     0,
     0
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     0,
     1,
     -1,
     1,
     -2,
     -1,
     0,
     -1,
     0
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     -1,
     1,
     0,
     0,
     1,
     -2,
     -1,
     -1,
     0
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     -1,
     0,
     1,
     1,
     0,
     -2,
     -1,
     -1,
     0
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     1,
     -2,
     1,
     -2,
     -2,
     1,
     -2,
     -2
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     1,
     -2,
     -2,
     1,
     1,
     -2,
     -2,
     -2
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     -2,
     1,
     1,
     -2,
     1,
     -2,
     -2,
     -2
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     -2,
     1,
     1,
     1,
     1,
     -2,
     -2,
     -2,
     -2
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     2,
     1,
     -3,
     1,
     -3,
     -2,
     -1,
     -1,
     0
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     2,
     1,
     -3,
     -2,
     -1,
     -1,
     0,
     1,
     -3
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     2,
     -2,
     -1,
     -1,
     0,
     1,
     -3,
     1,
     -3
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     1,
     2,
     -3,
     1,
     -3,
     -2,
     -1,
     0,
     -1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     1,
     2,
     -3,
     1,
     -2,
     -3,
     -1,
     -1,
     0
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     2,
     -3,
     -1,
     -2,
     -1,
     0,
     1,
     -3
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     1,
     1,
     -2,
     2,
     -3,
     -3,
     -1,
     -1,
     0
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     -1,
     0,
     1,
     2,
     -3,
     -2,
     -1,
     1,
     -3
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     2,
     2,
     -4,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     2,
     -1,
     -1,
     2,
     -4,
     -1,
     -1,
     -1,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     2,
     -1,
     -1,
     -1,
     -1,
     2,
     -4,
     -1,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     2,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     2,
     -4
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     -1,
     2,
     -1,
     2,
     -1,
     -4,
     -1,
     -1,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     -1,
     2,
     -1,
     -1,
     2,
     -1,
     -4,
     -1,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     -1,
     -1,
     2,
     2,
     -1,
     -1,
     -4,
     -1,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     4,
     -5,
     1,
     -5,
     -2,
     -2,
     1,
     -5
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     4,
     -5,
     -2,
     1,
     -5,
     -2,
     1,
     -5
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     1,
     -2,
     1,
     4,
     -5,
     -5,
     -2,
     1,
     -5
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     -2,
     1,
     1,
     4,
     -5,
     -2,
     -5,
     1,
     -5
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     5,
     2,
     -7,
     2,
     -7,
     -4,
     -1,
     -4,
     -1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     5,
     2,
     -7,
     -4,
     -1,
     2,
     -7,
     -4,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     5,
     -4,
     -1,
     2,
     -7,
     -4,
     -1,
     2,
     -7
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     5,
     -4,
     -1,
     2,
     -7,
     2,
     -7,
     -4,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     5,
     -4,
     -1,
     -4,
     -1,
     2,
     -7,
     2,
     -7
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     2,
     5,
     -7,
     2,
     -4,
     -7,
     -1,
     -4,
     -1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     2,
     5,
     -7,
     2,
     -7,
     -4,
     -1,
     -1,
     -4
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     2,
     5,
     -7,
     -1,
     -4,
     -4,
     -1,
     2,
     -7
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     2,
     5,
     -7,
     -4,
     2,
     -1,
     -7,
     -4,
     -1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     2,
     -4,
     2,
     5,
     -7,
     -1,
     -7,
     -4,
     -1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     2,
     2,
     -4,
     5,
     -7,
     -7,
     -1,
     -4,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     -1,
     -1,
     2,
     5,
     -7,
     -4,
     -4,
     2,
     -7
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     -4,
     5,
     -1,
     2,
     2,
     -7,
     -7,
     -4,
     -1
    ]
   },
   {
    "label": "PRIMITIVE",
    "tau": [
     -4,
     2,
     2,
     5,
     -1,
     -7,
     -7,
     -4,
     -1
    ]
   },
   {
    "label": "PADDING",
    "tau": [
     5,
     2,
     -7,
     -4,
     -1,
     -4,
     -1,
     2,
     -7
    ]
   }
  ]
 },
 "(4,10)": {
  "d": 10,
  "label_counts": {
   "FROZEN-CORE": 20,
   "PADDING": 59,
   "PRIMITIVE": 45,
   "STRUCTURAL": 1
  },
  "n": 4,
  "signatures": [
   {
    "labels": [
     "PADDING"
    ],
    "tau": [
     1,
     1,
     1,
     -3,
     1,
     -3,
     -3,
     -3,
     1,
     -3
    ]
   },
   {
    "labels": [
     "PRIMITIVE"
    ],
    "tau": [
     1,
     1,
     5,
     -7,
     5,
     -11,
     -7,
     -11,
     -3,
     -3
    ]
   },
   {
    "labels": [
     "PRIMITIVE"
    ],
    "tau": [
     0,
     1,
     1,
     -2,
     2,
     -4,
     -3,
     -3,
     -1,
     -1
    ]
   },
   {
    "labels": [
     "PADDING",
     "PRIMITIVE"
    ],
    "tau": [
     7,
     -1,
     3,
     -9,
     -5,
     -1,
     -9,
     -5,
     -1,
     -9
    ]
   },
   {
    "labels": [
     "PADDING"
    ],
    "tau": [
     1,
     -1,
     1,
     -1,
     1,
     -1,
     -3,
     -1,
     -1,
     -1
    ]
   },
   {
    "labels": [
     "PADDING"
    ],
    "tau": [
     1,
     0,
     0,
     -1,
     0,
     -1,
     -1,
     0,
     -1,
     -1
    ]
   },
   {
    "labels": [
     "PADDING"
    ],
    "tau": [
     3,
     -1,
     -1,
     -1,
     3,
     -5,
     -5,
     -5,
     -1,
     -5
    ]
   },
   {
    "labels": [
     "FROZEN-CORE",
     "PRIMITIVE"
    ],
    "tau": [
     5,
     1,
     1,
     -7,
     1,
     -7,
     -7,
     -3,
     -3,
     -3
    ]
   },
   {
    "labels": [
     "PRIMITIVE"
    ],
    "tau": [
     0,
     1,
     2,
     -3,
     3,
     -6,
     -4,
     -5,
     -2,
     -1
    ]
   },
   {
    "labels": [
     "PADDING"
    ],
    "tau": [
     1,
     -1,
     0,
     0,
     1,
     -2,
     -2,
     -1,
     0,
     -2
    ]
   },
   {
    "labels": [
     "FROZEN-CORE",
     "PADDING",
     "PRIMITIVE"
    ],
    "tau": [
     1,
     -1,
     3,
     -3,
     -1,
     1,
     -5,
     -3,
     -3,
     -1
    ]
   },
   {
    "labels": [
     "FROZEN-CORE",
     "PADDING",
     "PRIMITIVE"
    ],
    "tau": [
     0,
     0,
     1,
     -1,
     3,
     -4,
     -4,
     -3,
     -2,
     -2
    ]
   },
   {
    "labels": [
     "FROZEN-CORE",
     "PADDING",
     "PRIMITIVE"
    ],
    "tau": [
     1,
     2,
     0,
     -3,
     0,
     -3,
     -1,
     -2,
     -2,
     -1
    ]
   },
   {
    "labels": [
     "PRIMITIVE"
    ],
    "tau": [
     1,
     5,
     9,
     -15,
     1,
     -7,
     -15,
     -11,
     -7,
     -3
    ]
   }
  ]
 },
 "(4,9)": {
  "d": 9,
  "label_counts": {
   "FROZEN-CORE": 21,
   "PADDING": 14,
   "PRIMITIVE": 24,
   "STRUCTURAL": 1
  },
  "n": 4,
  "signatures": [
   {
    "labels": [
     "FROZEN-CORE",
     "PRIMITIVE"
    ],
    "tau": [
     3,
     1,
     1,
     -5,
     -3,
     -3,
     -1,
     -1,
     -1
    ]
   },
   {
    "labels": [
     "FROZEN-CORE"
    ],
    "tau": [
     1,
     1,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     1
    ]
   },
   {
    "labels": [
     "PADDING"
    ],
    "tau": [
     1,
     0,
     0,
     -1,
     0,
     -1,
     -1,
     0,
     -1
    ]
   },
   {
    "labels": [
     "PRIMITIVE"
    ],
    "tau": [
     1,
     1,
     1,
     -3,
     1,
     -3,
     -3,
     -3,
     1
    ]
   },
   {
    "labels": [
     "FROZEN-CORE",
     "PRIMITIVE"
    ],
    "tau": [
     7,
     -1,
     3,
     -9,
     -1,
     -9,
     -5,
     -5,
     -1
    ]
   },
   {
    "labels": [
     "FROZEN-CORE",
     "PRIMITIVE"
    ],
    "tau": [
     2,
     1,
     -2,
     -1,
     0,
     -3,
     -2,
     -1,
     0
    ]
   },
   {
    "labels": [
     "PADDING"
    ],
    "tau": [
     1,
     1,
     1,
     -3,
     -1,
     -1,
     -1,
     -1,
     -1
    ]
   },
   {
    "labels": [
     "FROZEN-CORE"
    ],
    "tau": [
     3,
     3,
     -1,
     -5,
     -1,
     -5,
     -1,
     -5,
     -1
    ]
   },
   {
    "labels": [
     "FROZEN-CORE"
    ],
    "tau": [
     3,
     1,
     -1,
     -3,
     -3,
     -3,
     -1,
     -1,
     1
    ]
   },
   {
    "labels": [
     "PRIMITIVE"
    ],
    "tau": [
     1,
     -1,
     0,
     0,
     1,
     -2,
     -2,
     -1,
     0
    ]
   },
   {
    "labels": [
     "PRIMITIVE"
    ],
    "tau": [
     2,
     0,
     1,
     -3,
     0,
     -3,
     -2,
     -1,
     -1
    ]
   },
   {
    "labels": [
     "FROZEN-CORE"
    ],
    "tau": [
     3,
     1,
     0,
     -4,
     -3,
     -2,
     -2,
     -1,
     0
    ]
   }
  ]
 },
 "(5,10)": {
  "d": 10,
  "label_counts": {
   "FROZEN-CORE": 45,
   "PADDING": 59,
   "PRIMITIVE": 56,
   "STRUCTURAL": 1
  },
  "n": 5,
  "signatures": [
   {
    "labels": [
     "PADDING"
    ],
    "tau": [
     -4,
     1,
     1,
     1,
     1,
     1,
     1,
     -4,
     -4,
     -4
    ]
   },
   {
    "labels": [
     "FROZEN-CORE"
    ],
    "tau": [
     2,
     2,
     2,
     -3,
     -3,
     -3,
     -3,
     -3,
     -3,
     2
    ]
   },
   {
    "labels": [
     "FROZEN-CORE",
     "PADDING"
    ],
    "tau": [
     -2,
     3,
     3,
     3,
     -7,
     -2,
     -2,
     -2,
     -7,
     -7
    ]
   },
   {
    "labels": [
     "PRIMITIVE"
    ],
    "tau": [
     4,
     4,
     9,
     -11,
     -6,
     -11,
     -6,
     -1,
     -16,
     -1
    ]
   },
   {
    "labels": [
     "FROZEN-CORE",
     "PADDING",
     "PRIMITIVE"
    ],
    "tau": [
     -2,
     3,
     -2,
     8,
     -7,
     8,
     -7,
     -17,
     -17,
     -12
    ]
   },
   {
    "labels": [
     "PRIMITIVE"
    ],
    "tau": [
     6,
     11,
     16,
     -4,
     -29,
     -24,
     -19,
     -14,
     -9,
     1
    ]
   },
   {
    "labels": [
     "FROZEN-CORE",
     "PADDING",
     "PRIMITIVE"
    ],
    "tau": [
     1,
     6,
     -4,
     6,
     -9,
     -4,
     1,
     -14,
     -9,
     -4
    ]
   },
   {
    "labels": [
     "FROZEN-CORE",
     "PADDING",
     "PRIMITIVE"
    ],
    "tau": [
     -2,
     8,
     -7,
     -2,
     -2,
     3,
     3,
     -7,
     -12,
     -12
    ]
   },
   {
    "labels": [
     "PADDING"
    ],
    "tau": [
     -3,
     2,
     2,
     2,
     -3,
     2,
     -3,
     -3,
     -3,
     -8
    ]
   },
   {
    "labels": [
     "FROZEN-CORE",
     "PADDING"
    ],
    "tau": [
     0,
     1,
     0,
     0,
     -1,
     0,
     -1,
     -1,
     0,
     -1
    ]
   },
   {
    "labels": [
     "PADDING",
     "PRIMITIVE"
    ],
    "tau": [
     -6,
     4,
     -1,
     4,
     -1,
     9,
     -6,
     -16,
     -16,
     -11
    ]
   },
   {
    "labels": [
     "FROZEN-CORE",
     "PRIMITIVE"
    ],
    "tau": [
     8,
     3,
     8,
     -12,
     -7,
     -7,
     -12,
     -2,
     -17,
     -2
    ]
   },
   {
    "labels": [
     "FROZEN-CORE",
     "PRIMITIVE"
    ],
    "tau": [
     -2,
     3,
     8,
     13,
     -22,
     -12,
     -7,
     -2,
     -17,
     -12
    ]
   },
   {
    "labels": [
     "PADDING",
     "PRIMITIVE"
    ],
    "tau": [
     2,
     7,
     -8,
     -3,
     2,
     12,
     -23,
     -18,
     -13,
     -8
    ]
   },
   {
    "labels": [
     "FROZEN-CORE"
    ],
    "tau": [
     6,
     6,
     1,
     -4,
     -9,
     -9,
     -9,
     -4,
     -4,
     1
    ]
   },
   {
    "labels": [
     "PADDING"
    ],
    "tau": [
     -6,
     4,
     -1,
     -1,
     -1,
     4,
     4,
     -6,
     -11,
     -11
    ]
   },
   {
    "labels": [
     "PRIMITIVE"
    ],
    "tau": [
     -2,
     4,
     -3,
     -1,
     0,
     1,
     2,
     -6,
     -5,
     -4
    ]
   },
   {
    "labels": [
     "PRIMITIVE"
    ],
    "tau": [
     6,
     11,
     16,
     -24,
     -19,
     -14,
     -9,
     1,
     -34,
     -4
    ]
   },
   {
    "labels": [
     "PRIMITIVE"
    ],
    "tau": [
     -8,
     -3,
     2,
     2,
     7,
     12,
     -8,
     -18,
     -18,
     -13
    ]
   },
   {
    "labels": [
     "PRIMITIVE"
    ],
    "tau": [
     4,
     9,
     9,
     -1,
     -21,
     -16,
     -11,
     -11,
     -6,
     -1
    ]
   }
  ]
 }
}
''')


def descriptor(tau: tuple[int, ...], particle_number: int) -> dict[str, object]:
    counts = Counter(tau)
    levels = tuple(sorted(counts))
    multiplicities = tuple(counts[level] for level in levels)
    allocations: list[tuple[tuple[int, ...], int]] = []
    for occupancy in product(*[range(multiplicity + 1) for multiplicity in multiplicities]):
        if sum(occupancy) != particle_number:
            continue
        if sum(level * count for level, count in zip(levels, occupancy)) != 0:
            continue
        dimension = math.prod(
            math.comb(multiplicity, count)
            for multiplicity, count in zip(multiplicities, occupancy)
        )
        allocations.append((tuple(occupancy), dimension))
    if not allocations:
        raise AssertionError("empty zero-grade module")
    return {
        "shape": tuple(sorted(tau)),
        "distinct_levels": len(levels),
        "fusion_rank": len(allocations),
        "zero_grade_dimension": sum(dimension for _occupancy, dimension in allocations),
    }


def full_system(payload: dict[str, object]) -> dict[str, object]:
    n = int(payload["n"])
    d = int(payload["d"])
    rows = payload["rows"]
    records = []
    for row in rows:
        item = descriptor(tuple(row["tau"]), n)
        item["label"] = row["label"]
        records.append(item)

    labels = Counter(record["label"] for record in records)
    by_label: dict[str, list[dict[str, object]]] = defaultdict(list)
    for record in records:
        by_label[record["label"]].append(record)

    primitive_shapes = {
        record["shape"] for record in records if record["label"] == "PRIMITIVE"
    }
    inherited_shapes = {
        record["shape"]
        for record in records
        if record["label"] in {"PADDING", "FROZEN-CORE"}
    }
    max_zdim = max(int(record["zero_grade_dimension"]) for record in records)
    return {
        "rows": len(records),
        "nonstructural_rows": len(records) - labels.get("STRUCTURAL", 0),
        "label_counts": dict(labels),
        "max_fusion_rank": max(int(record["fusion_rank"]) for record in records),
        "max_distinct_levels": max(int(record["distinct_levels"]) for record in records),
        "distinct_unoriented_signatures": len({record["shape"] for record in records}),
        "max_zero_grade_dimension": max_zdim,
        "ambient_dimension": math.comb(d, n),
        "max_zero_grade_density": list(
            (lambda value: (value.numerator, value.denominator))(
                Fraction(max_zdim, math.comb(d, n))
            )
        ),
        "mean_fusion_by_label": {
            label: list(
                (lambda value: (value.numerator, value.denominator))(
                    Fraction(
                        sum(int(record["fusion_rank"]) for record in selected),
                        len(selected),
                    )
                )
            )
            for label, selected in by_label.items()
        },
        "primitive_inherited_signature_overlap": len(
            primitive_shapes & inherited_shapes
        ),
    }


def compact_system(payload: dict[str, object]) -> dict[str, object]:
    n = int(payload["n"])
    d = int(payload["d"])
    records = []
    for row in payload["signatures"]:
        item = descriptor(tuple(row["tau"]), n)
        item["labels"] = set(row["labels"])
        records.append(item)
    nonstructural = [
        record for record in records if "STRUCTURAL" not in record["labels"]
    ]
    primitive = {
        record["shape"] for record in records if "PRIMITIVE" in record["labels"]
    }
    inherited = {
        record["shape"]
        for record in records
        if record["labels"] & {"PADDING", "FROZEN-CORE"}
    }
    max_zdim = max(int(record["zero_grade_dimension"]) for record in nonstructural)
    label_counts = {key: int(value) for key, value in payload["label_counts"].items()}
    return {
        "rows": sum(label_counts.values()),
        "nonstructural_rows": sum(
            value for key, value in label_counts.items() if key != "STRUCTURAL"
        ),
        "label_counts": label_counts,
        "max_fusion_rank": max(int(record["fusion_rank"]) for record in nonstructural),
        "max_distinct_levels": max(
            int(record["distinct_levels"]) for record in nonstructural
        ),
        "distinct_unoriented_signatures": len(nonstructural),
        "max_zero_grade_dimension": max_zdim,
        "ambient_dimension": math.comb(d, n),
        "max_zero_grade_density": list(
            (lambda value: (value.numerator, value.denominator))(
                Fraction(max_zdim, math.comb(d, n))
            )
        ),
        "primitive_inherited_signature_overlap": len(primitive & inherited),
    }


def main() -> int:
    expected = json.loads(RESULT.read_text())
    computed = {}
    for system, payload in INPUTS.items():
        computed[system] = (
            full_system(payload) if "rows" in payload else compact_system(payload)
        )

    comparison_keys = (
        "rows",
        "nonstructural_rows",
        "label_counts",
        "max_fusion_rank",
        "max_distinct_levels",
        "distinct_unoriented_signatures",
        "max_zero_grade_dimension",
        "ambient_dimension",
        "max_zero_grade_density",
        "primitive_inherited_signature_overlap",
    )
    for system, actual in computed.items():
        recorded = expected["systems"][system]
        for key in comparison_keys:
            if actual[key] != recorded[key]:
                raise AssertionError(
                    f"{system} {key}: {actual[key]!r} != {recorded[key]!r}"
                )
        if "mean_fusion_by_label" in actual:
            if actual["mean_fusion_by_label"] != recorded["mean_fusion_by_label"]:
                raise AssertionError(f"{system} mean fusion mismatch")

    score = expected["scorecard"]
    assert score["P1_distinct_tau_levels_at_most_2N"]["verdict"] == "FAIL"
    assert score["P1_distinct_tau_levels_at_most_2N"]["violating_systems"] == [
        "(3,10)",
        "(4,10)",
    ]
    for name in (
        "P2_zero_grade_mechanisms_compress_rows_by_at_least_2x",
        "P3_N3_primitive_mean_fusion_exceeds_inherited",
        "P4_primitive_and_inherited_mechanisms_overlap",
        "P5_N3_max_fusion_rank_at_most_d_minus_3",
        "P6_nonstructural_zero_grade_density_at_most_55_percent",
    ):
        assert score[name]["verdict"] == "PASS", name

    print("held-out Levi audit: PASS")
    print("P1: FAIL at (3,10), (4,10)")
    print("P2-P6: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
