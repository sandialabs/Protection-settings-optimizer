# -*- coding: utf-8 -*-
"""
load functions
"""

import os
import pygad
import networkx as nx
from copy import deepcopy

if __package__ in [None, '']:
    import src
    import src_OpenDSS
    import src_RONM
else:
    from .src import *
    from .src_OpenDSS import *
    from .src_RONM import *