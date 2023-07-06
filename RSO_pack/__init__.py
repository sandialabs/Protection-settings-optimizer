# -*- coding: utf-8 -*-
"""
load functions
"""

import os
import pygad
import networkx as nx
from copy import deepcopy

# Misc Functions
from .src.GenNxGraph import index_dict
from .src.Read_CSV_Functions import pdef
from .src.Read_CSV_Functions import abc2012

# Adapt functions
from .src.ADAPT_OPTV10 import runSettingsOptimizer
from .src.write2csv import writeSettingsToExcel
from .src.SEN import checkSensivity
from .src.OT_funs import OCIT,OCOT,OCTCC_Name,OCTCC_Num
from .src.Read_CSV_Functions import read_Relays_CSV_Data,read_Fault_CSV_Data

# OpenDSS Functions 
from .src_OpenDSS.OpenDSS_funs import getSysInfo
from .src_OpenDSS.OpenDSS_funs import simFaults
#from .src_OpenDSS.OpenDSS_funs import getBESSInfo
#from .src_OpenDSS.OpenDSS_funs import getGeneratorInfo
#from .src_OpenDSS.OpenDSS_funs import getPvInfo
#from .src_OpenDSS.OpenDSS_funs import getBusInfo
#from .src_OpenDSS.OpenDSS_funs import getTransformerInfo
#from .src_OpenDSS.OpenDSS_funs import getLineInfo
#from .src_OpenDSS.OpenDSS_funs import getFuseInfo
#from .src_OpenDSS.OpenDSS_funs import getRecloserInfo
#from .src_OpenDSS.OpenDSS_funs import getRelayInfo
from .src_OpenDSS.OpenDSS_funs import getLoadVI
from .src_OpenDSS.OpenDSS_funs import getLineVI
from .src_OpenDSS.OpenDSS_funs import getLoadInfo
from .src_OpenDSS.OpenDSS_funs import getFaultInfo
from .src_OpenDSS.OpenDSS_funs import getDeviceData

# RONM Functions
from .src_RONM.RONM_format import getRONMSysInfo
from .src_RONM.RONM_format import getRONMFaultData
from .src_RONM.RONM_format import get_Sw_Status
from .src_RONM.RONM_format import getRONMDeviceData