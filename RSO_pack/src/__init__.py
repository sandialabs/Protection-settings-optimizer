if __package__ in [None, '']:
    # Misc Functions
    from GenNxGraph import index_dict
    from Read_CSV_Functions import pdef
    from Read_CSV_Functions import abc2012

    # Adapt functions
    from ADAPT_OPTV10 import runSettingsOptimizer
    from write2csv import writeSettingsToExcel
    from SEN import checkSensivity
    from OT_funs import OCIT,OCOT,OCTCC_Name,OCTCC_Num
    from Read_CSV_Functions import read_Relays_CSV_Data,read_Fault_CSV_Data

else:
    # Misc Functions
    from .GenNxGraph import index_dict
    from .Read_CSV_Functions import pdef
    from .Read_CSV_Functions import abc2012

    # Adapt functions
    from .ADAPT_OPTV10 import runSettingsOptimizer
    from .write2csv import writeSettingsToExcel
    from .SEN import checkSensivity
    from .OT_funs import OCIT,OCOT,OCTCC_Name,OCTCC_Num
    from .Read_CSV_Functions import read_Relays_CSV_Data,read_Fault_CSV_Data