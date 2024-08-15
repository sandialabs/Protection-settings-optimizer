if __package__ in [None, '']:
    from OpenDSS_funs import getSysInfo
    from OpenDSS_funs import simFaults
    from OpenDSS_funs import getLoadVI
    from OpenDSS_funs import getLineVI
    from OpenDSS_funs import getLoadInfo
    from OpenDSS_funs import getFaultInfo
    from OpenDSS_funs import getDeviceData
    from OpenDSS_Example_funs import run
else:
    from .OpenDSS_funs import getSysInfo
    from .OpenDSS_funs import simFaults
    from .OpenDSS_funs import getLoadVI
    from .OpenDSS_funs import getLineVI
    from .OpenDSS_funs import getLoadInfo
    from .OpenDSS_funs import getFaultInfo
    from .OpenDSS_funs import getDeviceData
    from .OpenDSS_Example_funs import run
