if __package__ in [None, '']:
    from RONM_format import getRONMSysInfo
    from RONM_format import getRONMFaultData
    from RONM_format import get_Sw_Status
    from RONM_format import getRONMDeviceData
else:
    from .RONM_format import getRONMSysInfo
    from .RONM_format import getRONMFaultData
    from .RONM_format import get_Sw_Status
    from .RONM_format import getRONMDeviceData