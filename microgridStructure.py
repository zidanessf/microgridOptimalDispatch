from microgrid_Model import *
'''Initialize a special case of microgrid'''
class MicrogridCase:
    def __init__(self):
        microgrid_device = dict()
        microgrid_device['PV_1'] = PV()
        microgrid_device['ES_1'] = electricStorage()
        microgrid_device['ABSC_1'] = absorptionChiller()
        microgrid_device['Boiler_1'] = boiler()
        microgrid_device['CS_1'] = coldStorage()
        microgrid_device['AC_1'] = airConditioner()
        microgrid_device['GT_1'] = gasTurbine()
        microgrid_device['ut'] = utility()
        microgrid_device['inv'] = inverter()
        self.device = microgrid_device
        self.NumOfTime = 96
    def getKey(self,type):
        return [key for key in self.device.keys() if isinstance(self.device[key], type)]
