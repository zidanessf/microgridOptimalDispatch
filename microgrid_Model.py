class PV:
    def __init__(self, om = 0.0005 , output = list()):
        self.om = om
        self.output = output
        self.result = list()
    def show(self):
        return {"om":self.om}

class DRHeatLoad:
    def __init__(self,upper_bound=1000,lower_bound=-1000):
        self.upper_bound = upper_bound
        self.lower_bound = lower_bound

class electricStorage:
    def __init__(self, om = 0.005, Cbw = 0.00075 , capacity = 13000, SOCmin = 0.1, SOCmax = 0.9, SOCint = 0.1, Pmax_in = 1250, Pmax_out = 1250, efficiency = 1, selfRelease = 0.0025):
        self.om = om
        self.Cbw = Cbw
        self.capacity = capacity
        self.SOCmin = SOCmin
        self.SOCmax = SOCmax
        self.SOCint = SOCint
        self.Pmax_in = Pmax_in
        self.Pmax_out = Pmax_out
        self.efficiency = efficiency
        self.selfRelease = selfRelease
        self.maxDetP = self.Pmax_out * 0.5
        self.power_into = {}
        self.power_outof = {}
        self.energy = {}
        self.SOCnow = SOCint
    def show(self):
        return {"om":self.om,
                "Cbw":self.Cbw,
                "capacity":self.capacity,
                "SOCmin":self.SOCmin,
                "SOCmax":self.SOCmax,
                "SOCint":self.SOCint,
                "Pmax_out":self.Pmax_out,
                "Pmax_in":self.Pmax_in,
                "selfRelease":self.selfRelease,
                "efficiency":self.efficiency}


class absorptionChiller:
    def __init__(self, om = 0.00008, COP_htc = 0.8, COP_hth = 1, Hmin = 0, Hmax = 1000, ElecCost = 0.002):
        self.om = om
        self.COP_htc = COP_htc
        self.COP_hth = COP_hth
        self.Hmin = Hmin
        self.Hmax = Hmax
        self.ElecCost = ElecCost
        self.result = {}
    def show(self):
        return {"om":self.om,
                "COP_htc":self.COP_htc,
                "COP_hth":self.COP_hth}

class boiler:
    def __init__(self, om = 0.04, Pmax = 1000, Pmin = 50, efficiency = 0.85,maxDetP=500,ON_OFF_COST=100):
        self.om = om
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.efficiency = efficiency
        self.result = {}
        self.ON_OFF_COST = ON_OFF_COST
        self.maxDetP=maxDetP
    def show(self):
        return {"om":self.om,
                "Pmax":self.Pmax,
                "Pmin":self.Pmin,
                "efficiency":self.efficiency}

class heatStorage:
    def __init__(self, om = 0.04, capacity = 2000, Tmin = 0, Tmax = 0.95, Tint = 0.1, Hmax_in = 1500, Hmax_out = 1500, selfRelease = 0.003):
        self.om = om
        self.capacity = capacity
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.Tint = Tint
        self.Hmax_in = Hmax_in
        self.Hmax_out = Hmax_out
        self.selfRelease = selfRelease
        self.maxDetP = self.Hmax_in * 0.4
        self.result = {}
    def show(self):
        return {"om":self.om,
                "capacity":self.capacity,
                "Tmax":self.Tmax,
                "Tint":self.Tint,
                "Hmax_out":self.Hmax_out,
                "Hmax_in":self.Hmax_in}

class coldStorage:
    def __init__(self, om = 0.01, capacity = 3000, Tmin = 0.1, Tmax = 0.95, Tint = 0.1, Hin = 500, Hout = 500, Pmin = 0, Pmax = 500, maxDetP = 500, EER = 3 , efficiency = 1 , COP = 3 , Partition_in = 0.9, Partition_out = 0.9, mode = "并联" ,selfRelease = 0.002):
        self.om = om
        self.capacity = capacity
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.Tint = Tint
        self.Hin = Hin
        self.Hout = Hout
        self.Pmax = Pmax
        self.Pmin = Pmin
        self.EER = EER
        self.efficiency = efficiency
        self.COP = COP
        self.Partition_in = Partition_in
        self.Partition_out = Partition_out
        self.mode = mode
        self.maxDetP = maxDetP
        self.selfRelease = selfRelease
        self.result_electricity_tank = {}
        self.result_cold_tank = {}
        self.result_electricity_ref = {}
    def show(self):
        return {"om":self.om,
                "capacity":self.capacity,
                "Tmax":self.Tmax,
                "Tmin":self.Tmin,
                "Tint":self.Tint,
                "Hin":self.Hin,
                "Hout":self.Hout,
                "COP":self.COP,
                "Partition_out":self.Partition_out,
                "Partition_in":self.Partition_in,
                "EER":self.EER,
                "efficiency":self.efficiency}

class airConditioner:## P>=0 cooling, P<0 heating
    def __init__(self, om = 0.0097, Pmax = 500, Pmin = 0, EER = 4.3 , COP = 3.6, maxDetP = 100):
        self.om = om
        self.Pmax = Pmax
        self.EER = EER
        self.COP = COP
        self.Pmin = Pmin
        self.maxDetP = maxDetP
        self.result = {}
    def show(self):
        return {"om":self.om,
                "Pmax":self.Pmax,
                "EER":self.EER,
                "COP":self.COP}

class gasTurbine:
    def __init__(self, om = 0.063, Pmax = 1000, Pmin = 50, efficiency = 0.33, heat_recycle = 0.6,maxDetP=200,ON_OFF_COST=200,Cost=0):
        self.om = om
        self.Pmax = Pmax
        self.Pmin = Pmin
        self.efficiency = efficiency
        self.heat_recycle = heat_recycle
        self.low_heat_recycle = 1 - heat_recycle
        self.HER = (1 - efficiency)/efficiency
        self.maxDetP = maxDetP
        self.ON_OFF_COST = ON_OFF_COST
        self.Cost = Cost
        self.result = list()
    def show(self):
        return {"om":self.om,
                "Pmax":self.Pmax,
                "Pmin":self.Pmin,
                "efficiency":self.efficiency,
                "HER":self.HER,
                "heat_recycle":self.heat_recycle}


class utility:
    def __init__(self, buy_price = 0.8, sell_price = 0, gas_price = 0.349,steam_price = 348/996, CO2_utility = 0.997, CO2_user = 0.18):
        self.buy_price = buy_price
        self.sell_price = sell_price
        self.gas_price = gas_price
        self.result = {}
        self.steam_price = steam_price
        self.gas_utility = {}
        self.CO2_utility = CO2_utility
        self.CO2_user = CO2_user
        self.PCC = {'maxP':10000,
                    'maxH':10000}
    def show(self):
        return {"sell_price":self.sell_price,
                "gas_price":self.gas_price}

class inverter:
    def __init__(self, ac_dc = 1 , dc_ac = 1,maxP = 5000):
        self.ac_dc_efficiency = ac_dc
        self.dc_ac_efficiency = dc_ac
        self.maxP = maxP
        self.result = {}
    def show(self):
        return {"ac_dc_efficiency":self.ac_dc_efficiency,
                "dc_ac_efficiency":self.dc_ac_efficiency}
