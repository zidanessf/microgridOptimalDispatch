from importFile import *
from microgrid import *
import pandas as pd
import matplotlib.pyplot as plt
import datetime
optimalDispatch = ConcreteModel(name = 'IES_optimalDispatch')
'''Initialize a special case of microgrid'''
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
'''Load input data'''
microgrid_data = pd.read_excel('input.xlsx')
acLoad = microgrid_data['交流负荷'].tolist()
dcLoad =  microgrid_data['直流负荷'].tolist()
pv_output = microgrid_data['光伏出力'].tolist()
microgrid_device['ut'].buy_price = microgrid_data['电价'].tolist()
cold_load = microgrid_data['冷负荷'].tolist()
water_heat_load = microgrid_data['热水负荷'].tolist()
steam_heat_load = microgrid_data['蒸汽负荷'].tolist()
'''A general model and algorithm for microgrid optimal dispatch'''
'''define sets'''
N_T = 95
T = RangeSet(0,N_T)
step = 24/N_T
N_es = [key for key in microgrid_device.keys() if isinstance(microgrid_device[key],electricStorage)]
N_absc = [key for key in microgrid_device.keys() if isinstance(microgrid_device[key],absorptionChiller)]
N_bol = [key for key in microgrid_device.keys() if isinstance(microgrid_device[key],boiler)]
N_cs = [key for key in microgrid_device.keys() if isinstance(microgrid_device[key],coldStorage)]
N_ac = [key for key in microgrid_device.keys() if isinstance(microgrid_device[key],airConditioner)]
N_gt = [key for key in microgrid_device.keys() if isinstance(microgrid_device[key],gasTurbine)]
N_pv = [key for key in microgrid_device.keys() if isinstance(microgrid_device[key],PV)]
N_dev = [key for key in microgrid_device.keys()]
'''define variables'''
#electrical storage
optimalDispatch.es_power_in = Var(N_es,T,bounds=lambda mdl,i,T:(0,microgrid_device[i].Pmax_in))
optimalDispatch.es_power_out =Var(N_es,T,bounds=lambda mdl,i,T:(0,microgrid_device[i].Pmax_out))
optimalDispatch.es_energy = Var(N_es,T,bounds=lambda mdl,i,T:(microgrid_device[i].SOCmin*microgrid_device[i].capacity,microgrid_device[i].SOCmax*microgrid_device[i].capacity))
#absorption chiller
optimalDispatch.absc_heat_in = Var(N_absc,T,domain=PositiveReals)
#heat variables
optimalDispatch.high_heat = Var(T)
optimalDispatch.medium_heat = Var(T)
optimalDispatch.low_heat = Var(T)
#boiler
optimalDispatch.bol_power = Var(N_bol,T,bounds=lambda mdl,i,T:(0,microgrid_device[i].Pmax))
#cold storage
optimalDispatch.cs_power = Var(N_cs,T,bounds=lambda mdl,i,T:(0,microgrid_device[i].Pmax))
optimalDispatch.cs_cold_in = Var(N_cs,T,bounds=lambda mdl,i,T:(0,microgrid_device[i].Hin))
optimalDispatch.cs_cold_out = Var(N_cs,T,bounds=lambda mdl,i,T:(0,microgrid_device[i].Hout))
optimalDispatch.cs_cold_stored = Var(N_cs,T,bounds=lambda mdl,i,T:(microgrid_device[i].Tmin*microgrid_device[i].capacity,microgrid_device[i].Tmax*microgrid_device[i].capacity))
#air conditioner
optimalDispatch.ac_power = Var(N_ac,T,bounds=lambda mdl,i,T:(0,microgrid_device[i].Pmax))
#gas turbine
optimalDispatch.gt_power = Var(N_gt,T,bounds=lambda mdl,i,T:(microgrid_device[i].Pmin,microgrid_device[i].Pmax))
#inverter
optimalDispatch.inv_ac = Var(T,bounds=(0,microgrid_device['inv'].maxP)) #inv_ac > 0 means energy flows from inverter to ac side
optimalDispatch.inv_dc = Var(T,bounds=(0,microgrid_device['inv'].maxP)) #inv_dc > 0 means energy flows from inverter to dc side
#utility power
optimalDispatch.utility_power = Var(T,domain=PositiveReals)
#TODO add ac-dc inverter variables
'''define disjuncts(states)'''
'''Battery'''
def es_in_out_state(b,i,T,indicator):
    mdl = b.model()
    if indicator == 0:
        b.es_forbidden = Constraint(expr=mdl.es_power_in[i,T]==0)
    else:
        b.es_forbidden = Constraint(expr=mdl.es_power_out[i, T] == 0)
optimalDispatch.es_power_in_out_state = Disjunct(N_es,T,[0,1],rule=es_in_out_state)
def es_in_out(mdl,i,T):
    return [mdl.es_power_in_out_state[i,T,0],mdl.es_power_in_out_state[i,T,1]]
optimalDispatch.es_in_out = Disjunction(N_es,T,rule=es_in_out)
'''Cold Storage'''
def cs_cold_in_out_state(b,i,T,indicator):
    mdl = b.model()
    if indicator == 0:
        b.cs_forbidden = Constraint(expr=mdl.cs_cold_in[i, T] == 0)
    else:
        b.cs_forbidden = Constraint(expr=mdl.cs_cold_out[i, T] == 0)
optimalDispatch.cs_cold_in_out_state = Disjunct(N_cs,T,[0,1],rule=cs_cold_in_out_state)
def cs_in_out(mdl,i,T):
    return [mdl.cs_cold_in_out_state[i,T,0],mdl.cs_cold_in_out_state[i,T,1]]
optimalDispatch.cs_in_out = Disjunction(N_cs,T,rule=cs_in_out)
'''INVERTER'''

def inv_trans_state(b,T,indicator):
    mdl = b.model()
    if indicator == 0: # ac to dc
        b.acdc_state = Constraint(expr=mdl.inv_ac[T] == 0)
    else:
        b.acdc_state = Constraint(expr=mdl.inv_dc[T] == 0)
optimalDispatch.inv_trans_state = Disjunct(T,[0,1],rule = inv_trans_state)
def inv_ac2dc_dc2ac(mdl,T):
    return [mdl.inv_trans_state[T,0],mdl.inv_trans_state[T,1]]
optimalDispatch.inv_ac2dc_dc2ac = Disjunction(T,rule=inv_ac2dc_dc2ac)

'''define constraints'''
'''电功率平衡约束'''
def ACPowerBalance(mdl,t):
    power_supply = sum(mdl.gt_power[i,t] for i in N_gt) \
                   + mdl.utility_power[t] + mdl.inv_ac[t]
    power_demand =  sum(mdl.cs_power[i,t] for i in N_cs) \
                   + sum(mdl.ac_power[i,t] for i in N_ac) \
                   + acLoad[t] + (1/microgrid_device['inv'].ac_dc_efficiency)*mdl.inv_dc[t]
    return power_supply == power_demand
def DCPowerBalance(mdl,t):
    power_supply = sum(mdl.es_power_out[i,t] for i in N_es) + mdl.inv_dc[t] + pv_output[t]
    power_demand = dcLoad[t] + sum(mdl.es_power_in[i,t] for i in N_es) + (1/microgrid_device['inv'].dc_ac_efficiency)*mdl.inv_ac[t]
    return power_supply == power_demand
optimalDispatch.ACPowerBalance = Constraint(T,rule=ACPowerBalance)
optimalDispatch.DCPowerBalance = Constraint(T,rule=DCPowerBalance)
'''热功率平衡约束'''
H2M = 0.3
M2L = 0.3
def heatPowerBalance(mdl,t):
    heat_supply = sum(mdl.bol_power[i,t] for i in N_bol) \
                  + sum(microgrid_device[i].HER*microgrid_device[i].heat_recycle*mdl.gt_power[i,t] for i in N_gt)
    heat_demand = sum(mdl.absc_heat_in[i,t] for i in N_absc) \
                  + water_heat_load[t]
    return heat_supply >= heat_demand
optimalDispatch.heatPowerBalance = Constraint(T,rule=heatPowerBalance)
optimalDispatch.HPB1 = Constraint(T,rule=lambda mdl,t:mdl.high_heat[t]+sum(mdl.bol_power[n_bol,t] for n_bol in N_bol)>=steam_heat_load[t])
optimalDispatch.HPB2 = Constraint(T,rule=lambda mdl,t:mdl.medium_heat[t]==H2M*steam_heat_load[t]+sum(mdl.gt_power[n_gt,t]*microgrid_device[n_gt].HER*microgrid_device[n_gt].heat_recycle for n_gt in N_gt))
optimalDispatch.HPB3 = Constraint(T,rule=lambda mdl,t:mdl.medium_heat[t] >= sum(mdl.absc_heat_in[n_absc,t] for n_absc in N_absc))
optimalDispatch.HPB4 = Constraint(T,rule=lambda mdl,t:mdl.low_heat[t] == (M2L-1) * sum(mdl.absc_heat_in[n_absc,t] for n_absc in N_absc) + mdl.medium_heat[t])
optimalDispatch.HPB5 = Constraint(T,rule=lambda mdl,t:mdl.low_heat[t] >= water_heat_load[t])
#TODO fix bugs
'''冷功率平衡约束'''
def coldPowerBalance(mdl,t):
    cold_supply = sum(mdl.ac_power[i,t]*microgrid_device[i].EER for i in N_ac) \
                  +sum((mdl.cs_power[i,t]*microgrid_device[i].EER-mdl.cs_cold_in[i,t]) for i in N_cs)\
                  +sum(mdl.cs_cold_out[i,t] for i in N_cs)\
                  +sum(mdl.absc_heat_in[i,t] * microgrid_device[i].COP_htc for i in N_absc)
    cold_demand = cold_load[t]
    return cold_supply == cold_demand
optimalDispatch.coldPowerBalance = Constraint(T,rule=coldPowerBalance)
optimalDispatch.ChillerMoreThanColdIn = Constraint(T,N_cs,rule=lambda mdl,t,n_cs:mdl.cs_power[n_cs,t]*microgrid_device[n_cs].EER>=mdl.cs_cold_in[n_cs,t])
'''电池日平衡约束、自放电率、爬坡率约束'''
def batteryEnergyBalance(mdl,n_es,t):
    bat = microgrid_device[n_es]
    if t == 95:
        return mdl.es_energy[n_es,t] == mdl.es_energy[n_es,0]
    elif t == 0:
        return mdl.es_energy[n_es,t] == bat.SOCint * bat.capacity
    else:
        return mdl.es_energy[n_es,t] == mdl.es_energy[n_es,t - 1] * (1-bat.selfRelease) \
                                        + step * (bat.efficiency*mdl.es_power_in[n_es,t-1] - (1/bat.efficiency)*mdl.es_power_out[n_es,t-1])
optimalDispatch.batteryEnergyBalance = Constraint(N_es,T,rule=batteryEnergyBalance)
def batteryRampLimit(mdl,n_es,t):
    if t == 0:
        return Constraint.Skip
    else:
        return -microgrid_device[n_es].maxDetP <= (mdl.es_power_out[n_es,t] - mdl.es_power_in[n_es,t]) - (mdl.es_power_out[n_es,t-1] - mdl.es_power_in[n_es,t-1]) <= microgrid_device[n_es].maxDetP
optimalDispatch.batteryRampLimit = Constraint(N_es,T,rule=batteryRampLimit)
'''冰蓄冷日平衡约束、自放冷率、爬坡率约束'''
def coldStorageEnergyBalance(mdl,n_cs,t):
    ice = microgrid_device[n_cs]
    if t == 0:
        return mdl.cs_cold_stored[n_cs,t] == ice.Tint * ice.capacity
    elif t == 95:
        return mdl.cs_cold_stored[n_cs,t] == mdl.cs_cold_stored[n_cs,0]
    else:
        return mdl.cs_cold_stored[n_cs,t] == mdl.cs_cold_stored[n_cs,t-1]*(1-ice.selfRelease)\
                                             + step * (ice.efficiency*mdl.cs_cold_in[n_cs,t-1] - (1/ice.efficiency)*mdl.cs_cold_out[n_cs,t-1])
optimalDispatch.coldStorageEnergyBalance = Constraint(N_cs,T,rule=coldStorageEnergyBalance)
def coldStorageRampLimit(mdl,n_cs,t):
    if t == 0:
        return Constraint.Skip
    else:
        return -microgrid_device[n_cs].maxDetP <= (mdl.cs_cold_out[n_cs,t] - mdl.cs_cold_in[n_cs,t]) - (mdl.cs_cold_out[n_cs,t-1] - mdl.cs_cold_in[n_cs,t-1]) <= microgrid_device[n_cs].maxDetP
optimalDispatch.coldStorageRampLimit = Constraint(N_cs,T,rule=coldStorageRampLimit)

'''Define Objectives'''
def OM_Cost(mdl):
    om_cost = 0
    for id in N_dev:
        if(id in N_es):
            om_cost += microgrid_device[id].om * step * (sum(mdl.es_power_out[id,t] for t in T)+sum(mdl.es_power_in[id,t] for t in T))
        if(id in N_cs):
            om_cost += microgrid_device[id].om * step * (sum(mdl.cs_cold_out[id,t] for t in T)+sum(mdl.cs_cold_in[id,t] for t in T))
        if(id in N_absc):
            om_cost += microgrid_device[id].om * step * sum(mdl.absc_heat_in[id,t] for t in T)
        if(id in N_bol):
            om_cost += microgrid_device[id].om * step * sum(mdl.bol_power[id,t] for t in T)
        if(id in N_ac):
            om_cost += microgrid_device[id].om * step * sum(mdl.ac_power[id,t] for t in T)
    return om_cost
def Dep_Cost(mdl):
    return sum(sum(microgrid_device[id].Cbw * step * mdl.es_power_out[id,t] for id in N_es) for t in T)
def Fuel_Cost(mdl):
    fuel_cost = 0
    for id in N_gt:
        fuel_cost += (1/microgrid_device[id].efficiency) * microgrid_device['ut'].gas_price * step * sum(mdl.gt_power[id,t] for t in T)
    for id in N_bol:
        fuel_cost += (1/microgrid_device[id].efficiency) * microgrid_device['ut'].gas_price * step * sum(mdl.bol_power[id,t] for t in T)
    return  fuel_cost
def ElectricalFee(mdl):
    return step * sum(mdl.utility_power[t]*microgrid_device['ut'].buy_price[t] for t in T)
def HeatFee(mdl):
    return  step*sum(mdl.high_heat[t] for t in T)*microgrid_device['ut'].steam_price
def obj_Economical(mdl):
    return OM_Cost(mdl) + Dep_Cost(mdl) + Fuel_Cost(mdl) + ElectricalFee(mdl) + HeatFee(mdl)
optimalDispatch.obj_Economical = Objective(rule=obj_Economical)

'''Solve the Model'''
xfrm = TransformationFactory('gdp.chull')
xfrm.apply_to(optimalDispatch)
solver = SolverFactory('glpk')
solver.solve(optimalDispatch)
'''Retrieve the result'''
ts=pd.date_range('2017/8/28 00:00:00',periods=96,freq='15min')
df = pd.DataFrame()
df['交流负荷'] = pd.Series(acLoad)
df['直流负荷'] = pd.Series(dcLoad)
df['电价'] = pd.Series(microgrid_device['ut'].buy_price)
df['电网购电功率'] = pd.Series([value(optimalDispatch.utility_power[t]) for t in T])
df['电网购热功率'] = pd.Series([value(optimalDispatch.high_heat[t]) for t in T])
df['交流侧逆变器功率'] = pd.Series([value(optimalDispatch.inv_ac[t]) for t in T])
df['直流侧逆变器功率'] = pd.Series([value(optimalDispatch.inv_dc[t]) for t in T])
for es in N_es:
    df[es+'电池电量'] = pd.Series([value(optimalDispatch.es_energy[es, t]) for t in T])
    df[es+'电储能充电功率'] = pd.Series([value(optimalDispatch.es_power_in[es,t]) for t in T])
    df[es+'电储能放电功率'] = pd.Series([value(optimalDispatch.es_power_out[es, t]) for t in T])
for gt in N_gt:
    df[gt+'机组出力'] = pd.Series([value(optimalDispatch.gt_power[gt, t]) for t in T])
    df[gt+'余热锅炉热功率'] = pd.Series([value(optimalDispatch.gt_power[gt, t])*microgrid_device[gt].HER*microgrid_device[gt].heat_recycle for t in T])
df['光伏出力'] = pd.Series(pv_output)
for ac in N_ac:
    df[ac+'空调制冷耗电功率'] = pd.Series([value(optimalDispatch.ac_power[ac, t]) for t in T])
    df[ac+'空调制冷功率'] = df[ac+'空调制冷耗电功率']*microgrid_device[ac].EER
for cs in N_cs:
    df[cs+'冰蓄冷耗电功率'] = pd.Series([value(optimalDispatch.cs_power[cs, t]) for t in T])
    df[cs+'冰蓄冷储冷功率'] = pd.Series([value(optimalDispatch.cs_cold_in[cs, t]) for t in T])
    df[cs+'冰蓄冷供冷功率'] = pd.Series([value(optimalDispatch.cs_cold_out[cs, t]) for t in T])
    df[cs+'冰蓄冷制冷机直接供冷耗电功率'] = df[cs+'冰蓄冷耗电功率']*microgrid_device[cs].EER - df[cs+'冰蓄冷储冷功率']
    df[cs+'冰蓄冷储冷量'] = pd.Series([value(optimalDispatch.cs_cold_stored[cs, t]) for t in T])
for absc in N_absc:
    df[absc+'吸收式制冷机制冷功率'] = pd.Series([value(optimalDispatch.absc_heat_in[absc, t]) for t in T]) * microgrid_device[absc].COP_htc
for bol in N_bol:
    df[bol+'燃气锅炉热功率'] = pd.Series([value(optimalDispatch.bol_power[bol, t]) for t in T])
df['高品位热功率'] = pd.Series([value(optimalDispatch.high_heat[t]) for t in T])
df['中品位热功率'] = pd.Series([value(optimalDispatch.medium_heat[t]) for t in T])
df['低品位热功率'] = pd.Series([value(optimalDispatch.low_heat[t]) for t in T])
df.to_excel('output.xlsx')
plt.plot(df)
plt.show()


