from importFile import *
from pyomo.bilevel import *
from microgrid_Model import *
import pandas as pd
import copy
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
H2M = 0.2
def DayAheadModel(microgrid_data,case,T_range,mode):
    microgrid_device = case.device
    N_T = case.NumOfTime
    T = range(len(T_range))
    step = 24 / N_T
    N_es = case.getKey(electricStorage)
    N_absc = case.getKey(absorptionChiller)
    N_bol = case.getKey(boiler)
    N_cs = case.getKey(coldStorage)
    N_ac = case.getKey(airConditioner)
    N_gt = case.getKey(gasTurbine)
    N_pv = case.getKey(PV)
    N_bus = case.graph.nodes()
    N_branch = case.graph.edges()
    if case.type == 'Simple':
        microgrid_device['ut'].buy_price = microgrid_data['电价'][T[0]:T[-1]+ 1].tolist()
        acLoad = microgrid_data['交流负荷'][T[0]:T[-1] + 1].tolist()
        dcLoad = microgrid_data['直流负荷'][T[0]:T[-1] + 1].tolist()
        cold_load = microgrid_data['冷负荷'][T[0]:T[-1]+1].tolist()
        water_heat_load = microgrid_data['热水负荷'][T[0]:T[-1]+1].tolist()
        steam_heat_load = microgrid_data['蒸汽负荷'][T[0]:T[-1]+1].tolist()
        pv_output = microgrid_data['光伏出力'][T[0]:T[-1] + 1].tolist()
    if case.type == 'Graph':
        # acLoad = dict()
        # for node in case.graph.nodes():
        #     acLoad[node] = microgrid_data[str(node)+'节点交流负荷'][T[0]:T[-1]+1].tolist()
        steam_heat_load = len(T)*[0]
        water_heat_load = len(T)*[0]
        cold_load = len(T)*[0]
    '''A general model and algorithm for microgrid optimal dispatch'''
    '''define sets'''
    optimalDispatch = ConcreteModel(name='IES_optimalDispatch')
    wind_power_max = microgrid_data['风机出力上限'][T[0]:T[-1]+1].tolist()
    wind_power_min = microgrid_data['风机出力下限'][T[0]:T[-1]+1].tolist()
    optimalDispatch.wp = Var(N_pv,T,bounds=lambda mdl,n,t: (wind_power_min[t], wind_power_max[t]))
    '''This is the sub-problem'''
    eps = 0.001 #精度
    optimalDispatch.sub = SubModel()
    optimalDispatch.sub.T = T
    optimalDispatch.sub.T_range = T_range
    optimalDispatch.sub.input = microgrid_data
    optimalDispatch.sub.case = case
    '''define variables'''
    # electrical storage
    optimalDispatch.sub.es_power = Var(N_es, T, bounds=lambda mdl, i, T: (-microgrid_device[i].Pmax_in, microgrid_device[i].Pmax_out))
    optimalDispatch.sub.es_energy = Var(N_es, T, bounds=lambda mdl, i, T: (
    microgrid_device[i].SOCmin * microgrid_device[i].capacity,
    microgrid_device[i].SOCmax * microgrid_device[i].capacity))
    # absorption chiller
    optimalDispatch.sub.absc_heat_in = Var(N_absc, T, bounds=lambda mdl, i, T: (0, microgrid_device[i].Hmax))
    # heat variables
    if sum(steam_heat_load) + sum(water_heat_load) > 1:
        optimalDispatch.sub.buy_heat = Var(T,bounds = (0,microgrid_device['ut'].PCC['maxH']))
        optimalDispatch.sub.medium_heat = Var(T,bounds=(-10000,10000))
        optimalDispatch.sub.low_heat = Var(T,bounds=(-10000,10000))
    # boiler
    optimalDispatch.sub.bol_power = Var(N_bol, T)
    optimalDispatch.sub.bol_constraint1 = Constraint(N_bol,T,rule = lambda  mdl,i,t: mdl.bol_power[i,t] <= microgrid_device[i].Pmax)
    optimalDispatch.sub.bol_constraint2 = Constraint(N_bol, T,rule=lambda mdl, i, t: mdl.bol_power[i, t] >= microgrid_device[i].Pmin)
    # cold storage
    optimalDispatch.sub.cs_cold = Var(N_cs, T, bounds=lambda mdl, i, T: (-microgrid_device[i].Hin, microgrid_device[i].Hout))
    optimalDispatch.sub.cs_cold_stored = Var(N_cs, T, bounds=lambda mdl, i, T: (
    microgrid_device[i].Tmin * microgrid_device[i].capacity, microgrid_device[i].Tmax * microgrid_device[i].capacity))
    # air conditioner
    optimalDispatch.sub.ac_power = Var(N_ac, T, bounds=lambda mdl, i, T: (0, microgrid_device[i].Pmax))
    # gas turbine
    optimalDispatch.sub.gt_power = Var(N_gt, T)
    optimalDispatch.sub.gt_constraint1 = Constraint(N_gt,T,rule = lambda mdl,i,t: mdl.gt_power[i,t] <= microgrid_device[i].Pmax)
    optimalDispatch.sub.gt_constraint2 = Constraint(N_gt, T,rule=lambda mdl, i, t: mdl.gt_power[i, t] >= microgrid_device[i].Pmin)
    # Injection Power
    # optimalDispatch.sub.P_inj = Var(N_bus,T)
    if case.type == 'Simple':
        # inverter
        optimalDispatch.sub.inv_dc = Var(T)  # inv_dc > 0 means energy flows from inverter to dc side
        # utility power
        optimalDispatch.sub.utility_power = Var(T, bounds=(-10000,10000))
    '''define constraints'''
    '''电功率平衡约束'''
    def ACPowerBalance(mdl,t):
        power_supply = sum(mdl.gt_power[i, t] for i in N_gt) \
                       + mdl.utility_power[t] + optimalDispatch.wp[t]
        power_demand = 1.05*sum(mdl.ac_power[i, t] for i in N_ac) \
                       + acLoad[t] + (1 / microgrid_device['inv'].ac_dc_efficiency) * mdl.inv_dc[t]\
					   + sum(microgrid_device[i].ElecCost * mdl.absc_heat_in[i, t] for i in N_absc)
        return -eps <= power_supply - power_demand <= eps

    def DCPowerBalance(mdl, t):
        power_supply = sum(mdl.es_power[i, t] for i in N_es) + mdl.inv_dc[t] + pv_output[t]
        power_demand = dcLoad[t]
        return -eps <= power_supply - power_demand <= eps

    def PowerBalance(mdl,t):
        power_supply = sum(mdl.gt_power[i, t] for i in N_gt) + sum(optimalDispatch.wp[i,t] for i in N_pv)
        power_demand = sum(case.graph.node[node]['Load'][t] for node in case.graph.nodes())
        return -eps <= power_supply - power_demand <= eps
    if case.type == 'Simple':
        optimalDispatch.sub.ACPowerBalance = Constraint(T, rule=ACPowerBalance)
        optimalDispatch.sub.DCPowerBalance = Constraint(T, rule=DCPowerBalance)
    if case.type == 'Graph':
        optimalDispatch.sub.PowerBalance = Constraint(T, rule=PowerBalance)
    '''热功率平衡约束'''
    if sum(steam_heat_load) + sum(water_heat_load) > 1:
        H2M = 0.2
        optimalDispatch.sub.HPB1 = Constraint(T, rule=lambda mdl, t: -eps <= mdl.medium_heat[t] - mdl.buy_heat[t] + sum(
            mdl.bol_power[n_bol, t] for n_bol in N_bol) + sum(mdl.gt_power[n_gt, t] * microgrid_device[n_gt].HER * microgrid_device[n_gt].heat_recycle for n_gt in N_gt) <= eps)
        optimalDispatch.sub.HPB2 = Constraint(T,rule = lambda mdl,t:mdl.medium_heat[t] >= steam_heat_load[t])
        optimalDispatch.sub.HPB3 = Constraint(T, rule=lambda mdl, t: -eps <= mdl.low_heat[t] - (H2M) * steam_heat_load[t] <= eps)
        optimalDispatch.sub.HPB4 = Constraint(T, rule=lambda mdl, t: mdl.low_heat[t] + mdl.medium_heat[t] >= water_heat_load[t] + steam_heat_load[t] + sum(mdl.absc_heat_in[n_absc, t] for n_absc in N_absc))

    '''冷功率平衡约束'''
    def coldPowerBalance(mdl, t):
        cold_supply = sum(mdl.ac_power[i, t] * microgrid_device[i].EER for i in N_ac) \
                      + sum(mdl.cs_cold[i, t] for i in N_cs) \
                      + sum(mdl.absc_heat_in[i, t] * microgrid_device[i].COP_htc for i in N_absc)
        cold_demand = cold_load[t]
        return -eps <= cold_supply - cold_demand <= eps
    if sum(cold_load) > 1:
        optimalDispatch.sub.coldPowerBalance = Constraint(T, rule=coldPowerBalance)
    '''电池日平衡约束、自放电率、爬坡率约束'''

    def batteryEnergyBalance(mdl, n_es, t):
        bat = microgrid_device[n_es]
        if t == T[0]:
            return -eps <= mdl.es_energy[n_es, t] - bat.SOCnow * bat.capacity <= eps
        else:
            return -eps <= mdl.es_energy[n_es, t] - mdl.es_energy[n_es, t - 1] * (1 - bat.selfRelease) \
                                             - step * mdl.es_power[n_es, t - 1] <= eps

    optimalDispatch.sub.batteryEnergyBalance = Constraint(N_es, T, rule=batteryEnergyBalance)
    optimalDispatch.sub.batteryEnergyBalance0 = Constraint(N_es, rule=lambda mdl,n:-eps <= mdl.es_energy[n, T[-1]]* (1 - microgrid_device[n].selfRelease) \
                                             - step * mdl.es_power[n, T[-1]] - microgrid_device[n].SOCint * microgrid_device[n].capacity <= eps)


    '''冰蓄冷日平衡约束、自放冷率、爬坡率约束'''

    def coldStorageEnergyBalance(mdl, n_cs, t):
        ice = microgrid_device[n_cs]
        if t == 0:
            return -eps <= mdl.cs_cold_stored[n_cs, t] - ice.Tint * ice.capacity <= eps
        else:
            return -eps <= mdl.cs_cold_stored[n_cs, t] - mdl.cs_cold_stored[n_cs, t - 1] * (1 - ice.selfRelease) \
                                                  - step * mdl.cs_cold[n_cs, t - 1] <= eps

    optimalDispatch.sub.coldStorageEnergyBalance = Constraint(N_cs, T, rule=coldStorageEnergyBalance)
    optimalDispatch.sub.coldStorageEnergyBalance0 = Constraint(N_cs, rule=lambda mdl,n:-eps <= mdl.cs_cold_stored[n, T[-1]]* (1 - microgrid_device[n].selfRelease) \
                                             - step * mdl.cs_cold[n, T[-1]] - microgrid_device[n].capacity*microgrid_device[n].Tint <= eps)
    '''燃气轮机/锅炉爬坡率约束'''
    def gtRampLimit(mdl,n,t):
        if t == 0:
            return Constraint.Skip
        else:
            return -microgrid_device[n].maxDetP <= mdl.gt_power[n,t] - mdl.gt_power[n,t-1] <= microgrid_device[n].maxDetP
    def bolRampLimit(mdl,n,t):
        if t == 0:
            return Constraint.Skip
        else:
            return -microgrid_device[n].maxDetP <= mdl.bol_power[n,t] - mdl.bol_power[n,t-1] <= microgrid_device[n].maxDetP
    optimalDispatch.sub.gtRampLimit = Constraint(N_gt,T,rule=gtRampLimit)
    # optimalDispatch.sub.bolRampLimit = Constraint(N_bol,T,rule=bolRampLimit)

    '''线路潮流约束'''
    #TODO 线路潮流约束还没写好
    def PFlimit(mdl,nf,nt,t):
        # nf = line[0]
        # nt = line[1]
        R = case.graph.edge[nf][nt]['R']
        X = case.graph.edge[nf][nt]['X']
        limit = case.graph.edge[nf][nt]['Limit']
        if limit is None:
            return Constraint.Skip
        else:
            PF = 1/X * (sum(case.B_INV.item(nf,nb) * mdl.P_inj[nb,t] for nb in N_bus) - sum(case.B_INV.item(nt,nb) * mdl.P_inj[nb,t] for nb in N_bus))
            return -limit <= PF <= limit
    # optimalDispatch.sub.PFlimit = Constraint(N_branch,T,rule=PFlimit)
    '''注入功率约束（中间量）'''
    def Power_Injection(mdl,nb,t):
        m = mdl.model()
        Temp = 0.0
        for key,dev in case.graph.node[nb]['device'].items():
            if isinstance(dev,gasTurbine):
                Temp += mdl.gt_power[key,t]
            if isinstance(dev,PV):
                Temp += m.wp[key,t]
        Temp -= case.graph.node[nb]['Load'][t]
        return -eps <= mdl.P_inj[nb,t] - Temp <= eps
    # optimalDispatch.sub.Power_Injection =  Constraint(N_bus,T,rule=Power_Injection)
    '''Define Objectives'''

    def OM_Cost(mdl):
        om_cost = 0
        for id in microgrid_device.keys():
            if (id in N_absc):
                om_cost += microgrid_device[id].om * step * sum(mdl.absc_heat_in[id, t] for t in T)
            if (id in N_bol):
                om_cost += microgrid_device[id].om * step * sum(mdl.bol_power[id, t] for t in T)
            if (id in N_ac):
                om_cost += microgrid_device[id].om * step * sum(mdl.ac_power[id, t] for t in T)
        return om_cost

    def Fuel_Cost(mdl):
        fuel_cost = 0
        for id in N_gt:
            fuel_cost += (1 / microgrid_device[id].efficiency) * microgrid_device['ut'].gas_price * step * sum(
                mdl.gt_power[id, t] for t in T)
        for id in N_bol:
            fuel_cost += (1 / microgrid_device[id].efficiency) * microgrid_device['ut'].gas_price * step * sum(
                mdl.bol_power[id, t] for t in T)
        return fuel_cost

    def ElectricalFee(mdl):
        return step * sum(mdl.utility_power[t] * microgrid_device['ut'].buy_price[t] for t in T)

    def HeatFee(mdl):
        return step * sum(mdl.buy_heat[t] for t in T) * microgrid_device['ut'].steam_price
    def obj_Economical(mdl):
        return OM_Cost(mdl) + Fuel_Cost(mdl) + ElectricalFee(mdl) + HeatFee(mdl)
    def obj_Efficiency(mdl):
        return (Fuel_Cost(mdl)/2.3 * 1.2143 + 0.1229 * 0.25 *sum(mdl.utility_power[t] for t in mdl.T) + 3.6 * 0.3412 * 0.25 * sum(mdl.buy_heat[t] for t in mdl.T)) \
               / (sum(acLoad)+sum(dcLoad)+sum(cold_load)+sum(water_heat_load)+sum(steam_heat_load))
    def obj_simple(mdl):
        return sum(sum(0.25*microgrid_device[n_gt].Cost*(mdl.gt_power[n_gt,t]) for n_gt in N_gt) for t in T)
    def cost_per_15min(mdl,t):
        return sum(0.25*microgrid_device[n_gt].Cost*(mdl.gt_power[n_gt,t]) for n_gt in N_gt)
    optimalDispatch.sub.cost_per_15min = cost_per_15min
    optimalDispatch.sub.obj_Economical = obj_Economical
    optimalDispatch.sub.obj_Efficiency = obj_Efficiency
    optimalDispatch.sub.obj_simple = obj_simple
    optimalDispatch.sub.objective = Objective(rule=obj_simple)
    if mode == 'max':
        optimalDispatch.objective = Objective(rule=lambda mdl:  -obj_simple(mdl.sub))
    elif mode == 'min':
        optimalDispatch.objective = Objective(rule=lambda mdl: obj_simple(mdl.sub))
    return optimalDispatch
def retriveResult(microgrid_data,case,mdl):
    model = mdl.sub
    microgrid_device = case.device
    N_T = case.NumOfTime
    T = model.T
    step = 24 / N_T
    N_es = case.getKey(electricStorage)
    N_absc = case.getKey(absorptionChiller)
    N_bol = case.getKey(boiler)
    N_cs = case.getKey(coldStorage)
    N_ac = case.getKey(airConditioner)
    N_gt = case.getKey(gasTurbine)
    N_pv = case.getKey(PV)
    acLoad = microgrid_data['交流负荷'].tolist()
    dcLoad = microgrid_data['直流负荷'].tolist()
    pv_output = microgrid_data['光伏出力'].tolist()
    microgrid_device['ut'].buy_price = microgrid_data['电价'].tolist()
    cold_load = microgrid_data['冷负荷'].tolist()
    water_heat_load = microgrid_data['热水负荷'].tolist()
    steam_heat_load = microgrid_data['蒸汽负荷'].tolist()
    ts = pd.date_range('2017/8/28 00:00:00', periods=96, freq='15min')
    df = pd.DataFrame()
    df['交流负荷'] = pd.Series([acLoad[t] for t in T],index=T)
    df['直流负荷'] = pd.Series([dcLoad[t] for t in T],index=T)
    df['风机出力'] = pd.Series(value(mdl.wp[t] for t in T),index=T)
    df['冷负荷'] = microgrid_data['冷负荷'].loc[T]
    df['蒸汽负荷'] =  microgrid_data['蒸汽负荷'].loc[T]
    df['电价'] = pd.Series(microgrid_device['ut'].buy_price).loc[T]
    df['购电功率'] = pd.Series([value(model.utility_power[t]) for t in T],index=T)
    df['购热功率'] = pd.Series([value(model.buy_heat[t]) for t in T],index=T)
    df['交流侧逆变器功率'] = pd.Series([value(model.inv_ac[t]) for t in T],index=T)
    df['直流侧逆变器功率'] = pd.Series([value(model.inv_dc[t]) for t in T],index=T)
    for es in N_es:
        df[es + '电池电量'] = pd.Series([value(model.es_energy[es, t]) for t in T],index=T)
        df[es + '电储能充电功率'] = pd.Series([value(model.es_power_in[es, t]) for t in T],index=T)
        df[es + '电储能放电功率'] = pd.Series([value(model.es_power_out[es, t]) for t in T],index=T)
    for gt in N_gt:
        df[gt + '机组出力'] = pd.Series([value(model.gt_power[gt, t]) for t in T],index=T)
        df[gt + '余热锅炉中品位热功率'] = pd.Series(
            [value(model.gt_power[gt, t]) * microgrid_device[gt].HER * microgrid_device[gt].heat_recycle for t
             in T],index=T)
        df[gt + '余热锅炉低品位热功率'] = pd.Series(
            [value(model.gt_power[gt, t]) * microgrid_device[gt].HER * microgrid_device[gt].low_heat_recycle for t
             in T], index=T)
    df['光伏出力'] = pd.Series(pv_output)
    for ac in N_ac:
        df[ac + '空调制冷耗电功率'] = pd.Series([value(model.ac_power[ac, t]) for t in T],index=T)
        df[ac + '空调制冷功率'] = df[ac + '空调制冷耗电功率'] * microgrid_device[ac].EER
    for cs in N_cs:
        df[cs + '冰蓄冷耗电功率'] = pd.Series([value(model.cs_power[cs, t]) for t in T],index=T)
        df[cs + '冰蓄冷储冷功率'] = pd.Series([value(model.cs_cold_in[cs, t]) for t in T],index=T)
        df[cs + '冰蓄冷供冷功率'] = pd.Series([value(model.cs_cold_out[cs, t]) for t in T],index=T)
        df[cs + '冰蓄冷制冷机直接供冷耗电功率'] = df[cs + '冰蓄冷耗电功率'] * microgrid_device[cs].EER - df[cs + '冰蓄冷储冷功率']
        df[cs + '冰蓄冷储冷量'] = pd.Series([value(model.cs_cold_stored[cs, t]) for t in T],index=T)
    for absc in N_absc:
        df[absc + '吸收式制冷机制冷功率'] = pd.Series([value(model.absc_heat_in[absc, t]) for t in T],index=T) * \
                                  microgrid_device[absc].COP_htc
    for bol in N_bol:
        df[bol + '燃气锅炉热功率'] = pd.Series([value(model.bol_power[bol, t]) for t in T],index=T)
    df['中品位热功率'] = pd.Series([value(model.buy_heat[t]) for t in T],index=T)
    #df['低品位热功率'] = pd.Series([value(model.low_heat[t]) for t in T],index=T) + pd.Series([value(model.medium_heat[t]) for t in T],index=T)
    '''demond response'''
    try:
        if model.mode == 'E':
            df['期望电功率'] = pd.Series(model.P_ref)
        elif model.mode == 'H':
            df['期望热功率'] = pd.Series(model.H_ref)
            df['可调负荷增加热功率']=pd.Series([value(model.DRHeatLoad[t] for t in model.peak)],index=T)
    except Exception as e:
        pass
    return df

def df_sum(df,cols):
    newdf = pd.Series([0]*df.__len__(),index=df[cols[0]].index)
    for col in cols:
        newdf = newdf + df[col]
    return newdf
