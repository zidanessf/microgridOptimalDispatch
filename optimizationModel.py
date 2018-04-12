from importFile import *
from pyomo.bilevel import *
from microgrid_Model import *
import pandas as pd
import copy
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
H2M = 0.2
def DayAheadModel(microgrid_data,case,T_range):
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
    acLoad = microgrid_data['交流负荷'][T[0]:T[-1]+1].tolist()
    dcLoad = microgrid_data['直流负荷'][T[0]:T[-1]+1].tolist()
    pv_output = microgrid_data['光伏出力'][T[0]:T[-1]+1].tolist()
    microgrid_device['ut'].buy_price = microgrid_data['电价'][T[0]:T[-1]+1].tolist()
    cold_load = microgrid_data['冷负荷'][T[0]:T[-1]+1].tolist()
    water_heat_load = microgrid_data['热水负荷'][T[0]:T[-1]+1].tolist()
    steam_heat_load = microgrid_data['蒸汽负荷'][T[0]:T[-1]+1].tolist()
    '''A general model and algorithm for microgrid optimal dispatch'''
    '''define sets'''
    optimalDispatch = ConcreteModel(name='IES_optimalDispatch')
    optimalDispatch.T = T
    optimalDispatch.T_range = T_range
    optimalDispatch.input = microgrid_data
    optimalDispatch.case = case
    '''define variables'''
    # electrical storage
    optimalDispatch.es_power_in = Var(N_es, T, bounds=lambda mdl, i, T: (0, microgrid_device[i].Pmax_in))
    optimalDispatch.es_power_out = Var(N_es, T, bounds=lambda mdl, i, T: (0, microgrid_device[i].Pmax_out))
    #optimalDispatch.es_power_out_0 = Constraint(N_es,rule=lambda mdl,i: mdl.es_power_out[i,T[-1]] == 0)
    optimalDispatch.es_energy = Var(N_es, T, bounds=lambda mdl, i, T: (
    microgrid_device[i].SOCmin * microgrid_device[i].capacity,
    microgrid_device[i].SOCmax * microgrid_device[i].capacity))
    # absorption chiller
    optimalDispatch.absc_heat_in = Var(N_absc, T, bounds=lambda mdl, i, T: (0, microgrid_device[i].Hmax))
    # heat variables
    optimalDispatch.buy_heat = Var(T,bounds = (0,microgrid_device['ut'].PCC['maxH']))
    optimalDispatch.medium_heat = Var(T,bounds=(-10000,10000))
    optimalDispatch.low_heat = Var(T,bounds=(-10000,10000))
    # boiler
    optimalDispatch.bol_power = Var(N_bol, T)
    optimalDispatch.bol_state = Var(N_bol, T, within=Binary)
    optimalDispatch.bol_auxvar = Var(N_bol,T)
    optimalDispatch.bol_constraint1 = Constraint(N_bol,T,rule = lambda  mdl,i,t: mdl.bol_power[i,t] <= mdl.bol_state[i,t]*microgrid_device[i].Pmax)
    optimalDispatch.bol_constraint2 = Constraint(N_bol, T,rule=lambda mdl, i, t: mdl.bol_power[i, t] >= mdl.bol_state[i, t] *microgrid_device[i].Pmin)
    # cold storage
    optimalDispatch.cs_power = Var(N_cs, T, bounds=lambda mdl, i, T: (0, microgrid_device[i].Pmax))
    optimalDispatch.cs_cold_in = Var(N_cs, T, bounds=lambda mdl, i, T: (0, microgrid_device[i].Hin))
    optimalDispatch.cs_cold_out = Var(N_cs, T, bounds=lambda mdl, i, T: (0, microgrid_device[i].Hout))
    #optimalDispatch.cs_cold_out_0 = Constraint(N_cs, rule = lambda  mdl,i: mdl.cs_cold_out[i,T[-1]] == 0)
    optimalDispatch.cs_cold_stored = Var(N_cs, T, bounds=lambda mdl, i, T: (
    microgrid_device[i].Tmin * microgrid_device[i].capacity, microgrid_device[i].Tmax * microgrid_device[i].capacity))
    # air conditioner
    optimalDispatch.ac_power = Var(N_ac, T, bounds=lambda mdl, i, T: (0, microgrid_device[i].Pmax))
    # gas turbine
    optimalDispatch.gt_power = Var(N_gt, T)
    optimalDispatch.gt_state = Var(N_gt,T,within=Binary)
    optimalDispatch.gt_auxvar = Var(N_gt, T)
    optimalDispatch.gt_constraint1 = Constraint(N_gt,T,rule = lambda mdl,i,t: mdl.gt_power[i,t] <= mdl.gt_state[i,t]*microgrid_device[i].Pmax)
    optimalDispatch.gt_constraint2 = Constraint(N_gt, T,rule=lambda mdl, i, t: mdl.gt_power[i, t] >= mdl.gt_state[i, t] *microgrid_device[i].Pmin)
    # inverter
    optimalDispatch.inv_ac = Var(T, bounds=(
    0, microgrid_device['inv'].maxP))  # inv_ac > 0 means energy flows from inverter to ac side
    optimalDispatch.inv_dc = Var(T, bounds=(
    0, microgrid_device['inv'].maxP))  # inv_dc > 0 means energy flows from inverter to dc side
    # utility power
    optimalDispatch.utility_power = Var(T, domain=PositiveReals)

    '''define disjuncts(states)'''
    '''Battery'''
    def es_in_out(mdl,i,t):
        return complements(mdl.es_power_in[i,t] >= 0 , mdl.es_power_out[i,t] >= 0)
    optimalDispatch.es_in_out = Complementarity(N_es,T,rule=es_in_out)

    '''Cold Storage'''
    def cs_in_out(mdl,i,t):
        return complements(mdl.cs_cold_in[i,t] >= 0 , mdl.cs_cold_out[i,t] >= 0)
    optimalDispatch.cs_in_out = Complementarity(N_cs,T,rule=cs_in_out)

    '''INVERTER'''
    def inv_ac2dc_dc2ac(mdl,t):
        return complements(mdl.inv_ac[t] >= 0 , mdl.inv_dc[t] >= 0)
    optimalDispatch.inv_ac2dc_dc2ac = Complementarity(T,rule=inv_ac2dc_dc2ac)

    '''define constraints'''
    '''电功率平衡约束'''

    def ACPowerBalance(mdl, t):
        power_supply = sum(mdl.gt_power[i, t] for i in N_gt) \
                       + mdl.utility_power[t] + mdl.inv_ac[t]
        power_demand = 1.05*sum(mdl.cs_power[i, t] for i in N_cs) \
                       + 1.05*sum(mdl.ac_power[i, t] for i in N_ac) \
                       + acLoad[t] + (1 / microgrid_device['inv'].ac_dc_efficiency) * mdl.inv_dc[t]\
					   + sum(microgrid_device[i].ElecCost * mdl.absc_heat_in[i, t] for i in N_absc)
        return power_supply == power_demand

    def DCPowerBalance(mdl, t):
        power_supply = sum(mdl.es_power_out[i, t] for i in N_es) + mdl.inv_dc[t] + pv_output[t]
        power_demand = dcLoad[t] + sum(mdl.es_power_in[i, t] for i in N_es) + (1 / microgrid_device[
            'inv'].dc_ac_efficiency) * mdl.inv_ac[t]
        return power_supply == power_demand

    optimalDispatch.ACPowerBalance = Constraint(T, rule=ACPowerBalance)
    optimalDispatch.DCPowerBalance = Constraint(T, rule=DCPowerBalance)
    '''热功率平衡约束'''
    H2M = 0.2

    optimalDispatch.HPB1 = Constraint(T, rule=lambda mdl, t: mdl.medium_heat[t] == mdl.buy_heat[t] + sum(
        mdl.bol_power[n_bol, t] for n_bol in N_bol) + sum(mdl.gt_power[n_gt, t] * microgrid_device[n_gt].HER * microgrid_device[n_gt].heat_recycle for n_gt in N_gt))
    optimalDispatch.HPB2 = Constraint(T,rule = lambda mdl,t:mdl.medium_heat[t] >= steam_heat_load[t])
    optimalDispatch.HPB3 = Constraint(T, rule=lambda mdl, t: mdl.low_heat[t] == (H2M) * steam_heat_load[t])
    optimalDispatch.HPB4 = Constraint(T, rule=lambda mdl, t: mdl.low_heat[t] + mdl.medium_heat[t] >= water_heat_load[t] + steam_heat_load[t] + sum(mdl.absc_heat_in[n_absc, t] for n_absc in N_absc))
    # TODO 完善高中低品味热模型
    '''冷功率平衡约束'''

    def coldPowerBalance(mdl, t):
        cold_supply = sum(mdl.ac_power[i, t] * microgrid_device[i].EER for i in N_ac) \
                      + sum((mdl.cs_power[i, t] * microgrid_device[i].EER - mdl.cs_cold_in[i, t]) for i in N_cs) \
                      + sum(mdl.cs_cold_out[i, t] for i in N_cs) \
                      + sum(mdl.absc_heat_in[i, t] * microgrid_device[i].COP_htc for i in N_absc)
        cold_demand = cold_load[t]
        return cold_supply == cold_demand

    optimalDispatch.coldPowerBalance = Constraint(T, rule=coldPowerBalance)
    optimalDispatch.ChillerMoreThanColdIn = Constraint(T, N_cs, rule=lambda mdl, t, n_cs: mdl.cs_power[n_cs, t] *
                                                                                          microgrid_device[n_cs].EER >=
                                                                                          mdl.cs_cold_in[n_cs, t])
    '''电池日平衡约束、自放电率、爬坡率约束'''

    def batteryEnergyBalance(mdl, n_es, t):
        bat = microgrid_device[n_es]
        if t == T[0]:
            return mdl.es_energy[n_es, t] == bat.SOCnow * bat.capacity
        else:
            return mdl.es_energy[n_es, t] == mdl.es_energy[n_es, t - 1] * (1 - bat.selfRelease) \
                                             + step * (
            bat.efficiency * mdl.es_power_in[n_es, t - 1] - (1 / bat.efficiency) * mdl.es_power_out[n_es, t - 1])

    optimalDispatch.batteryEnergyBalance = Constraint(N_es, T, rule=batteryEnergyBalance)
    optimalDispatch.batteryEnergyBalance0 = Constraint(N_es, rule=lambda mdl,n:mdl.es_energy[n, T[-1]]* (1 - microgrid_device[n].selfRelease) \
                                             + step * (microgrid_device[n].efficiency * mdl.es_power_in[n, T[-1]] - (1 / microgrid_device[n].efficiency) * mdl.es_power_out[n, T[-1]]) == microgrid_device[n].SOCint * microgrid_device[n].capacity )

    def batteryRampLimit(mdl, n_es, t):
        if t == 0:
            return Constraint.Skip
        else:
            return -microgrid_device[n_es].maxDetP <= (mdl.es_power_out[n_es, t] - mdl.es_power_in[n_es, t]) - (
            mdl.es_power_out[n_es, t - 1] - mdl.es_power_in[n_es, t - 1]) <= microgrid_device[n_es].maxDetP

    optimalDispatch.batteryRampLimit = Constraint(N_es, T, rule=batteryRampLimit)
    '''冰蓄冷日平衡约束、自放冷率、爬坡率约束'''

    def coldStorageEnergyBalance(mdl, n_cs, t):
        ice = microgrid_device[n_cs]
        if t == 0:
            return mdl.cs_cold_stored[n_cs, t] == ice.Tint * ice.capacity
        else:
            return mdl.cs_cold_stored[n_cs, t] == mdl.cs_cold_stored[n_cs, t - 1] * (1 - ice.selfRelease) \
                                                  + step * (
            ice.efficiency * mdl.cs_cold_in[n_cs, t - 1] - (1 / ice.efficiency) * mdl.cs_cold_out[n_cs, t - 1])

    optimalDispatch.coldStorageEnergyBalance = Constraint(N_cs, T, rule=coldStorageEnergyBalance)
    optimalDispatch.coldStorageEnergyBalance0 = Constraint(N_cs, rule=lambda mdl,n:mdl.cs_cold_stored[n, T[-1]]* (1 - microgrid_device[n].selfRelease) \
                                             + step * (microgrid_device[n].efficiency * mdl.cs_cold_in[n, T[-1]] - (1 / microgrid_device[n].efficiency) * mdl.cs_cold_out[n, T[-1]]) == microgrid_device[n].capacity*microgrid_device[n].Tint)
    def coldStorageRampLimit(mdl, n_cs, t):
        if t == 0:
            return Constraint.Skip
        else:
            return -microgrid_device[n_cs].maxDetP <= (mdl.cs_cold_out[n_cs, t] - mdl.cs_cold_in[n_cs, t]) - (
            mdl.cs_cold_out[n_cs, t - 1] - mdl.cs_cold_in[n_cs, t - 1]) <= microgrid_device[n_cs].maxDetP

    optimalDispatch.coldStorageRampLimit = Constraint(N_cs, T, rule=coldStorageRampLimit)
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
    optimalDispatch.gtRampLimit = Constraint(N_gt,T,rule=gtRampLimit)
    optimalDispatch.bolRampLimit = Constraint(N_bol,T,rule=bolRampLimit)
    '''起停状态辅助约束'''
    def gtauxCons1(mdl,n,t):
        if t == 0:
            return mdl.gt_auxvar[n,t] == 0
        else:
            return mdl.gt_auxvar[n,t]>=mdl.gt_state[n,t]-mdl.gt_state[n,t-1]
    def gtauxCons2(mdl,n,t):
        if t == 0:
            return mdl.gt_auxvar[n,t] == 0
        else:
            return mdl.gt_auxvar[n,t]>=-(mdl.gt_state[n,t]-mdl.gt_state[n,t-1])
    def bolauxCons1(mdl,n,t):
        if t == 0:
            return mdl.bol_auxvar[n,t] == 0
        else:
            return mdl.bol_auxvar[n,t]>=mdl.bol_state[n,t]-mdl.bol_state[n,t-1]
    def bolauxCons2(mdl,n,t):
        if t == 0:
            return mdl.bol_auxvar[n, t] == 0
        else:
            return mdl.bol_auxvar[n,t]>=-(mdl.bol_state[n,t]-mdl.bol_state[n,t-1])
    optimalDispatch.gtauxCons1 = Constraint(N_gt,T,rule=gtauxCons1)
    optimalDispatch.gtauxCons2 = Constraint(N_gt, T, rule=gtauxCons2)
    optimalDispatch.bolauxCons1 = Constraint(N_bol,T,rule=bolauxCons1)
    optimalDispatch.bolauxCons2 = Constraint(N_bol, T, rule=bolauxCons2)
    '''Define Objectives'''

    def OM_Cost(mdl):
        om_cost = 0
        for id in microgrid_device.keys():
            if (id in N_es):
                om_cost += microgrid_device[id].om * step * (
                sum(mdl.es_power_out[id, t] for t in T) + sum(mdl.es_power_in[id, t] for t in T))
            if (id in N_cs):
                om_cost += microgrid_device[id].om * step * (
                sum(mdl.cs_cold_out[id, t] for t in T) + sum(mdl.cs_cold_in[id, t] for t in T))
            if (id in N_absc):
                om_cost += microgrid_device[id].om * step * sum(mdl.absc_heat_in[id, t] for t in T)
            if (id in N_bol):
                om_cost += microgrid_device[id].om * step * sum(mdl.bol_power[id, t] for t in T)
            if (id in N_ac):
                om_cost += microgrid_device[id].om * step * sum(mdl.ac_power[id, t] for t in T)
        return om_cost

    def Dep_Cost(mdl):
        return sum(sum(microgrid_device[id].Cbw * step * mdl.es_power_out[id, t] for id in N_es) for t in T)

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
    def StartShutdownFee(mdl):
        return sum(microgrid_device[n].ON_OFF_COST*sum(mdl.gt_auxvar[n,t] for t in T) for n in N_gt) + sum(microgrid_device[n].ON_OFF_COST*sum(mdl.bol_auxvar[n,t] for t in T) for n in N_bol)
    def obj_Economical(mdl):
        return OM_Cost(mdl) + Dep_Cost(mdl) + Fuel_Cost(mdl) + ElectricalFee(mdl) + HeatFee(mdl)+StartShutdownFee(mdl)
    optimalDispatch.obj_Economical = obj_Economical
    optimalDispatch.objective = Objective(rule=obj_Economical)
    '''sub problem'''
    optimalDispatch.sub = SubModel()
    optimalDispatch.sub.wp = Var(T)
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
def extendedResult(result):
    writer = pd.ExcelWriter('extenedResult.xlsx')
    sheet1 = pd.DataFrame()
    sheet1['电负荷'] = -result['交流负荷']-result['直流负荷']
    sheet1['购电功率'] = result['购电功率']
    sheet1['光伏出力'] = result['光伏出力']
    sheet1['电储能放电功率'] = df_sum(result,[col for col in result.columns if '电储能放电功率' in col])
    sheet1['燃气轮机发电功率'] = df_sum(result,[col for col in result.columns if '机组出力' in col])
    sheet1['冰蓄冷耗电功率'] = -df_sum(result,[col for col in result.columns if '冰蓄冷耗电功率' in col])
    sheet1['电储能充电功率'] = -df_sum(result,[col for col in result.columns if '电储能充电功率' in col])
    sheet1['空调制冷耗电功率'] = -df_sum(result,[col for col in result.columns if '空调制冷耗电功率' in col])
    sheet1['电价'] =result['电价']
    sheet1.to_excel(writer,sheet_name='电平衡优化调度结果')
    plt.figure(1)
    plt.rcParams['font.sans-serif'] = ['SimHei']
    load, = plt.plot(-sheet1['电负荷'],linewidth=3.0, linestyle='--', label='电负荷')
    sheet1colors = ['#f4f441','#42f486','#f412ee','#ff8000','#41b8f4','#408080','#7f41f4']
    plt.bar(result.index.values.tolist(),sheet1['购电功率'],color = '#f4f441')
    plt.bar(result.index.values.tolist(),sheet1['燃气轮机发电功率'],bottom=sheet1['购电功率'],color = '#42f486')
    plt.bar(result.index.values.tolist(),sheet1['电储能放电功率'],bottom=sheet1['燃气轮机发电功率']+sheet1['购电功率'],color = '#f442ee')
    plt.bar(result.index.values.tolist(), sheet1['光伏出力'], bottom=sheet1['燃气轮机发电功率'] + sheet1['购电功率']+sheet1['电储能放电功率'],
            color='#ff8000')
    plt.bar(result.index.values.tolist(),sheet1['电储能充电功率'],color = '#41b8f4')
    plt.bar(result.index.values.tolist(),sheet1['冰蓄冷耗电功率'],bottom=sheet1['电储能充电功率'],color = '#408080')
    plt.bar(result.index.values.tolist(),sheet1['空调制冷耗电功率'],bottom=sheet1['冰蓄冷耗电功率']+sheet1['电储能充电功率'],color ='#7f41f4' )
    first_legend = plt.legend([load],('电负荷',))
    ax = plt.gca().add_artist(first_legend)
    plt.legend([mpatches.Patch(color = c) for c in sheet1colors],['购电功率','燃气轮机发电功率','电储能放电功率','光伏出力','电储能充电功率','冰蓄冷耗电功率','空调制冷耗电功率'])
    plt.xlabel('时间')
    plt.ylabel('功率(kW)')
    plt.show()
    '''----------------华丽的分割线--------------------'''
    sheet2 = pd.DataFrame()
    sheet2['空调制冷功率'] = df_sum(result,[col for col in result.columns if '空调制冷功率' in col])
    sheet2['冰蓄冷供冷功率'] = df_sum(result,[col for col in result.columns if '冰蓄冷供冷功率' in col]) + df_sum(result,[col for col in result.columns if '冰蓄冷制冷机直接供冷耗电功率' in col])
    sheet2['吸收式制冷机制冷功率'] = df_sum(result,[col for col in result.columns if '吸收式制冷机制冷功率' in col])
    sheet2['冷负荷'] = -result['冷负荷']
    sheet2['电价'] = result['电价']
    sheet2.to_excel(writer, sheet_name='冷平衡优化调度结果')
    plt.figure(2)
    plt.rcParams['font.sans-serif'] = ['SimHei']
    load, = plt.plot(-sheet2['冷负荷'], linewidth=3.0, linestyle='--', label='冷负荷')
    sheet2colors = ['#f4f441', '#42f486', '#f442ee']
    #plt.stackplot(result.index.values.tolist(), sheet2['空调制冷功率'], sheet2['冰蓄冷供冷功率'], sheet2['吸收式制冷机制冷功率'],colors=sheet2colors)
    plt.bar(result.index.values.tolist(), sheet2['空调制冷功率'], color='#f4f441')
    plt.bar(result.index.values.tolist(), sheet2['冰蓄冷供冷功率'], bottom=sheet2['空调制冷功率'], color='#42f486')
    plt.bar(result.index.values.tolist(), sheet2['吸收式制冷机制冷功率'], bottom=sheet2['冰蓄冷供冷功率'] + sheet2['空调制冷功率'],color='#f442ee')
    first_legend = plt.legend([load], ('冷负荷',))
    ax = plt.gca().add_artist(first_legend)
    plt.legend([mpatches.Patch(color = c) for c in sheet1colors],['空调制冷功率','冰蓄冷供冷功率','吸收式制冷机制冷功率'])
    plt.xlabel('时间')
    plt.ylabel('功率(kW)')
    plt.show()
    '''----------------华丽的分割线--------------------'''
    sheet3 = pd.DataFrame()
    sheet3['冰蓄冷供冷功率'] = df_sum(result, [col for col in result.columns if '冰蓄冷供冷功率' in col])+ df_sum(result,[col for col in result.columns if '冰蓄冷制冷机直接供冷耗电功率' in col])
    sheet3['空调制冷功率'] = df_sum(result, [col for col in result.columns if '空调制冷功率' in col])
    sheet3['吸收式制冷机制冷功率'] = df_sum(result, [col for col in result.columns if '吸收式制冷机制冷功率' in col])
    sheet3['冷负荷'] = result['冷负荷']
    sheet3['电价'] = result['电价']
    sheet3.to_excel(writer, sheet_name='不考虑热品位冷平衡')
    '''----------------华丽的分割线--------------------'''
    sheet4 = pd.DataFrame()
    sheet4['购热功率'] = result['购热功率']
    sheet4['余热锅炉回收热功率'] = df_sum(result,[col for col in result.columns if '余热锅炉' in col])
    sheet4['吸收式制冷机耗热功率'] = df_sum(result, [col for col in result.columns if '吸收式制冷机制冷功率' in col])/0.8
    sheet4['蒸汽驱动负荷'] = -result['蒸汽负荷']
    sheet4['热负荷'] = sheet4['蒸汽驱动负荷']
    sheet4['电价'] = result['电价']
    sheet4.to_excel(writer, sheet_name='不考虑热品位热平衡')
    '''----------------华丽的分割线--------------------'''
    sheet5 = pd.DataFrame()
    sheet5['购热功率'] = result['购热功率']
    sheet5['余热锅炉中品位热功率'] = df_sum(result,[col for col in result.columns if '余热锅炉中品位热功率' in col])
    sheet5['余热锅炉低品位热功率'] =df_sum(result,[col for col in result.columns if '余热锅炉低品位热功率' in col])
    sheet5['蒸汽回收低品位热'] = H2M * result['蒸汽负荷']
    sheet5['中品位热功率'] = result['中品位热功率']
    #sheet5['低品位热功率'] = result['低品位热功率']
    sheet5['蒸汽驱动负荷'] = -result['蒸汽负荷']
    sheet5['热负荷'] = sheet5['蒸汽驱动负荷']
    sheet5['吸收式制冷机耗热功率'] = -df_sum(result, [col for col in result.columns if '吸收式制冷机制冷功率' in col]) / 0.8
    sheet5['吸收式制冷机制冷功率'] = df_sum(result, [col for col in result.columns if '吸收式制冷机制冷功率' in col])
    sheet5.to_excel(writer, sheet_name='考虑热品位热平衡')
    plt.figure(3)
    plt.rcParams['font.sans-serif'] = ['SimHei']
    load, = plt.plot(-sheet5['热负荷'], linewidth=3.0, linestyle='--', label='热负荷')
    sheet5colors = ['#f4f441', '#42f486', '#f442ee']
    #plt.stackplot(result.index.values.tolist(), sheet5['购热功率'], sheet5['余热锅炉中品位热功率'], sheet5['余热锅炉低品位热功率'],sheet5['蒸汽回收低品位热'],sheet5['吸收式制冷机耗热功率'],
                  #colors=sheet5colors)
    plt.bar(result.index.values.tolist(), sheet5['购热功率'], color='#f4f441')
    plt.bar(result.index.values.tolist(), sheet5['余热锅炉中品位热功率'], bottom=sheet5['购热功率'], color='#42f486')
    first_legend = plt.legend([load], ('蒸汽热负荷',))
    ax = plt.gca().add_artist(first_legend)
    plt.legend([mpatches.Patch(color=c) for c in sheet5colors], ['购热功率', '余热锅炉中品位热功率'])
    plt.xlabel('时间')
    plt.ylabel('功率(kW)')
    plt.show()
    '''----------------华丽的分割线--------------------'''
    sheet6 = pd.DataFrame()
    sheet6['蒸汽回收低品位热'] = H2M * result['蒸汽负荷']
    sheet6['余热锅炉低品位热功率'] = df_sum(result, [col for col in result.columns if '余热锅炉低品位热功率' in col])
    sheet6['吸收式制冷机耗热功率'] = -df_sum(result, [col for col in result.columns if '吸收式制冷机制冷功率' in col])/0.8
    sheet6['电价'] = result['电价']
    sheet6.to_excel(writer, sheet_name='吸收式制冷机耗热情况')
    plt.figure(4)
    plt.rcParams['font.sans-serif'] = ['SimHei']
    load, = plt.plot(-sheet6['吸收式制冷机耗热功率'], linewidth=3.0, linestyle='--', label='吸收式制冷机耗热功率')
    sheet6colors = ['#f4f441']
    plt.bar(result.index.values.tolist(), sheet6['蒸汽回收低品位热'], color='#f4f441')
    first_legend = plt.legend([load], ('ABSC耗热功率',))
    ax = plt.gca().add_artist(first_legend)
    plt.legend([mpatches.Patch(color=c) for c in sheet6colors], ['蒸汽回收低品位热'])
    plt.xlabel('时间')
    plt.ylabel('功率(kW)')
    plt.show()
    return
def responseModel(mdl,case,peak,amount,mode):
    model = copy.deepcopy(mdl)
    tmp = ConcreteModel()
    N_T = case.NumOfTime
    T = model.T
    microgrid_data = model.input
    step = 24 / N_T
    k1 = 1
    k2 = 1000
    peak = [t - model.T_range[0] for t in peak]
    model.peak = peak
    model.P_ref = list()
    model.H_ref = list()
    model.P_0 = [value(model.utility_power[t]) for t in T]
    model.H_0 = [value(model.buy_heat[t]) for t in T]
    if mode == 'E':
        for t in T:
            if t in peak:
                model.P_ref.append(value(model.utility_power[t]) - amount[t - peak[0]])
            else:
                model.P_ref.append(8000)
        model.H_ref = [1000]*len(T)
    elif mode == 'H':
        for t in T:
            if t in peak:
                model.H_ref.append(value(model.buy_heat[t]) + amount[t - peak[0]])
            else:
                model.H_ref.append(1000)
        model.P_ref = [8000]*len(T)
    model.DRHeatLoad = Var(peak,bounds=(case.device['DR_Heat_Load'].lower_bound,case.device['DR_Heat_Load'].upper_bound))
    steam_heat_load = microgrid_data['蒸汽负荷'].tolist()
    water_heat_load = microgrid_data['热水负荷'].tolist()
    # '''热损失惩罚函数'''
    # def wastingHeatPenalty(mdl):
    #     return 100000*sum(mdl.medium_heat[t] - steam_heat_load[t] for t in mdl.T)
    # model.wastingHeatPenalty = wastingHeatPenalty
    ''''更新目标函数'''
    def obj_response(mdl):
        if mode == 'E':
            return step * sum((mdl.utility_power[t] - mdl.P_ref[t]) for t in peak)
        elif mode == 'H':
            return step * sum((mdl.H_ref[t] - mdl.buy_heat[t]) for t in peak)
    tmp.obj = Objective(expr=obj_response(model))
    model.objective.set_value(k1*model.objective.expr + k2*tmp.obj.expr)

    ''''更新约束条件'''
    if mode == 'E':
        model.res_curve_u = Constraint(peak, rule=lambda mdl, t: mdl.utility_power[t] - mdl.P_ref[t] >= 0)
        model.pcc_limit = Constraint(set(T) - set(peak), rule=lambda mdl, t: mdl.utility_power[t] <= mdl.P_ref[t])
        model.heat_limit = Constraint(T, rule = lambda mdl,t: mdl.buy_heat[t] >= mdl.H_ref[t])
        model.eq_power = Constraint(peak, rule=lambda mdl,t: (mdl.P_0[peak[0]]-mdl.P_ref[peak[0]])*(mdl.utility_power[t]-mdl.P_0[t]) \
                                                             == (mdl.P_0[t]-mdl.P_ref[t])*(mdl.utility_power[peak[0]]-mdl.P_0[peak[0]]))
    elif mode == 'H':
        model.res_curve_u = Constraint(peak, rule=lambda mdl, t: mdl.buy_heat[t] - mdl.H_ref[t] <= 0)
        model.heat_limit = Constraint(set(T) - set(peak), rule=lambda mdl, t: mdl.buy_heat[t] >= mdl.H_ref[t])
        model.pcc_limit = Constraint(T, rule=lambda mdl, t: mdl.utility_power[t] <= mdl.P_ref[t])
        model.eq_power = Constraint(peak, rule=lambda mdl, t: (mdl.H_0[peak[0]] - mdl.H_ref[peak[0]]) * (
            mdl.buy_heat[t] - mdl.H_0[t]) \
                                                              == (mdl.H_0[t] - mdl.H_ref[t]) * (mdl.buy_heat[peak[0]] - mdl.H_0[peak[0]]))
        del model.HPB2
        del model.HPB2_index
        del model.HPB3
        del model.HPB3_index
        # del model.HPB4
        # del model.HPB4_index
        N_absc = model.case.getKey(absorptionChiller)
        # def low_heat_enough_or_not(b, t, indicator):
        #     m = b.model()
        #     if indicator == 0:  # low heat is not enough
        #         b.low_heat_state = Constraint(expr=m.low_heat[t] <= water_heat_load[t] + sum(m.absc_heat_in[n_absc, t] for n_absc in N_absc))
        #         b.HPB2 = Constraint(T,rule=lambda mdl, t: m.medium_heat[t] + m.low_heat[t]==
        #                                                       steam_heat_load[t] + m.DRHeatLoad[t] + water_heat_load[t] + sum(m.absc_heat_in[n_absc, t] for n_absc in N_absc) if t in peak else
        #         m.medium_heat[t] == steam_heat_load[t] + water_heat_load[t] + sum(m.absc_heat_in[n_absc, t] for n_absc in N_absc))
        #     else:
        #         b.low_heat_state = Constraint(expr=m.low_heat[t] >= water_heat_load[t] + sum(m.absc_heat_in[n_absc, t] for n_absc in N_absc))
        #         b.HPB2 = Constraint(T,rule=lambda mdl, t: m.medium_heat[t] + m.low_heat[t] == steam_heat_load[t] + m.DRHeatLoad[t] if t in peak else
        #         m.medium_heat[t] >= steam_heat_load[t])
        # model.low_heat_enough_or_not = Disjunct(T, [0, 1], rule=low_heat_enough_or_not)
        #
        # def low_heat_enough_disjunct(mdl, t):
        #     return [mdl.low_heat_enough_or_not[t, 0], mdl.low_heat_enough_or_not[t, 1]]
        #
        # model.low_heat_enough_disjunct = Disjunction(T, rule=low_heat_enough_disjunct)
        model.HPB2_1 = Constraint(T,rule = lambda mdl,t:mdl.medium_heat[t] >= steam_heat_load[t] + mdl.DRHeatLoad[t] if t in peak else
                                  mdl.medium_heat[t] >= steam_heat_load[t])
        model.HPB2_2 = Constraint(T,rule = lambda mdl,t:mdl.medium_heat[t] <= steam_heat_load[t] + mdl.DRHeatLoad[t] + sum(mdl.absc_heat_in[n_absc, t] for n_absc in N_absc) if t in peak else
                                  Constraint.Skip)
        model.HPB3 = Constraint(T,rule=lambda mdl,t:mdl.low_heat[t] == H2M*(steam_heat_load[t] + mdl.DRHeatLoad[t]) if t in peak else
        model.low_heat[t] == H2M * (steam_heat_load[t]))

    model.mode = mode
    xfrm = TransformationFactory('gdp.chull')
    xfrm.apply_to(model)
    if mode == 'E':
        solver = SolverFactory('glpk')
        solver.solve(model)
    elif mode == 'H':
        try:
            solver = SolverFactory('gurobi')
            solver.solve(model)
        except Exception:
            solver = SolverManagerFactory('neos')
            solver.solve(model, solver=SolverFactory('cplex'))
    return model

def getMaxAmount(mdl,case,peak,amount,mode):
    model = responseModel(mdl,case,peak,amount,mode)
    if mode == 'E':
        MaxAmount = [model.P_0[t-model.T[0]] - value(model.utility_power[t]) for t in peak]
    elif mode == 'H':
        MaxAmount = [- model.H_0[t-model.T[0]] + value(model.buy_heat[t]) for t in peak]
    else:
        MaxAmount = 0
    return (model,MaxAmount)

def AddDayInSubModel(mdl,t,microgrid_data, case):
    eps = mdl.eps
    N_es = case.getKey(electricStorage)
    N_absc = case.getKey(absorptionChiller)
    N_bol = case.getKey(boiler)
    N_cs = case.getKey(coldStorage)
    N_ac = case.getKey(airConditioner)
    N_gt = case.getKey(gasTurbine)
    N_pv = case.getKey(PV)
    N_ss = list(set(N_bol)|set(N_gt))
    T = mdl.T
    step = mdl.step
    device = case.device
    L = range(4)
    '''MPC Variables'''
    sub = SubModel()
    sub.detP_gt = Var(N_gt,L,bounds=lambda b, i, s: (- device[i].maxDetP, device[i].maxDetP))
    sub.detP_bol = Var(N_gt,L,bounds=lambda b, i, s: (- device[i].maxDetP, device[i].maxDetP))
    sub.detP_ac = Var(N_ac,L,bounds=lambda b, i, s: (- device[i].maxDetP, device[i].maxDetP))
    sub.detP_es = Var(N_es,L)
    sub.detP_cs = Var(N_cs, L)
    sub.detC_cs = Var(N_cs, L)
    sub.es_energy = Var(N_es,L,bounds=lambda b,i,s: (device[i].capacity*device[i].SOCmin,device[i].capacity*device[i].SOCmax))
    sub.cs_cold_stored = Var(N_cs,L,bounds=lambda b,i,s: (device[i].capacity*device[i].Tmin,device[i].capacity*device[i].Tmax))
    sub.detP_ut = Var(L)
    sub.detP_DER = Var(L)
    sub.ESSocAuxVar = Var(N_es,L)
    sub.CSSocAuxVar = Var(N_cs, L)
    '''cold storage constraints'''
    sub.detC_csLimit1 = Constraint(N_cs,L,rule = lambda b,i,s: 0 <= mdl.cs_cold_out[i,t+s] + sub.detC_cs[i,s] <= 2 * mdl.cs_cold_out[i,t+s] + eps)
    sub.detC_csLimit1 = Constraint(N_cs,L,rule = lambda b,i,s: mdl.cs_cold_out[i,t+s] + sub.detC_cs[i,s] <= )

    '''SOC conditions'''
    prev = getattr(mdl.sub,'MPC_'+str(t-1))
    def batteryEnergyBalance(b, n_es, s):
        bat = device[n_es]
        if s == 0:
            if t == 0:
                return -eps <= b.es_energy[n_es, s] - bat.SOCnow * bat.capacity <= eps
            else:
                return -eps <= b.es_energy[n_es, s] - prev.es_energy[n_es,1] <= eps
        else:
            return -eps <= mdl.es_energy[n_es, s] - mdl.es_energy[n_es, s - 1] * (1 - bat.selfRelease) \
                                             - step * mdl.es_power[n_es, s - 1] <= eps

    sub.batteryEnergyBalance = Constraint(N_es, L, rule=batteryEnergyBalance)
    def coldPowerBalance(b, n_cs, s):
        cs = device[n_cs]
        if s == 0:
            if t == 0:
                return -eps <= b.cs_cold_stored[n_cs, s] - cs.SOCnow * cs.capacity <= eps
            else:
                return -eps <= b.cs_cold_stored[n_cs, s] - prev.cs_cold_stored[n_cs,1] <= eps
        else:
            return -eps <= b.cs_cold_stored[n_cs, s] - b.cs_cold_stored[n_cs, s - 1] * (1 - cs.selfRelease) \
                                             - step * b.cs_power[n_cs, s - 1] <= eps
    sub.coldPowerBalance = Constraint(N_cs, L, rule=coldPowerBalance)

    def obj_MinSocError(mdl):
        obj = 0
        for i in N_es:
            obj += sum(mdl.ESSocAuxVar[i,t] for t in L)
        for i in N_cs:
            obj += sum(mdl.CSSocAuxVar[i, t] for t in L)
        return obj
    sub.obj_MinSocError = obj_MinSocError
    sub.objective = Objective(rule=lambda mdl: obj_MinSocError(mdl))
    sub.ESCsocAuxCons1 = Constraint(N_es,L,rule = lambda b,i,s:b.ESSocAuxVar[i,s]>=
                                                                 (b.es_energy[i,s] - mdl.es_energy[i,t+s])/device[i].capacity)
    sub.ESsocAuxCons2 = Constraint(N_es, L, rule=lambda b, i,s: b.ESSocAuxVar[i, s] >=
                                                                - (b.es_energy[i, s] - mdl.es_energy[i, t + s]) / case.device[i].capacity)
    sub.CSCsocAuxCons1 = Constraint(N_cs,L,rule = lambda b,i,s:b.CSSocAuxVar[i,s]>=
                                                                 (b.cs_cold_stored[i,s] - mdl.cs_cold_stored[i,t+s])/device[i].capacity)
    sub.CSsocAuxCons2 = Constraint(N_cs, L, rule=lambda b, i,s: b.CSSocAuxVar[i, s] >=
                                                                - (b.cs_cold_stored[i, s] - mdl.cs_cold_stored[i, t + s]) / case.device[i].capacity)

    setattr(mdl.sub,'MPC_'+str(t),sub)

def df_sum(df,cols):
    newdf = pd.Series([0]*df.__len__(),index=df[cols[0]].index)
    for col in cols:
        newdf = newdf + df[col]
    return newdf
