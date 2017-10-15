from importFile import *
from microgrid_Model import *
import pandas as pd
def responseModel(model,case,peak,reduceAmount,powerLimit):
    lambda1 = 100
    lambda2 = 1
    T = case.NumOfTime
    load_before = [value(model.utility_power[t]) for t in T]
    responseStrategy = model
    return responseStrategy