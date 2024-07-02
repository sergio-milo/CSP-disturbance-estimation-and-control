#!/usr/bin/env python
# coding: utf-8

# In[ ]:



# function that returns dx/dt that then will be integrated 
def ThermoCouple(x,t,u,Time_costant_sensor):
    tau_T=Time_costant_sensor
    dxdt = (-x + u)/tau_T   #Espressione di G(s)=1/(1+sT)
    return dxdt
