import math
import numpy as np


class Predictor():
    def __init__(self, prop_gain=0, int_gain=0, dt=0):
        self.state = 0
        self.int_gain = int_gain
        self.prop_gain = prop_gain
        self.dt = dt

        
        
        
        
    def update_parameters(self, prop_gain, int_gain):
        self.prop_gain = prop_gain
        self.int_gain = int_gain
    # Integrator discretized with Tustin: out = gain*dt/2*(z+1)/(z-1) * in
    # def step_integrator_tustin(self, input):
    #     integral_action = self.state + 0.5*self.int_gain*self.dt*input
    #     self.state = self.state + self.int_gain*self.dt*input
    #     return integral_action
    
    # Integrator discretized with Forward Euler: out(t+1) = out(t) + gain*dt*in(t)
    def step(self, input):
        output = self.state + self.int_gain*self.dt*input + self.prop_gain*input
        self.state = output
        return output
    
    def set_state(self, state):
        self.state = state