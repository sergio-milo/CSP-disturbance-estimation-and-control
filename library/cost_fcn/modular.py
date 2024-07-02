import numpy as np
import matplotlib.pyplot as plt
from pyfmi import load_fmu
import scipy.io
from library.predictor.predictor import Predictor
from library.Noisy_meas.noisy import noisy
from library.Noisy_meas.ThermoCouple import ThermoCouple
from tqdm import tqdm
import time
import random
from scipy.integrate import odeint




class SimulationError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

def Non_negative(x):
    if x>=0:
        return x
    else:
        return 0

def cost_function_tot(mat_file1, mat_file2, mat_file3, fmu_file, tube_fmu, k_p, q_p, k_i, q_i, Time_costant_sensor, noise_level, debug=False):
    cost1 = cost_function(mat_file1, fmu_file, tube_fmu, k_p, q_p, k_i, q_i, Time_costant_sensor, noise_level, debug=False)
    print('First cloud:', cost1)
    cost2 = cost_function(mat_file2, fmu_file, tube_fmu, k_p, q_p, k_i, q_i, Time_costant_sensor, noise_level, debug=False)
    print('Second cloud:', cost2)
    cost3 = cost_function(mat_file3, fmu_file, tube_fmu, k_p, q_p, k_i, q_i, Time_costant_sensor, noise_level, debug=False)
    print('Third cloud:', cost3)
    total_cost= cost1+cost2+cost3
    return total_cost 


def general(mat_file1, mat_file2, mat_file3, fmu_file, tube_fmu, k_p, q_p, k_i, q_i, Time_costant_sensor, noise_level, debug=False):
    try:
        cost = cost_function_tot(mat_file1, mat_file2, mat_file3, fmu_file, tube_fmu, k_p, q_p, k_i, q_i, Time_costant_sensor, noise_level, debug=False)
        return cost
    except Exception as e:
        default_cost=1000
        print(e)
        return default_cost

def cost_function(mat_file, fmu_file, tube_fmu, k_p, q_p, k_i, q_i, Time_costant_sensor, noise_level, debug=False):

    # Load the FMU model
    model = load_fmu(fmu_file)

    
    CSP_tube_model=[]

    # Load the FMU model
    for j in range (0,7):
        CSP_tube_model.append(load_fmu(tube_fmu))

    ##------Getting information from the FMU-----##
    num_tot_volumes = model.get_variable_start('Nw')
    num_panels = model.get_variable_start('Np')
    num_termometers = num_panels + 1
    num_volumes_per_panel = model.get_variable_start('Nnodes') - 1  # 21 nodes for 20 panels

    ##--------Importing Data------##

    tmp_mat_data = scipy.io.loadmat(mat_file)  # The mass flow rate computed in Dymola with the PI controller
    w_ = np.array(tmp_mat_data['flow'])

    flux_table = np.array(tmp_mat_data['heatmap'])
    time_s = np.array(tmp_mat_data['time'])

    tmp_t_out_meas = np.array(tmp_mat_data['temperature'])
    t_out_meas = np.zeros([num_termometers, tmp_t_out_meas.shape[1]])
    for i in range(0, num_termometers):
        t_out_meas[i, :] = tmp_t_out_meas[i * num_volumes_per_panel, :]

    Tstart = float(time_s[0])
    Tend = float(time_s[-1])

    dt = round((Tend - Tstart) / len(time_s), 2)

    
        
    for j in range (0,7):
        CSP_tube_model[j].set('tin', 563.15)
        for n in range(0, num_volumes_per_panel):
            CSP_tube_model[j].set(f"d[{n+1}]", np.mean(flux_table[j*num_volumes_per_panel:(j+1)*num_volumes_per_panel, 0], axis=0))
        CSP_tube_model[j].initialize()

    ##-------------FMU SIMULATION---------------##
    P_gain = 0.0
    I_gain = 0.0
    Time_costant_sensor = 2
    noise_level = 1 

    d_meas = np.zeros([num_tot_volumes, flux_table.shape[1]])
    d_meas_mean=np.zeros([num_panels,flux_table.shape[1]])
    
    d_estim_panels = np.zeros([num_panels, flux_table.shape[1]+1])
    d_estim = np.zeros([7, num_volumes_per_panel, flux_table.shape[1]+1])

    t_out_sim = np.zeros([num_termometers, flux_table.shape[1]])
    first_order_deriv_estim = np.zeros([num_panels, flux_table.shape[1]])
    first_order_deriv_meas = np.zeros([num_panels, flux_table.shape[1]])
    
    t_error = np.zeros([num_panels, flux_table.shape[1]])
    d_error = np.zeros([num_panels, flux_table.shape[1]])

    for i in range(0, num_tot_volumes):
        d_meas[i] = flux_table[i]

    time = Tstart    
   

    simulation_results = []
    time_sim = []

    d_distribution = []

    t_out_tube = np.zeros([num_termometers, flux_table.shape[1]])
    t_out_tube_sim = np.zeros([num_termometers, flux_table.shape[1]])

    for j in range(0, 7):
        for n in range(0, num_volumes_per_panel):
            d_estim[j][n][:-1] = np.mean(d_meas[j*num_volumes_per_panel:(j+1)*num_volumes_per_panel, 0], axis=0)
            
    for j in range(0,7):
        d_meas_mean[j][:]=np.mean(d_meas[j*num_volumes_per_panel:(j+1)*num_volumes_per_panel, :], axis=0)

    P_proportional = P_gain
    I_integrator = I_gain
    correction_rule = []

    for i in range(0, num_panels):
        correction_rule.append(Predictor(prop_gain=P_proportional, int_gain=I_integrator, dt=dt))

    for j in range(0,num_panels):    
        correction_rule[j].set_state(state=np.mean(d_meas[j*num_volumes_per_panel:(j+1)*num_volumes_per_panel, 0], axis=0))

    #progress_bar = tqdm(total=len(time_sim), position=0, ncols=100, desc="Simulation Progress")
    i = 0

    t_out_sim_sens = np.zeros([num_termometers, flux_table.shape[1]])
    t_out_meas_sens = np.zeros([num_termometers, flux_table.shape[1]])
    t_out_meas_noisy = np.zeros([num_termometers, flux_table.shape[1]])
    t_out_tube_sens = np.zeros([num_termometers, flux_table.shape[1]])

    P_value = np.zeros(flux_table.shape[1])
    I_value = np.zeros(flux_table.shape[1])

    Tmin = 563.15

    z0 = np.ones(num_termometers) * Tmin
    p0 = np.ones(num_termometers) * Tmin
    t0 = np.ones(num_termometers) * Tmin

    k_p = k_p
    q_p = q_p

    k_i = k_p
    q_i = q_p

    while time <= Tend:
        P_value[i] = k_p * w_[i, 0] + q_p
        I_value[i] = k_i * w_[i, 0] + q_i

        for j in range(0, num_panels):
            correction_rule[j].update_parameters(P_value[i], I_value[i])

        for j in range(0, num_termometers):
            if time > 0:
                tspan = [time - dt, time]
                p = odeint(ThermoCouple, p0[j], tspan, args=(t_out_meas[j][i], Time_costant_sensor,))
                t_out_meas_sens[j][i] = p[1]
                p0[j] = p[1]

                t_out_meas_noisy[j][i] = noisy(t_out_meas_sens[j][i], noise_level)

        for j in range(0, num_panels):

            CSP_tube_model[j].time = time
            CSP_tube_model[j].set('w', w_[i, 0])
            for n in range(0, num_volumes_per_panel):
                CSP_tube_model[j].set(f"d[{n+1}]", d_estim[j][n][i])

            CSP_tube_model[j].set('tin', t_out_meas_noisy[j,i])
            status4 = CSP_tube_model[j].do_step(time, dt, True)

        if status4 != 0:
            raise ValueError("There are some errors in the Step of the CSP FMI ")

        for j in range(0, 7):
            t_out_sim[j+1][i] = CSP_tube_model[j].get(f"flow2DFV_H.T[21]")
            if time > 0:
                tspan = [time - dt, time]
                t = odeint(ThermoCouple, t0[j], tspan, args=(t_out_sim[j+1][i], Time_costant_sensor,))
                t_out_sim_sens[j+1][i] = t[1]
                t0[j] = t[1]


        t_out_sim_sens[0][i] = t_out_meas[0][i]
        t_out_sim[0][i] = t_out_meas[0][i]
        t_out_tube[0][i] = t_out_meas[0][i]
        t_out_tube_sim[0][i] = t_out_meas[0][i]

        for j in range(0, num_panels):

            if i > 50:
                #t_error[j, i] = (t_out_meas_noisy[j+1][i] - t_out_tube_sens[j+1][i]) - (t_out_sim_sens[j+1, i] - t_out_tube_sens[j+1][i])
                t_error[j, i] = (t_out_meas_noisy[j+1][i]) - (t_out_sim[j+1, i])
                d_estim_panels[j, i+1] = Non_negative(correction_rule[j].step(t_error[j,i]))
                for n in range (0, num_volumes_per_panel): 
                    d_estim[j][n][i+1] = d_estim_panels[j, i+1]

        time_sim.append(time)
        i = i + 1
        time = time + dt
        #progress_bar.update(1)
        #progress_bar.set_postfix({"Percentage": f"{round((time / Tend),3) * 100:.1f}%"}, refresh=True)

    d_estim = d_estim[:][:][:-1]
    d_estim_panels = d_estim_panels[:, :-1]
    
    d_error= d_meas_mean-d_estim_panels 
    cost_d=np.sqrt( np.mean( d_error[:,300:]**2 ))

    cost = np.sqrt( np.mean( t_error[:,300:]**2 ))
    
    #cost = cost + cost_d
    cost = cost_d
    
    if np.isnan(cost) or np.isinf(cost) or cost > 1e6:
        cost = 1e3

    if debug == True:
        return cost, k_integrator, d_meas, d_estim, d_estim_panels, t_out_sim, t_out_meas, w_, num_panels, num_termometers, num_tot_volumes, num_volumes_per_panel, time_sim, first_order_deriv_estim, first_order_deriv_meas
    else:
        return cost
    #progress_bar.close()
