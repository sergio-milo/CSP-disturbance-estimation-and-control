import random

def noisy_meas(t_out,noise_level):
    
    noisy_meas=[]
    num_measurements = t_out.shape[1]
    num_sensor = t_out.shape[0]
    
    for j in range (0, num_sensor):
        noisy_meas.append([])
        for  i in range (0, num_measurements):
            noisy_meas[j].append(t_out[j][i] + random.uniform(-noise_level, noise_level))
            
    return noisy_meas