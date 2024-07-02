import random

def noisy(t_out,noise_level):
    noisy=t_out + random.uniform(-noise_level, noise_level)
    return noisy