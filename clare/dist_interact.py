# Suppress warnings
import warnings 
warnings.simplefilter('ignore')

import pyrcel as pm
import numpy as np
import calculations as clc
import aerosol_plots as aerplt
import matplotlib.pyplot as plt

# environmental variables
P0 = 1e5 # Initial Pressure, Pa
T0 = 280.   # Initial Temperature, K
S0 = -0.05  # Initial Supersaturation, 1-RH (95% here)

# regime-determining variables
w = 1.5 # updraft velocity (m/s)
N = 3e3 # total particle number

# initial lognormal distribution variables
# aer1 = ammonium sulfate aerosol (Chen, 2018)
name1 = 'sulfate'
mu1 = 0.06 # mean (um)
sig = 1.5 # geometric standard deviation
kappa1 = 0.61 # hygroscopicity
bins = 100
sulf = pm.AerosolSpecies(name1, pm.Lognorm(mu=mu1, sigma=sig, N=N), kappa=kappa1, bins=bins)

# aer2 = sea salt aerosol (Zieger, 2017)
name2 = 'sea salt'
mu2 = 0.2 # mean (um)
kappa2 = 1.1 # hygroscopicity
ss = pm.AerosolSpecies(name2, pm.Lognorm(mu=mu2, sigma=sig, N=N), kappa=kappa2, bins=bins)

initial_aerosols = [sulf, ss]

dt = 1.0 # timestep (s)
h_end = 10 # end altitude (km)
t_end = h_end*1e3/w # end time (s)

model = pm.ParcelModel(initial_aerosols, w, T0, S0, P0, console=False, accom=0.5)
parcel_trace, aerosol_traces = model.run(t_end, dt, solver='cvode')
sulf_arr = aerosol_traces[name1].values
ss_arr = aerosol_traces[name2].values

aerplt.dist_int(w, N, [sulf.Nis,ss.Nis], [sulf_arr,ss_arr], [name1,name2])