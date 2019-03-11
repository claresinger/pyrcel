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

# initial lognormal distribution variables
mu1 = 0.06 # mean (um)
sig1 = 1.5 # geometric standard deviation
kappa1 = 0.61 # hygroscopicity
N1 = 1e3 # total particle number
sulf = pm.AerosolSpecies('sulfate', pm.Lognorm(mu=mu1, sigma=sig1, N=N1), kappa=kappa1, bins=100)

mu2 = 0.6 # mean (um)
sig2 = 1.5 # geometric standard deviation
kappa2 = 0.1 # hygroscopicity
N2 = 1e3 # total particle number
ss = pm.AerosolSpecies('sea salt', pm.Lognorm(mu=mu2, sigma=sig2, N=N2), kappa=kappa2, bins=100)

initial_aerosols = [sulf, ss]

w = 1.5 # updraft velocity (m/s)
dt = 1.0 # timestep (s)
t_end = 1e3/w # end time (s)... 5 km simulation 

model = pm.ParcelModel(initial_aerosols, w, T0, S0, P0, console=False, accom=0.5)
parcel_trace, aerosol_traces = model.run(t_end, dt, solver='cvode')
sulf_arr = aerosol_traces['sulfate'].values
ss_arr = aerosol_traces['sea salt'].values

aerplt.dist_int(w, [N1,N2], [sulf.Nis,ss.Nis], [sulf_arr,ss_arr])