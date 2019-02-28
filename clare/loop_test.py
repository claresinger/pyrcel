# Suppress warnings
import warnings 
warnings.simplefilter('ignore')

import pyrcel as pm
import numpy as np
import calculations as clc
import aerosol_plots as aerplt
import matplotlib.pyplot as plt

# environmental variables
P0 = 100000. # Pressure, Pa
T0 = 280.   # Temperature, K
S0 = -0.15  # Supersaturation, 1-RH (85% here)

# initial lognormal distribution variables
mu1 = 0.1 # mean
kappa1 = 0.5 # hygroscopicity

# variables to change to test aerosol- vs updraft-limited regime
N1list = [1e3,1e4,1e5,1e6,1e7] # total number
wlist = [1e0,1e1,1e2] # updraft speed, m/s
# N1list = [1e2,1e5] # total number
# wlist = [1e0] # updraft speed, m/s

ws = list()
N1s = list()
Nis = list()
aer_arrays = list()

print("w","\t","Na")
for i in np.arange(len(N1list)):
    N1 = N1list[i]
    for j in np.arange(len(wlist)):
        w = wlist[j]

        print(w,"\t",N1,"\t",w/N1)

        dt = 1.0 # timestep, seconds
        t_end = 5e3/w # end time, seconds... 5 km simulation 

        name = 'aer1'
        aer1 = pm.AerosolSpecies(name, pm.Lognorm(mu=mu1, sigma=2.0, N=N1), kappa=kappa1, bins=200)
        Ni = aer1.Nis
        initial_aerosols = [aer1]

        model = pm.ParcelModel(initial_aerosols, w, T0, S0, P0, console=False, accom=0.3)
        parcel_trace, aerosol_traces = model.run(t_end, dt, solver='cvode')
        aer_array = aerosol_traces[name].values

        ws.append(w)
        N1s.append(N1)
        Nis.append(Ni)
        aer_arrays.append(aer_array)

        aerplt.eps_ev(w,N1,Ni,aer_array)
        aerplt.dist_ev(w,N1,Ni,aer_array)

aerplt.mult_eps_ev(ws,N1s,Nis,aer_arrays)
aerplt.mult_eps_static(ws,N1s,Nis,aer_arrays)