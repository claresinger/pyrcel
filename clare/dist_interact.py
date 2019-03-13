# Suppress warnings
import warnings 
warnings.simplefilter('ignore')

import pyrcel as pm
import numpy as np
import calculations as clc
import aerosol_plots as aerplt
import matplotlib.pyplot as plt

def main():
    # aerosol type (citation)
    # prop = [name, mean radius (um), hygroscopicity=kappa]

    # ammonium sulfate aerosol (Chen et al., 2018)
    sulf = ['sulfate', 0.06, 0.61]
    # sea salt aerosol (Zieger et al., 2017)
    seasalt = ['sea salt', 0.2, 1.1]
    # mineral dust (Johnson & Osborne, 2011), (Koehler et al., 2009)
    mindust = ['mineral dust',1.3,0.045]
    # black carbon (Wu et al., 2017), (Liu et al., 2013)
    blkcarb = ['black carbon', 0.213, 0.09]
    # fresh biomass burning (Yi, 2018)
    freshburn = ['fresh biomass burning',0.074,0.21]
    
    # # aged biomass burning (Yi, 2018)
    # ageburn = ['aged biomass burning',,0.29]
    # # nitrate
    # nit = ['nitrate']
    # # organic 
    # org = ['misc. organic',,0.15]
    
    # run experiment on known aerosol varieties
    runtest(sulf,seasalt)
    runtest(sulf,mindust)
    runtest(sulf,blkcarb)
    runtest(sulf,freshburn)

    # # test
    # prop2 = ['not sea salt',0.2,0.1]
    # runtest(sulf,prop2)

    mu0 = 0.1 # (um), typical range is 10 - 1000 (nm), so 0.01 - 1.0 (um)
    kappa0 = 0.48 # typical range is 0.15 - 1.2
    prop1 = ['name',mu0,kappa0]
    n = 6
    mus = np.logspace(-2,1,n)
    kappas = 1.5*np.logspace(-2,0,n)
    for i in np.arange(n):
        prop2 = ['name',mus[i],kappas[i]]
        runtest(prop1,prop2)

def runtest(prop1, prop2):
    # environmental variables
    P0 = 1e5 # Initial Pressure, Pa
    T0 = 280.   # Initial Temperature, K
    S0 = -0.15  # Initial Supersaturation, 1-RH (85% here)

    # regime-determining variables
    w = 6.0 # updraft velocity (m/s)
    N = w*2.0e3 # total particle number

    # initial lognormal distribution variables
    # shared variables
    sig = 1.5 # geometric standard deviation
    bins = 200 # number of bins to track (maybe choose this in relation to sigma)

    [name1,mu1,kappa1] = prop1
    aer1 = pm.AerosolSpecies(name1, pm.Lognorm(mu=mu1, sigma=sig, N=N), kappa=kappa1, bins=bins)

    [name2,mu2,kappa2] = prop2
    aer2 = pm.AerosolSpecies(name2, pm.Lognorm(mu=mu2, sigma=sig, N=N), kappa=kappa2, bins=bins)

    initial_aerosols = [aer1, aer2]

    dt = 1.0 # timestep (s)
    h_end = 5e3 # end altitude (m)
    t_end = h_end/w # end time (s)

    model = pm.ParcelModel(initial_aerosols, w, T0, S0, P0, console=False, accom=0.5)
    parcel_trace, aerosol_traces = model.run(t_end, dt, solver='cvode')
    aer1_arr = aerosol_traces[name1].values
    aer2_arr = aerosol_traces[name2].values

    aerplt.dist_int(w, N, [aer1.Nis,aer2.Nis], [aer1_arr,aer2_arr], [name1,name2], [kappa1,kappa2], [mu1,mu2])
    aerplt.dist_compete(w, N, [aer1.Nis,aer2.Nis], [aer1_arr,aer2_arr], [name1,name2], [kappa1,kappa2], [mu1,mu2])


if __name__ == "__main__":
    main()