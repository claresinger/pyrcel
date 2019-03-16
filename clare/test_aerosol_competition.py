# Suppress warnings
import warnings
warnings.simplefilter('ignore')

import numpy as np
import pyrcel as pm
import matplotlib.pyplot as plt
import aerosol_plots as aerplt
import calculations as clc
import parcel_anim as panim

def main():
    # max height for simulation
    h = 1e3
    
    # create aerosol list to test on 
    #aerosol_list = get_real_aerosols()
    aerosol_list = generate_aerosols()
    print(aerosol_list)
    
    # create list of all aerosol tests
    init_aerosols_list = list()
    parcel_trace_list = list()
    aerosol_trace_list = list()

    # run parcel for individual aerosols
    for aerosol in aerosol_list:
        initial_aerosols, parcel_trace, aerosol_trace = runparcel([aerosol], h)
        init_aerosols_list.append(initial_aerosols)
        parcel_trace_list.append(parcel_trace)
        aerosol_trace_list.append(aerosol_trace)

#     # run parcel for pairs of aerosols
#     midi = int(len(aerosol_list)/2)
#     print(midi)
#     mid = aerosol_list[int(len(aerosol_list)/2)]
#     print(mid)
#     for i,aerosol in enumerate(aerosol_list):
#         if i != midi:
#             initial_aerosols, parcel_trace, aerosol_trace = runparcel([mid, aerosol], h)
#             init_aerosols_list.append(initial_aerosols)
#             parcel_trace_list.append(parcel_trace)
#             aerosol_trace_list.append(aerosol_trace)
            
    # run experiment on known aerosol varieties
    # test all combinations
    for i in np.arange(len(aerosol_list)):
        for j in np.arange(len(aerosol_list)):
            if j > i:
                initial_aerosols, parcel_trace, aerosol_trace = runparcel([aerosol_list[i], aerosol_list[j]], h)
                init_aerosols_list.append(initial_aerosols)
                parcel_trace_list.append(parcel_trace)
                aerosol_trace_list.append(aerosol_trace)
    
    # i,j = 3,0
    # ptrace = parcel_trace_list[i]
    # atraces = aerosol_trace_list[i]
    # keys = list(atraces.keys())
    # key = keys[j]
    # atrace = aerosol_trace_list[i][key]
    # initaer = init_aerosols_list[i][j]
    # panim.make_anim(ptrace,atrace,initaer,key)
    
def runparcel(aerosols, h):
    npop = len(aerosols)

    # environmental variables
    P0 = 1e5  # Initial Pressure, Pa
    T0 = 280.   # Initial Temperature, K
    S0 = -0.15  # Initial Supersaturation, 1-RH (85% here)

    # regime-determining variables
    w = 6.0  # updraft velocity (m/s)
    Ntot = w*2.0e3 # total particle number
    N = Ntot/npop # particle number per population

    # initial lognormal distribution variables
    # shared variables
    sig = 1.5  # geometric standard deviation
    bins = 200 # number of bins to track (maybe choose this in relation to sigma)

    initial_aerosols = list()
    for i in np.arange(npop):
        prop = aerosols[i]
        [name,mu,kappa] = prop
        aer = pm.AerosolSpecies(name, pm.Lognorm(
            mu=mu, sigma=sig, N=N), kappa=kappa, bins=bins)
        initial_aerosols.append(aer)

    dt = 1.0  # timestep (s)
    h_end = h  # end altitude (m)
    t_end = h_end/w  # end time (s)

    model = pm.ParcelModel(initial_aerosols, w, T0, S0,
                           P0, console=False, accom=0.5)
    parcel_trace, aerosol_traces = model.run(t_end, dt, solver='cvode')

    return(initial_aerosols, parcel_trace, aerosol_traces)

def get_real_aerosols():
    # aerosol type (citation)
    # prop = [name, mean radius (um), hygroscopicity=kappa]

    # ammonium sulfate aerosol (Chen et al., 2018)
    sulf = ['sulfate', 0.06, 0.61]
    # sea salt aerosol (Zieger et al., 2017)
    seasalt = ['sea-salt', 0.2, 1.1]
    # mineral dust (Johnson & Osborne, 2011), (Koehler et al., 2009)
    mindust = ['mineral-dust', 1.3, 0.045]
    # black carbon (Wu et al., 2017), (Liu et al., 2013)
    blkcarb = ['black-carbon', 0.213, 0.09]
    # fresh biomass burning (Yi, 2018)
    freshburn = ['fresh-biomass-burning', 0.074, 0.21]
    # urban (Kreidenweis & Asa-Awuku, 2014).
    urban = ['urban', 0.2, 0.2]
    # aged biomass burning (Malm et al., 2005), (Kreidenweis & Asa-Awuku, 2014) [(Yi, 2018) = 0.29]
    ageburn = ['aged-biomass-burning', 0.33, 0.04]

    # # nitrate
    # nit = ['nitrate',,]
    # # organic (Kreidenweis & Asa-Awuku, 2014).
    # org = ['misc-organic',,0.15]

    aerosol_list = list([sulf, seasalt, mindust, blkcarb,
                        urban, freshburn, ageburn])

    #aerosol_list = list([sulf, seasalt, freshburn, ageburn])

    return(aerosol_list)

def generate_aerosols():
    # aerosol type (citation)
    # prop = [name, mean radius (um), hygroscopicity=kappa]

    mu = 0.1
    aerosol_list = list()
    for kappa in np.linspace(0,1.5,9):
        aer = ['', mu, kappa]
        aerosol_list.append(aer)

    # ammonium sulfate aerosol (Chen et al., 2018)
    sulf = ['sulfate', mu, 0.61]
    # sea salt aerosol (Zieger et al., 2017)
    seasalt = ['sea-salt', mu, 1.1]
    # mineral dust (Johnson & Osborne, 2011), (Koehler et al., 2009)
    mindust = ['mineral-dust', mu, 0.045]
    # black carbon (Wu et al., 2017), (Liu et al., 2013)
    blkcarb = ['black-carbon', mu, 0.09]
    # fresh biomass burning (Yi, 2018)
    bburn = ['biomass-burning', mu, 0.21]

    aerosol_list = list([sulf, seasalt, mindust, blkcarb, bburn])

    return(aerosol_list)

if __name__ == "__main__":
    main()