# Suppress warnings
import warnings
warnings.simplefilter('ignore')

import matplotlib.pyplot as plt
import aerosol_plots as aerplt
import calculations as clc
import numpy as np
import pyrcel as pm

def main():
    run_real_pops()
    # run_ensemble(6)


def run_ensemble(n):
    # aerosol type (citation)
    # prop = [name, mean radius (um), hygroscopicity=kappa]
    mu1 = 0.1  # (um), typical range is 10 - 1000 (nm), so 0.01 - 1.0 (um)
    kappa1 = 0.48  # typical range is 0.15 - 1.2
    prop1 = ['name1', mu1, kappa1]
    mus = np.logspace(-2, 1, n)
    kappas = 1.5*np.logspace(-2, 0, n)
    for i in np.arange(n):
        mu2 = mus[i]
        kappa2 = kappas[i]
        prop2 = ['name2', mu2, kappa2]
        runtest(prop1, prop2)


def run_real_pops():
    # aerosol type (citation)
    # prop = [name, mean radius (um), hygroscopicity=kappa]

    # ammonium sulfate aerosol (Chen et al., 2018)
    sulf = ['sulfate', 0.06, 0.61]
    # sea salt aerosol (Zieger et al., 2017)
    seasalt = ['sea salt', 0.2, 1.1]
    # mineral dust (Johnson & Osborne, 2011), (Koehler et al., 2009)
    mindust = ['mineral dust', 1.3, 0.045]
    # black carbon (Wu et al., 2017), (Liu et al., 2013)
    blkcarb = ['black carbon', 0.213, 0.09]
    # fresh biomass burning (Yi, 2018)
    freshburn = ['fresh biomass burning', 0.074, 0.21]
    # urban (Kreidenweis & Asa-Awuku, 2014).
    urban = ['urban', 0.2, 0.2]
    # aged biomass burning (Malm et al., 2005), (Kreidenweis & Asa-Awuku, 2014) [(Yi, 2018) = 0.29]
    ageburn = ['aged biomass burning', 0.33, 0.04]

    # # nitrate
    # nit = ['nitrate']
    # # organic (Kreidenweis & Asa-Awuku, 2014).
    # org = ['misc. organic',,0.15]

    aerosol_list = list([sulf, seasalt, mindust, blkcarb,
                        urban, freshburn, ageburn])
#     # test each aerosol type seperately
#     for i in np.arange(len(aerosol_list)):
#         pop = aerosol_list[i]
#         runtest_single(pop)

    # run experiment on known aerosol varieties
    # test all combinations
    for i in np.arange(len(aerosol_list)):
        for j in np.arange(len(aerosol_list)):
            if j > i:
                print(aerosol_list[i][0], aerosol_list[j][0])
                print(aerosol_list[i][1], aerosol_list[j][1])
                runtest(aerosol_list[i], aerosol_list[j])


def runtest(prop1, prop2):
    # environmental variables
    P0 = 1e5  # Initial Pressure, Pa
    T0 = 280.   # Initial Temperature, K
    S0 = -0.15  # Initial Supersaturation, 1-RH (85% here)

    # regime-determining variables
    w = 6.0  # updraft velocity (m/s)
    N = w*2.0e3  # total particle number

    # initial lognormal distribution variables
    # shared variables
    sig = 1.5  # geometric standard deviation
    # number of bins to track (maybe choose this in relation to sigma)
    bins = 200

    [name1, mu1, kappa1] = prop1
    aer1 = pm.AerosolSpecies(name1, pm.Lognorm(
        mu=mu1, sigma=sig, N=N), kappa=kappa1, bins=bins)

    [name2, mu2, kappa2] = prop2
    aer2 = pm.AerosolSpecies(name2, pm.Lognorm(
        mu=mu2, sigma=sig, N=N), kappa=kappa2, bins=bins)

    initial_aerosols = [aer1, aer2]

    dt = 1.0  # timestep (s)
    h_end = 5e2  # end altitude (m)
    t_end = h_end/w  # end time (s)

    model = pm.ParcelModel(initial_aerosols, w, T0, S0,
                           P0, console=False, accom=0.5)
    parcel_trace, aerosol_traces = model.run(t_end, dt, solver='cvode')
    aer1_arr = aerosol_traces[name1].values
    aer2_arr = aerosol_traces[name2].values

    aerplt.dist_int(w, N, [aer1.Nis, aer2.Nis], [aer1_arr, aer2_arr], [
                    name1, name2], [kappa1, kappa2], [mu1, mu2])
    aerplt.dist_compete(w, N, [aer1.Nis, aer2.Nis], [aer1_arr, aer2_arr], [
                        name1, name2], [kappa1, kappa2], [mu1, mu2])


def runtest_single(prop1):
    # environmental variables
    P0 = 1e5  # Initial Pressure, Pa
    T0 = 280.   # Initial Temperature, K
    S0 = -0.15  # Initial Supersaturation, 1-RH (85% here)

    # regime-determining variables
    w = 6.0  # updraft velocity (m/s)
    N = w*2.0e3  # total particle number

    # initial lognormal distribution variables
    # shared variables
    sig = 1.5  # geometric standard deviation
    # number of bins to track (maybe choose this in relation to sigma)
    bins = 200

    [name1, mu1, kappa1] = prop1
    aer1 = pm.AerosolSpecies(name1, pm.Lognorm(
        mu=mu1, sigma=sig, N=N), kappa=kappa1, bins=bins)

    initial_aerosols = [aer1]

    dt = 1.0  # timestep (s)
    h_end = 5e2  # end altitude (m)
    t_end = h_end/w  # end time (s)

    model = pm.ParcelModel(initial_aerosols, w, T0, S0,
                           P0, console=False, accom=0.5)
    parcel_trace, aerosol_traces = model.run(t_end, dt, solver='cvode')
    aer1_arr = aerosol_traces[name1].values

    #aerplt.dist_int(w, N, [aer1.Nis], [aer1_arr], [name1], [kappa1], [mu1])
    aerplt.dist_compete(w, N, [aer1.Nis], [aer1_arr], [name1], [kappa1], [mu1])


if __name__ == "__main__":
    main()
