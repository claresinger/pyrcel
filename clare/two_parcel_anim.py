# Suppress warnings
import warnings 
warnings.simplefilter('ignore')

import pyrcel as pm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.animation as animation

def main():
    h = 1e3

    w = 6.0 # updraft velocity (m/s)
    #N = w*2.0e1 # total particle number, aerosol-limited
    #savestr = "aerosol-limited"
    #N = w*2.0e3 # total particle number, transition
    #savestr = "transition"
    N = w*2.0e5 # total particle number, updraft-limited
    savestr = "updraft-limited"
    ratio = w/N
    initial_aerosols, parcel_trace, aerosol_traces = run_two_component_parcel(h, N)
    
    keys = list(aerosol_traces.keys())
    make_anim(parcel_trace, aerosol_traces, initial_aerosols, keys, ratio, savestr)

def make_anim(parcel_trace, aerosol_traces, initial_aerosols, keys, ratio, savestr):
    T = parcel_trace["T"].values
    z = parcel_trace["z"].values
    N1 = initial_aerosols[0].Nis
    maxN1 = np.max(N1)
    r1 = aerosol_traces[keys[0]].values
    N2 = initial_aerosols[1].Nis
    maxN2 = np.max(N2)
    r2 = aerosol_traces[keys[1]].values

    frames = np.shape(z)[0]

    # initialize figure and axes
    fig = plt.figure(figsize=(10,6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1,4]) 
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    fig.suptitle('$W/N_a = {:.1e}$'.format(ratio), y=0.95, fontsize=14)
    
    # initialize parcel profile plot
    Tprof, = ax1.plot(T[0:1],z[0:1],'ro-')
    ax1.grid(True, linestyle = '-', color = '0.5')
    ax1.set_xlim([np.min(T), np.max(T)])
    ax1.set_ylim([np.min(z), np.max(z)])
    ax1.set_xlabel("Temp. (K)")
    ax1.set_ylabel("Alt. (m)")
    
    # initialize aerosol distribution plot
    dist1, = ax2.semilogx(r1[0,:], N1/maxN1, 'bo-',label=keys[0])
    dist2, = ax2.semilogx(r2[0,:],N2/maxN2, 'ko-',label=keys[1])
    ax2.grid(True, linestyle = '-', color = '0.55')
    ax2.set_xlim([np.min([r1,r2]), np.max([r1,r2])])
    ax2.set_ylim([0,1])
    ax2.set_xlabel(r"Droplet radius ($\mu$m)")
    ax2.set_ylabel("Number")
    ax2.legend(loc=2)
    ax2.yaxis.tick_right()

    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=30)

    anim = animation.FuncAnimation(fig, _update_plot, fargs = (fig, Tprof, T, z, dist1, dist2, r1, N1, maxN1, r2, N2, maxN2), frames=frames)
    hstr = "{:.1f}".format(np.max(z)/1e3)
    anim.save("anims/"+savestr+"_two_component_parcel_"+hstr+".mp4", writer=writer, dpi=300)

def _update_plot(i, fig, Tprof, T, z, dist1, dist2, r1, N1, maxN1, r2, N2, maxN2):
    # update parcel plot
    Tprof.set_xdata(T[0:i+1])
    Tprof.set_ydata(z[0:i+1])
    # update aerosol distribution plot
    dist1.set_xdata(r1[i,:])
    dist1.set_ydata(N1/maxN1)
    dist2.set_xdata(r2[i,:])
    dist2.set_ydata(N2/maxN2)
    # return
    return Tprof, dist1, dist2

def run_two_component_parcel(h, N):
    # environmental variables
    P0 = 1e5 # Initial Pressure, Pa
    T0 = 280.   # Initial Temperature, K
    S0 = -0.15  # Initial Supersaturation, 1-RH (85% here)

    # regime-determining variables
    w = 6.0 # updraft velocity (m/s)

    # initial lognormal distribution variables
    # shared variables
    sig = 1.5 # geometric standard deviation
    bins = 200 # number of bins to track (maybe choose this in relation to sigma)
    
    name1 = "$\kappa = 0.1$"
    mu1 = 0.1
    kappa1 = 0.1
    aer1 = pm.AerosolSpecies(name1, pm.Lognorm(mu=mu1, sigma=sig, N=N), kappa=kappa1, bins=bins)

    name2 = "$\kappa = 1.0$"
    mu2 = 0.1
    kappa2 = 1.0
    aer2 = pm.AerosolSpecies(name2, pm.Lognorm(mu=mu2, sigma=sig, N=N), kappa=kappa2, bins=bins)
    
    initial_aerosols = [aer1,aer2]

    dt = 1.0 # timestep (s)
    h_end = h # end altitude (m)
    t_end = h_end/w # end time (s)

    model = pm.ParcelModel(initial_aerosols, w, T0, S0, P0, console=False, accom=0.5)
    parcel_trace, aerosol_traces = model.run(t_end, dt, solver='cvode')

    return(initial_aerosols, parcel_trace, aerosol_traces)

if __name__ == "__main__":
    main()