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
    h = 2e3
    initial_aerosols, parcel_trace, aerosol_traces = run_default_parcel(h)
    
    #for i,key in enumerate(aerosol_traces.keys()):
    i = 0
    key = "default"
    init_aer = initial_aerosols[i]
    aerosol_trace = aerosol_traces[key]
    make_anim(parcel_trace, aerosol_trace, init_aer, key)

def make_anim(parcel_trace, aerosol_trace, init_aer, name):
    T = parcel_trace["T"].values
    z = parcel_trace["z"].values
    N = init_aer.Nis
    r = aerosol_trace.values
    frames = np.shape(z)[0]

    # initialize figure and axes
    fig = plt.figure(figsize=(10,6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1,4]) 
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    
    # initialize parcel profile plot
    Tprof, = ax1.plot(T[0:1],z[0:1],'ro-')
    ax1.grid(True, linestyle = '-', color = '0.5')
    ax1.set_xlim([np.min(T), np.max(T)])
    ax1.set_ylim([np.min(z), np.max(z)])
    ax1.set_xlabel("Temp. (K)")
    ax1.set_ylabel("Alt. (m)")
    
    # initialize aerosol distribution plot
    maxN = np.max(N)
    dist, = ax2.semilogx(r[0,:], N/maxN, 'bo-')
    ax2.grid(True, linestyle = '-', color = '0.55')
    ax2.set_xlim([np.min(r), np.max(r)])
    ax2.set_ylim([0,1])
    ax2.set_xlabel(r"Droplet radius ($\mu$m)")
    ax2.set_ylabel("Number")
    ax2.yaxis.tick_right()

    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=30)

    anim = animation.FuncAnimation(fig, _update_plot, fargs = (fig, Tprof, dist, T, z, r, N, maxN), frames=frames)
    hstr = "{:.1f}".format(np.max(z)/1e3)
    anim.save("anims/"+name+"_parcel_"+hstr+".mp4", writer=writer, dpi=300)

def _update_plot(i, fig, Tprof, dist, T, z, r, N, maxN):
    # update parcel plot
    Tprof.set_xdata(T[0:i+1])
    Tprof.set_ydata(z[0:i+1])
    # update aerosol distribution plot
    dist.set_xdata(r[i,:])
    dist.set_ydata(N/maxN)
    # return
    return Tprof, dist,

def run_default_parcel(h):
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
    
    name1 = 'default'
    mu1 = 0.1
    kappa1 = 0.5
    aer1 = pm.AerosolSpecies(name1, pm.Lognorm(mu=mu1, sigma=sig, N=N), kappa=kappa1, bins=bins)
    initial_aerosols = [aer1]

    dt = 1.0 # timestep (s)
    h_end = h # end altitude (m)
    t_end = h_end/w # end time (s)

    model = pm.ParcelModel(initial_aerosols, w, T0, S0, P0, console=False, accom=0.5)
    parcel_trace, aerosol_traces = model.run(t_end, dt, solver='cvode')

    return(initial_aerosols, parcel_trace, aerosol_traces)

if __name__ == "__main__":
    main()