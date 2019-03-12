import numpy as np
import matplotlib.pyplot as plt
import calculations as clc

def eps_ev(w,Na,Ni,aer_array):
    endt = np.shape(aer_array)[0]
    eps = np.zeros(endt)
    
    for t in np.arange(endt):
        r = aer_array[t,:]
        eps[t] = clc.calc_eps(r,Ni)
    
    plt.figure(figsize=(10,6))
    fs = 12

    plt.plot(np.arange(endt),eps,'ko-')
    
    title = "$w/N_a = {:.1e}$".format(w/Na)
    plt.title(title,fontsize=fs)
    plt.xlabel("Time, $t$ (s)",fontsize=fs)
    plt.ylabel(r"$\epsilon = \sigma / \mu$",fontsize=1.2*fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.tight_layout()
    
    savename = "./figs/eps_r_{:.1e}_w_{:.1e}_Na_{:.1e}.png".format(w/Na,w,Na)
    plt.savefig(savename,dpi=300)

def mult_eps_ev(ws,Nas,Nis,aer_arrays):
    n = len(ws)
    plt.figure(figsize=(10,6))
    fs = 12

    # loop through w/Na pairs, calculate eps and plot
    for i in np.arange(n):
        w = ws[i]
        Na = Nas[i]
        Ni = Nis[i]
        aer_array = aer_arrays[i]
        
        endt = np.shape(aer_array)[0]
        eps = np.zeros(endt)

        for t in np.arange(endt):
            r = aer_array[t,:]
            eps[t] = clc.calc_eps(r,Ni)

        plt.plot(np.arange(endt),eps,'o-',label="{:.1e}".format(w/Na))
    
    # set plot properties
    plt.xlabel("Time, $t$ (s)",fontsize=fs)
    plt.ylabel(r"$\epsilon = \sigma / \mu$",fontsize=1.2*fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.legend(title="$w/N_a$",fontsize=fs)
    plt.tight_layout()

    # save plot
    savename = "./figs/eps_ev_all.png"
    plt.savefig(savename,dpi=300)

def mult_eps_static(ws,Nas,Nis,aer_arrays):
    n = len(ws)
    plt.figure(figsize=(10,6))
    fs = 12

    ratios = list()
    # loop through w/Na pairs, calculate eps and plot
    for i in np.arange(n):
        w = ws[i]
        Na = Nas[i]
        Ni = Nis[i]
        aer_array = aer_arrays[i]
        
        endt = np.shape(aer_array)[0]
        eps = np.zeros(endt)

        for t in np.arange(endt):
            r = aer_array[t,:]
            eps[t] = clc.calc_eps(r,Ni)

        Deps = eps[-1] - eps[0]
        ratio = w/Na
        ratios.append(ratio)
        
        plt.semilogx(ratio, Deps,'o',label="$w={:.1e}$, $N_a={:.1e}$".format(w,Na))
    
    plt.semilogx([min(ratios),max(ratios)], [0,0], 'k--')

    # set plot properties
    plt.title("Condensational braodening/narrowing of aerosol population",fontsize=fs)
    plt.xlabel("$w/N_a$",fontsize=1.2*fs)
    plt.ylabel(r"$\epsilon = \sigma / \mu$",fontsize=1.2*fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.legend(fontsize=fs)
    plt.tight_layout()

    # save plot
    savename = "./figs/eps_static_all.png"
    plt.savefig(savename,dpi=300)


def dist_ev(w,Na,Ni,aer_array):
    endt = np.shape(aer_array)[0]
    tlist = np.geomspace(1,endt,9,endpoint=True).astype(int)-1

    plt.figure(figsize=(10,6))
    fs = 12

    for t in tlist:
        r = aer_array[t,:]
        plt.semilogx(r*1e6, Ni*1e-6, linestyle='solid', marker='.', label="aerosol, t="+str(t))

    title = "$w/N_a = {:.1e}$".format(w/Na)
    plt.title(title,fontsize=fs)
    plt.xlabel(r"Aerosol wet radius ($\mu$m)",fontsize=fs)
    plt.ylabel(r"Aerosl number conc. (cm$^{-3})$",fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.legend(loc='upper left',fontsize=fs)
    plt.tight_layout()
    
    savename = "./figs/dist_ev_r_{:.1e}_w_{:.1e}_Na_{:.1e}.png".format(w/Na,w,Na)
    plt.savefig(savename,dpi=300)

def dist_int(w,Na,Ni,aer_array,name):
    plt.figure(figsize=(10,6))
    fs = 12
    col = ['blue','limegreen','gold','red','azure','pink']
    
    npop = len(aer_array)
    endt = np.shape(aer_array[0])[0]
    tlist = np.geomspace(1,endt,9,endpoint=True).astype(int)-1

    for t in tlist:
        x = clc.calc_int(aer_array[0],aer_array[1],t,Ni)
        [intr0,intr1] = x
        #print(t,intr0,intr1)

        for i in np.arange(npop):
            N = Ni[i]
            maxN = np.max(N)
            arr = aer_array[i]
            r = arr[t,:]
            
            if t == tlist[-1]:
                plt.semilogx(r*1e6, N/maxN, '-', lw=2, color=col[i], label=name[i])
                if i == 0:
                    plt.semilogx(r[intr0]*1e6, N[intr0]/maxN, 'ro')
                if i == 1:
                    plt.semilogx(r[intr1]*1e6, N[intr1]/maxN, 'ro')
            elif t == tlist[0]:
                plt.semilogx(r*1e6, N/maxN, '-', color=col[i], lw=1)
            else:
                al = 0.2+t/tlist[-1]
                plt.semilogx(r*1e6, N/maxN, '-', color=col[i], alpha=al, lw=0.5)
    
    plt.title("Parcel with two aerosol populations -- $w/N_a = {:.1e}$".format(w/Na),fontsize=fs)
    plt.xlabel(r"Aerosol wet radius ($\mu$m)",fontsize=fs)
    plt.ylabel(r"Aerosol number (normalized)",fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.legend(loc='upper left',fontsize=fs)
    plt.tight_layout()
    
    savename = "sulfate_seasalt.png"
    plt.savefig(savename, dpi=300)

def dist_compete(w,Na,Ni,aer_array,name):
    plt.figure(figsize=(10,6))
    fs = 12
    col = ['blue','limegreen','gold','red','azure','pink']
    
    npop = len(aer_array)
    endt = np.shape(aer_array[0])[0]

    for i in np.arange(npop):
        v = clc.calc_tot_vol(aer_array[i],Ni[i])
        grow = np.diff(v)/v[0:-1]
        plt.loglog(np.arange(1,endt),grow*100.0,'-',color=col[i],label=name[i])

    plt.title("Parcel with two aerosol populations -- $w/N_a = {:.1e}$".format(w/Na),fontsize=fs)
    plt.xlabel(r"Time (s)",fontsize=fs)
    plt.ylabel(r"Percentage change in population volume ($\Delta V/V$)",fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks((0.1,1,10),('0.1%','1%','10%'),fontsize=fs)
    plt.legend(loc='upper left',fontsize=fs)
    plt.tight_layout()
    #plt.xlim((0,200))
    
    savename = "growth_rates.png"
    plt.savefig(savename, dpi=300)