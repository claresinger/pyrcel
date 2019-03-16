import numpy as np
import matplotlib.pyplot as plt

def calc_eps(r,N):
    # find full-max
    maxN = np.max(N)
    maxi = np.where(N == maxN)
    maxR = r[maxi]
    
    # find half-max
    x = np.abs(N - maxN/2)
    hmi = np.concatenate((np.where(x == sorted(x)[0]), np.where(x == sorted(x)[1]))).flatten()
    hmR = r[hmi]
    #hmN = np.ones(len(hmR))*maxN/2
   
    # don't plot
    # plt.figure(figsize=(8,6))
    # plt.semilogx(r,N,'ko-',alpha=0.2)
    # plt.plot(maxR,maxN,'ro')
    # plt.plot(hmR,hmN,'bo')
    # plt.semilogx(r,N,'k-',alpha=0.2)
    # plt.show()
    
    # calculate mean, fwhm (sig), and eps = sig/mu
    mu = float(maxR)
    fwhm = np.abs(hmR[1]-hmR[0])
    eps = fwhm/mu
    
    return(eps)

def calc_int(rx,ry,t,Ni):
    if rx[t,0] < ry[t,0]:
        r1, r2 = rx, ry 
    else:
        r1, r2 = ry, rx
    intr1 = np.argwhere(np.diff(np.sign(r1[t,:] - np.flip(r2[t,:])))).flatten()[0]
    intr2 = np.argwhere(np.diff(np.sign(np.flip(r1[t,:]) - r2[t,:]))).flatten()[0] + 1
    
    n1 = np.sum(Ni[0][intr1:-1])
    n2 = np.sum(Ni[1][0:intr2])
    
    v1 = np.sum((4*np.pi/3)*1e18*(r1[t,intr1:-1])**3*(Ni[0][intr1:-1]))
    v2 = np.sum((4*np.pi/3)*1e18*(r2[t,0:intr2])**3*(Ni[1][0:intr2]))

    #print(t,n1,n2,v1,v2)

    if rx[t,0] < ry[t,0]:
        x = [intr1,intr2]
    else:
        x = [intr2,intr1]
    return(x)

def calc_tot_vol(r,N):
    vols = (4*np.pi/3)*1e18*r**3*N
    if r.size > N.size:
        v = np.sum(vols,axis=1)
    else:
        v = np.sum(vols)
    return(v)