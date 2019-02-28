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
    hmN = np.ones(len(hmR))*maxN/2
   
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