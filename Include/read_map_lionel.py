import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt
import numpy.ma as ma
from scipy.ndimage import gaussian_filter
from scipy import ndimage
import scipy.signal as sp
import os
################################################################
file = ""
file_ini = ""
fold = ""
##############################################################################
fontsize = 16

# tracing parameters
thr = 50           #threshold below which this is not a filament anymore

def trace(MAT, ini):
    """
    Trace the trajectories of filaments
    """
    n = len(MAT[:,0])
    MX = []
    L = [[],[],[],[]]
    D = []
    I = []
    ntip = 7
    for i in np.arange(ini,n,1):
        row = MAT[i,:]
        row = gaussian_filter(row, sigma = 2)
        mx, properties = sp.find_peaks(row, prominence=1, width=1)

        if len(mx) == 0:
            if i>ini + ntip:
                I.append(np.mean(MAT[i-ntip,MX[-ntip]]))
            else:
                I.append(0)
            break
            
        
        if len(properties["widths"]) >0:
            D.append(np.mean(properties["widths"]))
        MX.append(mx)

    return [MX, D, I]
def meaure(DATA, DIAM, INT): 
    """
    Get Characteristics of the fialments. Length, number. Diameter not working yet
    """
    N_TOT = []
    L_TOT = []
    D_TOT = []
    I_TOT = []
    for FIL in DATA:
        L = 0
        N = np.zeros(5)
        for i in range(len(FIL)):
            n = len(FIL[i])-1
            if n > -1:
                N[n] += 1
        N_TOT.append(np.sum(N  > 1))
        
        for i in range(len(FIL)-1):
            for ind1 in FIL[i+1]:
                for ind0 in FIL[i]:
                    if 3 > np.abs(ind0 - ind1):
                        L += np.sqrt(1 + (float(ind0 - ind1))**2)
                        break
        L_TOT.append(L)               
    for D in DIAM:
        D_TOT.append([np.mean(D), np.std(D)])
    for I in INT:
        I_TOT.append(np.mean(I))
    return [L_TOT, N_TOT, D_TOT, I_TOT]
def plot_step(DATA,DIAM,INT,X,Y,BB, file):
    FIL = DATA[-1]
    fig = plt.figure(constrained_layout=True, figsize=(5,7))
    gs = fig.add_gridspec(5,1, height_ratios=[4,1,1,1,1])#,wspace=0.7,hspace=0.7) #  width_ratios=[10,1],
    
    #########################################################
    # Density with tracing
    #ax.set_facecolor((0.5, 0., 0.))
    ext = [X[0,0], X[0,-1], Y[0,0], Y[-1,0]]
    ax = fig.add_subplot(gs[0, 0])
    map = ax.imshow(BB, interpolation='bilinear',
                        origin='lower', extent=ext)#,
    map.set_clim(vmin=thr, vmax=thr+100)

    for i in np.arange(0,len(FIL),4):
        row = FIL[i]
        ii = i + ini
        for p in row:
            ax.plot(X[ii,p],Y[ii,p], ".r", linewidth = 0.1, markersize = 0.5)

    
    ax.set_xlim(-3,3)
    ax.set_ylim(-1,3.8)
    for iii in range(2):
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
    
    #####################################################
    # Counts
    L,N,D,I = meaure(DATA, DIAM, INT)
    ax = fig.add_subplot(gs[1, 0])
    ax.plot(L, linewidth = 2)
    ax.set_ylabel("Length")
    ax.set_xlim(0,16)
    plt.locator_params(axis='y', nbins=3)
    for iii in range(2):
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)    
    
    ax = fig.add_subplot(gs[2, 0])
    ax.plot(N, linewidth = 2)
    ax.set_ylabel("Number")
    ax.set_xlim(0,16)
    plt.locator_params(axis='y', nbins=3)
    for iii in range(2):
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)

    #####################################################
    # Diameter
    D = np.array(D) 
    ax = fig.add_subplot(gs[3, 0])
    #ax.plot(range(len(D[:,0])), D[:,0], linewidth = 2)
    ax.errorbar(range(len(D[:,0])), D[:,0],yerr=D[:,1], linewidth = 2)
    ax.set_ylabel("Diameter")
    ax.set_xlim(0,16)
    plt.locator_params(axis='y', nbins=3)
    for iii in range(2):
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
    #####################################################
    # Density
    I = np.array(I) 
    ax = fig.add_subplot(gs[4, 0])
    ax.plot(I, linewidth = 2)
    ax.set_ylabel("Intensity tip")
    ax.set_xlim(0,16)
    plt.locator_params(axis='y', nbins=3)
    for iii in range(2):
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
    plt.savefig(fold + file+"_plot" + ".png")
    plt.savefig("../Png-"+nameCase+"/Filament"+file[:-4]+".png")
##################################
# List files
##################################

#Init file
files = os.listdir()
for f in files:
    if "ini" in f:
        file_ini = f
        break
sim0 = loadmat(fold+file_ini)
nameCase = str(sim0["caseName"][0])
nx = int(sim0["nx"][0][0])
ini = int(nx*0.27)  #Where tracing starts
DOM = np.array(sim0["domDef"])        
DOM = DOM[::-1,:]
##############################
# Mask
DOM = np.array(1-DOM.reshape((nx,nx))).astype('int')
DOM2 = ndimage.binary_erosion(DOM).astype('int')
DOM = DOM - DOM2
mask = ma.masked_where(DOM>0, DOM)

##############################
# Read and process other files
FIL = []
DIAM = []
I = []
for f in files:
    
    if ( "ini" in f ) or not(".mat" in f):
        pass
    elif f.endswith(".mat"):
        file = f
        ##################################################
        #read file
        sim = loadmat(fold+file)
        
        B = np.array(sim["b"])
        X = np.array(sim["X"])
        Y = np.array(sim["Y"])
        B = B.reshape((nx,nx))
        B = B[::-1,:]
        Y = Y[::-1,:]
        
        # Uniform background
        B = (B > thr) * B + (1-(B > thr)) * thr

        #################################################
        #Domain
        BB = B
        #BB = ma.masked_array(B,mask)


        #################################################
        #Trace Filaments
        DATA = trace(B,ini)
        FIL.append(DATA[0])
        DIAM.append(DATA[1])
        I.append(DATA[2])
        plot_step(FIL,DIAM,I, X,Y,BB, file)
        ######################################################


