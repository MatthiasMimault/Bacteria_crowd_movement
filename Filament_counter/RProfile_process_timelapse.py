#!/usr/bin/python
# -*- coding: latin-1 -*-
"""	 Main authors: L. Dupuy
	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License along
	with this program; if not, write to the Free Software Foundation, Inc.,
	51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA."""
import os
import os.path
import numpy as np

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import minimize
from matplotlib.colors import LightSource, Normalize
from scipy.ndimage import gaussian_filter
from scipy.signal import medfilt
from scipy.ndimage.filters import minimum_filter
import pickle
from I_to_CFU import get_param

# TODO
folder_base = os.path.dirname(os.path.realpath(__file__))+"\\" #"C:\\WORK\\PROGRAMMING\\SENSOIL\\FILAMENTS_PROFILING\\"
folder = folder_base + "DATA\\"
folder_plot = folder_base + "PLOT\\"
folder_list = os.listdir(folder)

fontsize = 16	

# Scaling
scale = 95.    # (1 mm = 95 pixels).
# Convert I to CFU
coeff = get_param()
a_I = coeff[0]#7E4 
b_I = coeff[1]#-1.5E7
n_sampling = 10.
threshold = np.log10(1680*a_I + b_I)              # tip
threshold_diameter = np.log10(1670*a_I + b_I) 
threshold_plt = np.log10(1700*0.97*a_I + b_I)     # for the plots	I = 1E-07 CFU + 399,75
                                                  # CFU = I 1E7 - 4 E9
threshold_plt2 = np.log10(1900*a_I + b_I)
ltip = 8                   # length of the tip section in pixels
dict_timing = {}


# Output Data structure
V_TIP_PERCOL = []
V_BASE_PERCOL = []
L_TOT_PERCOL = []
I_TOT_PERCOL = []
D_TOT_PERCOL = []

V_TIP_WATER = []
V_BASE_WATER = []
L_TOT_WATER = []
I_TOT_WATER = []
D_TOT_WATER = []

V_TIP_MS = []
V_BASE_MS = []
L_TOT_MS = []
I_TOT_MS = []
D_TOT_MS = []


# Various functions to analyse the data and fit the Gaussian Ridge
def get_threshold(data):
    return threshold_diameter#th    
def diameter_fit(p):
    if len(p) == 3:
        pp = [ p[0] , 1.5*np.sqrt(p[1])/scale , p[2] ]
    else:
        pp = p
    return pp
def build_dict(file):
    f = open(file,'r')
    i=0
    dict = {}
    for line in f:
        i+=1
        if i>1:
            row = line.split(";")
            id = row[0] + "_" + row[1]
            dict[id] = float(row[5].replace(',','.'))
    return dict        
def RBF(p, X):
	c1, s1, m1, c2 = p
	Y2 = c1*np.exp(-((X-m1)**2)/s1) + c2 
	return Y2
def Root_Function(p,L,X):
    a1, a2, a3 = p

    Y2 = a1*np.exp(-((X)**2)/(a2)) + a3 #*np.exp(b1*L)
    return Y2
def profile_error(p, data):
	X = np.arange(len(data)) - len(data)/2
	return np.sum((RBF(p,X) - data)**2)

def root_error(p, data):
	L,X,D = data
	return np.sum((Root_Function(p,L,X) - D)**2)
	
def optimise_profile(data,x0):
	res = minimize(profile_error, x0, method='Nelder-Mead', tol=1e-4, args = (data,))
	return res.x

def optimise_root(data,x0):
	res = minimize(root_error, x0, method='Nelder-Mead', tol=1e-4, args = (data,))
	return res.x	

# Fit each profile independently
def fit_root_step(data):	
	P_sim = []
	for row in data:
		row = np.array(row)
		x0 = np.array([np.max(row),len(row)*3,0,np.mean(0.5*row[:5] + 0.5*row[-5:])])
		p = optimise_profile(np.array(row), x0)
		#print np.abs(x0[:2] - p[:2]), " / ", x0
		if np.sum(np.abs(x0[:2] - p[:2]) > 150*x0[:2])>0:
			P_sim.append(np.array([np.nan]*4))
		else:
			P_sim.append(p)
	return np.array(P_sim)

# Fit the whole root in one step assuming evolution function for parameters
def fit_root_global(D,L):
    #x = np.arange(D[0,0], D[-1,0]+0.00001, (D[-1,0] - D[0,0])/float(len(D[:,0]) - 1) )
    y = ( np.arange(len(D[0,:])) - (len(D[0,:]))/2 ) #/scale
    X, Y = np.meshgrid(L, y)	

    x0 = np.array([0.01,500,threshold_diameter])

    p = optimise_root([X,Y,np.transpose(D)],x0)
    Z2 = Root_Function(p,X,Y)
    return [np.transpose(Z2),p]

# call plot for debubugging	
def plt_root(data):
	D_sim = []
	P_sim = []
	
	P_sim = fit_root(data)
	for p in P_sim:
		n = len(data[0])
		D_sim.append(RBF(p,np.arange( n )- n/2))

	D = np.array(data)
	fig = plt.figure()


	x = np.arange(D[0,0], D[-1,0]+0.00001, (D[-1,0] - D[0,0])/float(len(D[:,0]) - 1) )
	y = np.arange(len(D[0,:]) - 1) - (len(D[0,:]) - 1)/2
	X, Y = np.meshgrid(x, y)
	Z = np.transpose(D[:,1:])# - np.mean(D[0,1:])

	cmap = cm.coolwarm#plt.cm.copper

	plt.imshow(Z, interpolation='bilinear',
                origin='lower', extent=[x[0], x[-1], y[0], y[-1]],
                vmax=abs(Z).max(), vmin=-abs(Z).max())#, cmap=cm.coolwarm)

	plt.show()
    
def save_plt_timelapse(data, data_centre, data_tip, data_fit, DIAM, title, dir): 
    n = len(data)

    directory = dir + "\\res\\"# + title + "_res\\"
    if not os.path.exists(directory):
        os.makedirs(directory)
        
    directory = directory + title + "\\"
    if not os.path.exists(directory):
        os.makedirs(directory)

    n = len(data)
    m = 6           # approx number of steps       
    inc = max(1,int(n/m+0.5))
    nsub = len(np.arange(0,n,inc))        
    fig = plt.figure(constrained_layout=True, figsize=(5,5)) 
    gs = fig.add_gridspec(1, 2, width_ratios=[10,1], wspace=0.05)

    ii = 0
    
    ###### - Save time lapse of the iamge of the filament
    inc = 1
        
    for i in np.arange(0,n,inc):
        ii+=1
        L = data[i][1]
        diam = np.array(DIAM[i])
        p = np.array(diameter_fit(data_fit[i][2]))*0+0.1
        p = diam
        p = (p<0.5)*p + (p>=5)*0.5
        mtip = len(data_tip[i][1])
        D = data[i][0]
        
        y = np.arange(len(D[0,:])) - (len(D[0,:]))/2

        ax = fig.add_subplot(gs[0, 0])
        tip = data[i][2]
        base = data[i][3]
        Z1 = np.transpose(np.array(data[i][0]))
        map = plt.imshow(Z1, interpolation='bilinear',
                    origin='lower', extent=np.array([L[-1], L[0], -(y[-1]-y[0])*(ii-1) /scale, -(y[-1]-y[0])*(ii) /scale]),cmap=plt.get_cmap('viridis'))#,
        plt.clim(threshold_plt, threshold_plt2)
        ax = plt.gca()
        shift = (y[-1]-y[0])*(-ii+0.5)/scale

        # CODE BELOW ADDS LINES TO INDICATE THE TIP OF THE FILAMENT.
        if tip-mtip>0 and (len(p) >0):
            ax.plot([L[0]-L[tip],L[0]-L[tip]],[-p[1]/2 + shift, p[1]/2  + shift] ,'r', linewidth = 3)
            ax.plot([L[0]-L[tip-mtip],L[0]-L[tip-mtip]],[-p[1]/2 +  shift, p[1]/2 + shift] ,'r:', linewidth = 3)
           
        
        
        for iii in range(2):
            for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize)
            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize)
    # scale bar
    shift += (y[-1]-y[0]) / scale / 2. + 0.1
    ax.plot([abs(L[0]) - 0.3, abs(L[0]) - 1.3] , [shift - 1.8, shift - 1.8], 'w-', linewidth = 3)
    
    ax.axis('off')
    ax.invert_yaxis()
    ax = fig.add_subplot(gs[:, 1])
   
    
    fig.colorbar(map, cax=ax, aspect=50)
    plt.savefig(directory + "RAW" + ".png")
    plt.close()

    ###### - Time lapse of profiles 
    m = 8
    IND = np.arange(0,n,max(int(n/m),1))
    fig = plt.figure(constrained_layout=True, figsize=(5,5)) 
    gs = fig.add_gridspec(7, 1,hspace=1.0)    
    ax1 = fig.add_subplot(gs[:2, 0])
    ax2 = fig.add_subplot(gs[2:, 0])
    for i in IND:
        w = len(data[i][0][:,0])
        h = len(data[i][0][0,:])
        #c = cm.viridis(int(float(i)/float(n)*256.))
        c = cm.viridis(int(float(i+1)/float(n)*256.))
        #nipy_spectral
        #jet
        Y = medfilt(data[i][0], kernel_size=3)
        Y1 = gaussian_filter(Y[:,int(h/2)], sigma = 6)
        ax1.plot(Y1, color = c, linewidth = 2)
        
        Y2 = gaussian_filter(Y[int(w/2),:], sigma = 3)
        if len(Y[0,:])>100:
            Y2 = gaussian_filter(Y[int(w/2),30:-30], sigma = 3)
        ax2.plot(Y2, color = c, linewidth = 2)
        ax=ax1
        for iii in range(2):
            for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize)
            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize)
        ax=ax2
        for iii in range(2):
            for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize)
            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize)
    plt.savefig(directory + "Profile" + ".png")
    plt.close()

    # ###### - IMSHOW - tip and centred data
    # for i in np.arange(0,n):
        # mtip = len(data_tip[i][1])
        # tip = data[i][2]
        # base = data[i][3]
        
        # L = data[i][1]
        # x = np.linspace(L[0], L[-1], len(L))
        # y = np.linspace(-len(D[0,:])/2, len(D[0,:])/2., len(D[0,:]))
        # X,Y = np.meshgrid(x, y)
        # D = data_centre[i][0]
        # y = np.arange(len(D[0,:])) - (len(D[0,:]))/2
        # m = len(D)
        # dl = L[0] - L[1]

        # Z1 = np.transpose(np.array(D))
        # cmap = cm.coolwarm#plt.cm.copper
        # map = plt.imshow(Z1, interpolation='bilinear',
                    # origin='lower', extent=[L[0], L[-1], y[0]/scale, y[-1]/scale])#,
        # plt.colorbar(orientation="horizontal")
        
        # plt.plot([L[tip],L[tip]],[-20/scale,20/scale],'r')

        # if tip-mtip>0:
            # plt.plot([L[tip-mtip],L[tip-mtip]],[-20/scale,20/scale],'r:')
        # plt.plot([L[base],L[base]],[-20/scale,20/scale],'y')
        # map.axes.get_yaxis().set_visible(False)
        # plt.axis('equal')
        # plt.gca().set_aspect('equal', 'box')
        # #plt.clim(threshold_plt, threshold_plt)
        # plt.savefig(directory + "CENTRED_t"+ str(i) + ".png")
        # plt.close()    

    # for i in np.arange(0,n):
        # L = data_fit[i][1]
        # if len(L)>2: 
            # D = data_fit[i][0]
            # y = np.arange(len(D[0,:])) - (len(D[0,:]))/2

            # m = len(D)
            # dl = L[0] - L[1]
            # Z1 = np.transpose(np.array(D))
            # cmap = cm.coolwarm#plt.cm.copper
            # map = plt.imshow(Z1, interpolation='bilinear',
                        # origin='lower', extent=[0, L[-1], y[0]/scale, y[-1]/scale])#,
            # plt.colorbar(orientation="horizontal")
            # map.axes.get_yaxis().set_visible(False)
            # plt.gca().set_aspect('equal', 'box')
            # #plt.clim(threshold_plt, threshold_plt2) 
            # plt.savefig(directory + "TIP_t"+ str(i) + ".png")
            # plt.close()      
    
def plt_timelapse(data, data_centre, data_tip, data_fit, title, dir): 
    n = len(data)
    m = 7           # approx number of steps       
    inc = max(1,int(n/m+0.5))
    nsub = len(np.arange(0,n,inc))
    

    fig = plt.figure(constrained_layout=True) #constrained_layout=True
    
    gs = fig.add_gridspec(nsub, 3, width_ratios=[10,1,10], wspace=0.05)
    
    ii = 0
    for i in np.arange(0,n,inc):
        ii+=1
        L = data[i][1]
        mtip = len(data_tip[i][1])
        D = data[i][0]
        y = np.arange(len(D[0,:])) - (len(D[0,:]))/2
        y = y/scale
        ax = fig.add_subplot(gs[ii-1, 0])
        tip = data[i][2]
        base = data[i][3]
        Z1 = np.transpose(np.array(data[i][0]))

        map = plt.imshow(Z1, interpolation='bilinear',
                    origin='lower', extent=np.array([L[0], L[-1], y[0], y[-1]]),
                    cmap=plt.get_cmap('viridis'))#,
                    #vmax=abs(Z1).max(), vmin=-abs(Z1).max())#, cmap=cm.coolwarm)


        ax.plot([L[tip],L[tip]],[-20/scale,20/scale],'r', linewidth = 2)
        if tip-mtip>0:
            ax.plot([L[tip-mtip],L[tip-mtip]],[-20/scale,20/scale],'r:', linewidth = 2)
        ax.plot([L[base],L[base]],[-20/scale,20/scale],'y', linewidth = 2)
        map.axes.get_yaxis().set_visible(False)
        if ii<nsub:
            map.axes.get_xaxis().set_visible(False)
        
        plt.clim(threshold_plt, threshold_plt2)
        for iii in range(2):
            for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize)
            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize)
    ax = fig.add_subplot(gs[:, 1])
    #ax.
    fig.colorbar(map, cax=ax, aspect=50)
    
    m = 8
    IND = np.arange(0,n,max(int(n/m),1))
    for i in IND:
        ax1 = fig.add_subplot(gs[:int(m/4), 2])
        ax2 = fig.add_subplot(gs[int(m/4)+1:, 2])
        w = len(data[i][0][:,0])
        h = len(data[i][0][0,:])
        #c = cm.viridis(int(float(i)/float(n)*256.))
        c = cm.viridis(int(float(i+1)/float(n)*256.))
        #nipy_spectral
        #jet
        Y = medfilt(data[i][0], kernel_size=3)
        Y1 = gaussian_filter(Y[:,int(h/2)], sigma = 6)
        ax1.plot(Y1, color = c, linewidth = 2)
        Y2 = gaussian_filter(Y[int(w/2),:], sigma = 3)
        ax2.plot(Y2, color = c, linewidth = 2)
        ax=ax1
        for iii in range(2):
            for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize)
            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize)
        ax=ax2
        for iii in range(2):
            for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize)
            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize)
    # fig = plt.figure(2)
    # fig.set_figheight(nsub*10.)
    # plt.title("Centred")
    # ii=0
    # for i in np.arange(0,n,inc):
        # ii+=1
        # mtip = len(data_tip[i][1])
        # tip = data[i][2]
        # base = data[i][3]
        
        # L = data[i][1]
        # x = np.linspace(L[0], L[-1], len(L))
        # y = np.linspace(-len(D[0,:])/2, len(D[0,:])/2., len(D[0,:]))
        # X,Y = np.meshgrid(x, y)
        # D = data_centre[i][0]
        # y = np.arange(len(D[0,:])) - (len(D[0,:]))/2
        # m = len(D)
        # dl = L[0] - L[1]
        # plt.subplot(nsub, 1, ii)
        # Z1 = np.transpose(np.array(D))
        # cmap = cm.coolwarm#plt.cm.copper
        # map = plt.imshow(Z1, interpolation='bilinear',
                    # origin='lower', extent=[L[0], L[-1], y[0], y[-1]])#,
        # #plt.contourf(X,Y,Z1)
        # plt.plot([L[tip],L[tip]],[-20,20],'r')

        # if tip-mtip>0:
            # plt.plot([L[tip-mtip],L[tip-mtip]],[-20,20],'r:')
        # plt.plot([L[base],L[base]],[-20,20],'y')
        # map.axes.get_yaxis().set_visible(False)
        # #plt.axis('equal')
        # plt.clim(threshold_plt, threshold_plt2)
    
    
    
    

     
        
    # fig = plt.figure(3)
    # fig.set_figheight(nsub*10.)
    # plt.title("Fit")
    # ii=0
    # for i in np.arange(0,n,inc):
        # ii+=1
        # L = data_fit[i][1]
        # if len(L)>2:
            
            # D = data_fit[i][0]
            # y = np.arange(len(D[0,:])) - (len(D[0,:]))/2

            # m = len(D)
            # dl = L[0] - L[1]
            # plt.subplot(nsub, 1, ii)
            # Z1 = np.transpose(np.array(D))
            # cmap = cm.coolwarm#plt.cm.copper
            # map = plt.imshow(Z1, interpolation='bilinear',
                        # origin='lower', extent=[0, L[-1], y[0], y[-1]])#,

            # map.axes.get_yaxis().set_visible(False)
            # #plt.axis('equal')
            # plt.clim(threshold_plt, threshold_plt2) 
    # #plt.savefig(dir + "\\" + title + "_dyn_map_tip.png")
              
    plt.show()
    #plt.close()  
    
def export_param(p,file):
	a1, b1, a2, b2, a3, b3 = p

	f = open(file,'w')
	f.write("peak concentration at tip\t" + str(a1)+"\n")
	f.write("decay exudation with L\t" + str(b1)+"\n")
	f.write("spread of exudate\t" + str(a2)+"\n")
	f.write("basal background value\t" + str(a3)+"\n")
	f.write("background decline with L\t" + str(b3)+"\n")		
	f.close()
def export_profile(I,file):

	f = open(file,'w')
	f.write("Integrated intensity profile\n")
	f.write("Distance from the tip , Integrated Intensity\n")
	for i in range(len(I)):
		f.write(str(i) + "," + str(I[i])+"\n")
	f.close()
		
def get_base_tip(data):
    i = len(data)-1
    '''sig = 21#51
    ini = 0
    
    row = medfilt(data[i][ini:],kernel_size=sig)
    density = np.max(row,axis=0)
    while (density < threshold) and i>0:
        i-=1
        row = medfilt(data[i][ini:], kernel_size=sig)
        density = np.max(row,axis=0)
    #plt.plot(row)
    #plt.show()        

    
    j = 0
    density = np.max(medfilt(data[j][ini:],kernel_size=sig),axis=0)
    while (density < threshold) and j<len(data)-2:
        j+=1
        row = medfilt(data[j][ini:], kernel_size=sig)
        density = np.max(row,axis=0)  

    # adjust
    if i<len(data)-1: i+=1        
    if j>0: j-=1'''
    return [len(data)-1,0]

def get_diameter(TIP):
    IJ =[]
    im=0
    jm=0
    for k in range(len(TIP)):
        row =TIP[k,:]
        i = int(len(row)/2)
        j = int(len(row)/2) 
        density = row[i]
        while (density > threshold_diameter) and i>0:
            i-=1
            density = row[i]

        density = row[j]
        while (density > threshold_diameter) and j<len(row)-2:
            j+=1
            density = row[j]
        #print ("     i ", i )
        #print ("     j ", j )
        IJ.append([i,j])
        jm += float(i) / float(len(TIP))
        im += float(j) / float(len(TIP))
    return IJ #[im,jm]
def centre_root(data):
	D2 = []					# Centred data
	I  = []					# integrated data
	B  = []					# background data
	Dtest = []              # jsut some test
	for row in data:
		# fit the profile
		row = np.array(row)[1:]
		row2 = gaussian_filter(row, sigma=5)
		D2.append(row2)
		I.append(np.sum(np.array(row2)-np.mean(0.5*row[:5] + 0.5*row[-5:])))
		B.append([])#p[3])
	# Integrated data
	D2 = np.array(D2)

	return [D2,I, B]    
def centre_root1(data):
	D2 = []					# Centred data
	I  = []					# integrated data
	B  = []					# background data
	Dtest = []              # jsut some test
	for row in data:
		# fit the profile
		row = np.array(row)[1:]
		row2 = gaussian_filter(row, sigma=10)
		#row2 = medfilt(row, kernel_size=11)
		D2.append(row2)
		I.append(np.sum(np.array(row2)-np.mean(0.5*row[:5] + 0.5*row[-5:])))
		B.append([])#p[3])
	# Integrated data
	D2 = np.array(D2)
	D2 = gaussian_filter(D2, sigma = 3)
	return [D2,I, B]
def centre_root2(data):
    D2 = []					# Centred data
    for row in data:
        # fit the profile
        row = np.array(row)[1:]
        row1 = np.concatenate((np.ones(20)*threshold_plt,row,np.ones(20)*threshold_plt),axis=0)

        row2 = gaussian_filter(row1, sigma=5)
        imax = np.max(range(len(row2)) * (row2>= np.max(row2)))
        di = imax - int(len(row2)/2+0.5)
        if di>0:
            row2 = np.concatenate((row1[di:],np.ones(di)*threshold_plt),axis=0)
        elif di<0:
            row2 = np.concatenate((np.ones(-di)*threshold_plt,row1[:di]),axis=0)
        #row2 = gaussian_filter(row2, sigma=10)
        # plt.plot(row2)
        # plt.plot(gaussian_filter(row2, sigma=10))
        # plt.plot(gaussian_filter(row2, sigma=5))
        # plt.show()
        
        D2.append(row2)
        #print (len(D2), "//", len(D2[0]))
        #print(np.shape(np.array(D2)[0]))
    D2 = gaussian_filter(np.array(D2), sigma=1) #[20:-20,:]
    return D2[:,20:-20]
def sort_file_list(file_list):
    file_list2 = {}
    for f in file_list:
        if (f.split("."))[-1] == "txt" and not("Log" in f):
            f1 = f.split("F")
            f2 = f1[1].split(".")
            file_list2[f] = int(f2[0])
    file_dict = {k: v for k, v in sorted(file_list2.items(), key=lambda item: item[1])}
    return file_dict


g=open(folder +"result.csv","w")
g.write("sep=;\n")	
sep=";"
format = "%10.3f"

dict_timing = build_dict(folder + 'metadata.csv')
for ff in folder_list:                       # loop over experiments
    root = folder + ff
    if os.path.isdir(root): 
        n = len(os.listdir(root))
        #**********************************
        g.write(str(ff)+sep+"\n")
        #**********************************
        print (ff)
        for j in range(n):                      # loop over filaments within one image
            fold = root + "\\filament "+ str(j+1) +"\\"
            if os.path.exists(fold):# and "20200510WT_percol" in fold and "filament 2" in fold: 
                
                start_time = dict_timing.get(ff + "_filament " + str(j+1))
                
                #**********************************
                g.write("filament "+ str(j+1) +sep+"\n")
                #**********************************            
                file_list = os.listdir(fold)
                file_list = sort_file_list(file_list).keys()
                D_TOT = []
                D_TIP = []
                D_CENTRE = []
                D_FIT = []
                TIP_tot = []
                DIAM = []
                t=0
                for file in file_list:                  # time lapses of the same filament growing
                    if (file.split("."))[-1] == "txt":
                        t+=1
                        ltip_i = ltip
                        f = open(fold + file,"r")
                        DATA = []
                        
                        i=-1
                        for line in f:
                            i+=1
                            if i>0:
                                row = []
                                cols = line.split("\t")
                                for k in range(len(cols)):
                                    col = cols[k]
                                    if k == 0:                        # distance along the root
                                        row.append(float(col)/scale)
                                    else:                            # Bacteria CFU
                                        if float(col)*a_I + b_I > 0:
                                            row.append(np.log10(float(col)*a_I + b_I))
                                        else:
                                            row.append(np.log10(1630.*a_I + b_I))
                                #if len(row) > 100:
                                #    row[1:] =minimum_filter(np.array(row)[1:], size = 10)
                                DATA.append(row)
                        f.close()

                        #DATA = DATA[::-1]
                        
                        # Varius types of centering
                        D0, I, B = centre_root(DATA)
                        D1, I, B = centre_root1(DATA)
                        D2 = centre_root2(D1)
                        D0 = np.array(DATA)
                        if t==1:
                            th = get_threshold(D1)
                            
                        
                        # Get only the apical part of the filament
                        TP = get_base_tip(D1)#ATA)
                        TIP_tot.append(TP)

                        nbase = max(TP[0]- ltip_i,0)
                        nbase = max(TP[1],nbase)
                        ltip_i = TP[0] - nbase
                        BASE_TIP = D2[nbase:TP[0],:]
                        BASE_TIP = BASE_TIP[::-1,:]

                        DIAM.append(get_diameter(D1[nbase:TP[0],:]))
                        
                        
                        
                        D_TOT.append([D0, np.array(DATA)[:,0], TP[0], TP[1]])
                        
                        Ltip = np.arange(ltip_i)*np.abs(D0[-1,0] - D0[-2,0])  
                        D_TIP.append([BASE_TIP,Ltip])

                        if ltip_i>6:
                            FIT = fit_root_global(BASE_TIP, Ltip)
                            D_FIT.append([FIT[0], Ltip, FIT[1]])
                        else:
                            D_FIT.append([[], [], []])
                        D_CENTRE.append([D2,D0[:,0]])
                        
                    
                #plt_timelapse(D_TOT, D_CENTRE, D_TIP, D_FIT, str(ff)+"_filament "+ str(j+1), root)
                save_plt_timelapse(D_TOT, D_CENTRE, D_TIP, D_FIT, DIAM, str(ff)+"_filament "+ str(j+1), root)
                #**********************************
                # Export data
                
                if start_time != None:
                    T = start_time
                    g.write("time (h)"+sep)
                    for i in range(len(D_TOT)):
                        L = D_TOT[i][1]
                        g.write(format%(T) +sep)
                        T += 0.5
                    g.write(sep+"\n") 
                    
                    
                g.write("tip (mm)"+sep)
                row =[]
                for i in range(len(D_TOT)):
                    L = D_TOT[i][1]
                    g.write(format%(L[D_TOT[i][2]]) +sep)
                    row.append(L[D_TOT[i][2]])
                g.write(sep+"\n") 
                
                
                Vtip = np.gradient(row,0.5)
                g.write("tip velocity (mm/h)"+sep)
                for v in Vtip:
                    g.write(format%(v) +sep)
                g.write(sep+"\n")  
                g.write("max tip velocity (mm/h)"+sep)
                g.write(format%(np.max(-Vtip)) + "\n")
                
                g.write("base (mm)"+sep)
                rowb =[]
                for i in range(len(D_TOT)):
                    Lb = D_TOT[i][1]
                    g.write(format%(Lb[D_TOT[i][3]]) +sep)
                    rowb.append(Lb[D_TOT[i][3]])
                g.write(sep+"\n") 
               

                
                #######################################################################
                Vbase = np.gradient(rowb,0.5)
                g.write("Base velocity (mm/h)"+sep)
                for v in Vbase:
                    g.write(format%(v) +sep)
                g.write(sep+"\n")  
                g.write("max base velocity (mm/h)"+sep)
                g.write(format%(np.max(-Vbase)) + "\n")
                Vbase = np.gradient(rowb,0.5)

                # g.write("Diameter (mm)"+sep)
                
                # Fit
                g.write("Fit parameters"+sep+"\n")

                g.write("Bacterial concentration (CFU / ml)" +sep)
                A1 = []
                for i in range(len(D_FIT)):
                    if len(D_FIT[i][1]) > 0:
                        a1, a2, a3 = diameter_fit(D_FIT[i][2])
                        a1 += a3
                        g.write(format%(a1)+sep)
                        A1.append(np.sum(np.sum(D_TOT[i][0])))

                    else:
                        A1.append(np.sum(np.sum(D_TOT[i][0])))
                        
                g.write(sep+"\n")
                g.write("Median (CFU / ml)" +sep)
                if len(A1) > 2:
                    a1 = np.percentile(A1, 0.5)
                    g.write(format%(a1)+sep)
                else: g.write("not available"+sep)
                g.write(sep+"\n") 
               
                g.write("Diameter (mm)"+sep)
                A2 = []
                A22 = []
                for i in range(len(D_FIT)):
                    if len(D_FIT[i][1]) > 0:
                        a1, a2, a3 = diameter_fit(D_FIT[i][2])  
                        g.write(format%(a2)+sep)
                        A2.append(a2)
                        A22.append(a2)
                    else:
                        A22.append(-1)
                g.write(sep+"\n")
                if len(A2) > 2:
                    g.write("Median (mm)" +sep)
                    a2 = np.percentile(A2, 0.5)
                    g.write(format%(a2)+sep)
                
                else: g.write("not available"+sep)
                g.write(sep+"\n")             
                #********************************** 
                if "water" in ff:
                    V_TIP_WATER.append(Vtip)
                    V_BASE_WATER.append(Vbase)
                    I_TOT_WATER.append(A1)
                    D_TOT_WATER.append(A22)
                    L_TOT_WATER.append(np.array(row) - np.array(rowb))
                elif "mspercol" in ff:
                    V_TIP_MS.append(Vtip)
                    V_BASE_MS.append(Vbase)
                    I_TOT_MS.append(A1)
                    D_TOT_MS.append(A22)
                    L_TOT_MS.append(np.array(row) - np.array(rowb))                        
                else:
                    V_TIP_PERCOL.append(Vtip)
                    V_BASE_PERCOL.append(Vbase)
                    I_TOT_PERCOL.append(A1)
                    D_TOT_PERCOL.append(A22)
                    L_TOT_PERCOL.append(np.array(row) - np.array(rowb))                   
                ############################################
g.close()
pickle.dump( V_TIP_PERCOL, open( folder_plot + "V_TIP_PERCOL.p", "wb" ) )
pickle.dump( V_BASE_PERCOL, open( folder_plot + "V_BASE_PERCOL.p", "wb" ) )
pickle.dump( L_TOT_PERCOL, open( folder_plot + "L_TOT_PERCOL.p", "wb" ) )
pickle.dump( D_TOT_PERCOL, open( folder_plot + "D_TOT_PERCOL.p", "wb" ) )
pickle.dump( I_TOT_PERCOL, open( folder_plot + "I_TOT_PERCOL.p", "wb" ) )

pickle.dump( V_TIP_WATER, open( folder_plot + "V_TIP_WATER.p", "wb" ) )
pickle.dump( V_BASE_WATER, open( folder_plot + "V_BASE_WATER.p", "wb" ) )
pickle.dump( L_TOT_WATER, open( folder_plot + "L_TOT_WATER.p", "wb" ) )
pickle.dump( D_TOT_WATER, open( folder_plot + "D_TOT_WATER.p", "wb" ) )
pickle.dump( I_TOT_WATER, open( folder_plot + "I_TOT_WATER.p", "wb" ) )

pickle.dump( V_TIP_MS, open( folder_plot + "V_TIP_MS.p", "wb" ) )
pickle.dump( V_BASE_MS, open( folder_plot + "V_BASE_MS.p", "wb" ) )
pickle.dump( L_TOT_MS, open( folder_plot + "L_TOT_MS.p", "wb" ) )
pickle.dump( D_TOT_MS, open( folder_plot + "D_TOT_MS.p", "wb" ) )
pickle.dump( I_TOT_MS, open( folder_plot + "I_TOT_MS.p", "wb" ) )

