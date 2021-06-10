import pickle
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline
import centre_rota
import I_to_CFU
type=2

if type == 1:
    V_TIP  = pickle.load( open( "V_TIP_PERCOL.p", "rb" ) )
    V_BASE = pickle.load( open( "V_BASE_PERCOL.p", "rb" ) )
    L_TOT  = pickle.load( open( "L_TOT_PERCOL.p", "rb" ) )
    D_TOT  = pickle.load( open( "D_TOT_PERCOL.p", "rb" ) )
    I_TOT  = pickle.load( open( "I_TOT_PERCOL.p", "rb" ) )
elif type == 2:
    V_TIP  = pickle.load( open( "V_TIP_WATER.p", "rb" ) )
    V_BASE = pickle.load( open( "V_BASE_WATER.p", "rb" ) )
    L_TOT  = pickle.load( open( "L_TOT_WATER.p", "rb" ) )
    D_TOT  = pickle.load( open( "D_TOT_WATER.p", "rb" ) )
    I_TOT  = pickle.load( open( "I_TOT_WATER.p", "rb" ) )

elif type == 3:
    V_TIP  = pickle.load( open( "V_TIP_MS.p", "rb" ) )
    V_BASE = pickle.load( open( "V_BASE_MS.p", "rb" ) )
    L_TOT  = pickle.load( open( "L_TOT_MS.p", "rb" ) )
    D_TOT  = pickle.load( open( "D_TOT_MS.p", "rb" ) )
    I_TOT  = pickle.load( open( "I_TOT_MS.p", "rb" ) )



prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

fontsize = 16	
n_ext = 3
lw=1.5
F_LEN_WATER = []
F_NUM_WATER = []
F_LEN_PERCOL = []
F_NUM_PERCOL = []
F_LEN_MS = []
F_NUM_MS = []
F_TIME_WATER = []
F_TIME_PERCOL = []
F_TIME_MS = []


sep = ";"
# Old total filament dataset
if False:
    f = open("total_length.csv","r")
    i=0
    for line in f:
        i+=1    
        if i > 1:
            row = line.split(sep)
            if "water" in  row[0]:
                F_LEN_WATER.append([])
            if "percol" in  row[0]:
                F_LEN_PERCOL.append([])  
            else:
                F_LEN_PERCOL.append([])

            is_read = 0
            for j in range(len(row)):
                col = row[j]
                if j >0:
                    if float(col) > 0:
                        is_read = 1
                if "water" in  row[0] and is_read == 1:
                    F_LEN_WATER[-1].append(float(col))
                elif "percol" in  row[0] and is_read == 1:
                    F_LEN_PERCOL[-1].append(float(col))
                elif "ms" in  row[0] and is_read == 1:
                    F_LEN_PERCOL[-1].append(float(col))
    f.close()  
    g = open("total_number.csv","r")
    i=0
    for line in g:
        i+=1    
        if i > 1:
            row = line.split(sep)
            if "water" in  row[0]:
                F_NUM_WATER.append([])
            if "percol" in  row[0]:
                F_NUM_PERCOL.append([])            

            is_read = 0
            for j in range(len(row)):
                col = row[j]
                if j >0:
                    if float(col) > 0:
                        is_read = 1
                if "water" in  row[0] and is_read == 1:
                    F_NUM_WATER[-1].append(float(col))
                elif "percol" in  row[0] and is_read == 1:
                    F_NUM_PERCOL[-1].append(float(col))
                elif "ms" in  row[0] and is_read == 1:
                    F_NUM_MS[-1].append(float(col))
    g.close()  
else:
    f = open("filament count all data.csv","r")
    i=0
    fil_name =""
    is_new = False
    for line in f:
        i+=1 
        row = line.split(sep)
        if fil_name!= row[0]:
            fil_name = row[0]
            is_new = True
            i=0
        else:
            is_new = False
        
        D = np.array(row[2:]).astype('float')


        if "water" in  row[0]:
            if i==0:
                F_TIME_WATER.append(D)
            elif i==2:
                F_NUM_WATER.append(D)
            elif i==3:
                F_LEN_WATER.append(D)
        if "percol" in  row[0]:
            if i==0:
                F_TIME_PERCOL.append(D)
            elif i==2:
                F_NUM_PERCOL.append(D)
            elif i==3:
                F_LEN_PERCOL.append(D)
        else:
            if i==0:
                F_TIME_MS.append(D)
            elif i==2:
                F_NUM_MS.append(D)
            elif i==3:
                F_LEN_MS.append(D)
    f.close()

fig = plt.figure(constrained_layout=False, figsize=(5,10)) 
gs = fig.add_gridspec(5,2, height_ratios=[1,1,1.5,1.5,1.5],wspace=0.7,hspace=0.7) #  width_ratios=[10,1],


# TOP CURVES FOR TOTAL LENGTH AND NUMBERS
ax = fig.add_subplot(gs[0, :]) 
if type == 1:
    for i in range(len(F_LEN_PERCOL)):
        Y = F_LEN_PERCOL[i]
        X = F_TIME_PERCOL[i] #np.arange(len(Y)) * 0.5
        ax.plot(X,np.array(Y), "-k",linewidth=lw*0.9)#,color=colors[0])
elif type == 2:        
    for i in range(len(F_LEN_WATER)):  
        Y = F_LEN_WATER[i]
        X = F_TIME_WATER[i]#np.arange(len(Y)) * 0.5
        ax.plot(X,Y, "-k",linewidth=lw*0.9)
elif type == 3:        
    for i in range(len(F_LEN_MS)):  
        Y = F_LEN_MS[i]
        X = F_TIME_MS[i]#np.arange(len(Y)) * 0.5
        ax.plot(X,Y, "-k",linewidth=lw*0.9)
#ax.set_ylim(0,27)
#ax.set_xlim(0, 25)
for iii in range(2):
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        
ax = fig.add_subplot(gs[1, :]) 
if type == 1:
    for i in range(len(F_NUM_PERCOL)):
        Y = F_NUM_PERCOL[i]
        X = F_TIME_PERCOL[i] #np.arange(len(Y)) * 0.5
        ax.plot(X,np.array(Y), "-k",linewidth=lw*0.9)#,color=colors[0])
elif type == 2:        
    for i in range(len(F_NUM_WATER)):  
        Y = F_NUM_WATER[i]
        X = F_TIME_WATER[i]#np.arange(len(Y)) * 0.5
        ax.plot(X,Y, "-k",linewidth=lw*0.9)
elif type == 3:        
    for i in range(len(F_LEN_MS)):  
        Y = F_NUM_MS[i]
        X = F_TIME_MS[i]#np.arange(len(Y)) * 0.5
        ax.plot(X,Y, "-k",linewidth=lw*0.9)
#ax.set_ylim(0,5) 
#ax.set_xlim(0, 25) 
for iii in range(2):
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)


########################################################################################################
# Velocity and Length
ax = fig.add_subplot(gs[2, 0]) 
for i in range(len(V_TIP)):
    Y1 = (-V_TIP[i] + V_BASE[i]).tolist()
    L0 = np.abs(L_TOT[i][0])
    Y1.insert(0,-L0/0.5)
    for k in range(n_ext-1):
        Y1.insert(0,0.)

    
    Y1  = -np.array(Y1)
    Y1  =  np.abs(L_TOT[i])  
    Y1 = gaussian_filter( Y1 , sigma = n_ext*0.2)    
    T = 0.5*np.arange(len(Y1)) #- n_ext*0.5
    ax.plot(T, Y1, "-",linewidth=lw)
    
ax.set_xlim(0, 11)

for iii in range(2):
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)



ax = fig.add_subplot(gs[2, 1])         
for i in range(len(V_TIP)):
    I = I_TOT[i]
    I = gaussian_filter(I, sigma = n_ext*0.2)
    #I = I-I.min()
    I = I_to_CFU.convert(I)#np.log10(
    T = 0.5*np.arange(len(I)) #- n_ext*0.5    
    ax.plot(T, I, "-",linewidth=lw)
    ax.set_yscale('log')
    
ax.set_xlim(0, 11)


for iii in range(2):
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)   

# PHASE LIKE DIAGRAM
ax = fig.add_subplot(gs[3:, :]) 
for i in range(len(V_TIP)):
    np.abs(V_TIP[i] - V_BASE[i])
    Y1 = gaussian_filter(np.array(V_TIP[i] - V_BASE[i]),sigma = 0.5) #[2:]
    Y2 = gaussian_filter(np.array(I_TOT[i]),sigma = 0.7) #[2:]
    
    Y1 = Y1[Y2>7.2]
    Y2 = Y2[Y2>7.2]
    if len(Y1)>2 and np.std(Y1)>0:
        Y1 /= (np.max(Y1) - np.min(Y1))*0.5
        Y2 /= (np.max(Y2) - np.min(Y2))*0.5

        # interpolation for smoothing
        dY1 = np.diff(Y1)
        dY2 = np.diff(Y2)
        L = np.sqrt(np.power(dY1,2) + np.power(dY2,2))
        L = np.cumsum(L)
        L = np.insert(L,0,0)
        
        f1 = interp1d(L, Y1, kind='quadratic')#'cubic')
        f2 = interp1d(L, Y2, kind='quadratic')#'cubic')
        T  = np.arange(L[0], L[-1], (L[-1] - L[0])/50.)
        T  = np.insert(T,-1,L[-1])
        
        YY1 = f1(T)
        YY2 = f2(T)
        O = centre_rota.centre_rota_trajectory(YY1,YY2)
        YY1 -= np.median(O[:,0])
        YY2 -= np.median(O[:,1])

        pl, = ax.plot(YY1, YY2, "-",linewidth=lw)
        dYY1 = np.diff(YY1)
        dYY2 = np.diff(YY2)
        L = np.sqrt(np.power(dYY1,2) + np.power(dYY2,2))*100.
        for j in np.arange(0,len(dYY1),5):
            ax.arrow(YY1[j], YY2[j], dYY1[j]/L[j], dYY2[j]/L[j], head_width=0.05, head_length=0.1, fc=pl.get_color(), ec=pl.get_color())
    else:
        print ("Not included: ", i)

ax.locator_params(nbins=5)
ax.set_xlim(-2, 2)
ax.set_ylim(-1.5, 1.5)
for iii in range(2):
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

# PHASE LIKE DIAGRAM II
plt.figure(2)
ax =plt.gca()
nn = 4
is_cycle = [1,3,5,2,6,9] #
is_multi =[4,7,8, 0,10] #
is_cycle = [0,1,3,5,2,4]
k=-1
for i in is_cycle[:nn]:
    k+=1
    np.abs(V_TIP[i] - V_BASE[i])
    Y1 = gaussian_filter(np.array(V_TIP[i] - V_BASE[i]),sigma = 0.5) #[2:]
    #Y1 = gaussian_filter(np.array(V_TIP[i]),sigma = 1.5) #[2:]
    Y2 = gaussian_filter(np.array(I_TOT[i]),sigma = 0.7) #[2:]
    
    Y1 = Y1[Y2>7.2]
    Y2 = Y2[Y2>7.2]
    if len(Y1)>2 and np.std(Y1)>0:
        Y1 /= (np.max(Y1) - np.min(Y1))*0.5
        Y2 /= (np.max(Y2) - np.min(Y2))*0.5

        # interpolation for smoothing
        dY1 = np.diff(Y1)
        dY2 = np.diff(Y2)
        L = np.sqrt(np.power(dY1,2) + np.power(dY2,2))
        L = np.cumsum(L)
        L = np.insert(L,0,0)
        
        f1 = interp1d(L, Y1, kind='quadratic')#'cubic')
        f2 = interp1d(L, Y2, kind='quadratic')#'cubic')
        #f1 = CubicSpline(L, Y1)#, kind='cubic')
        #f2 = CubicSpline(L, Y2)#, kind='cubic')

        T  = np.arange(L[0], L[-1], (L[-1] - L[0])/50.)
        T  = np.insert(T,-1,L[-1])
        
        YY1 = f1(T)
        YY2 = f2(T)
        O = centre_rota.centre_rota_trajectory(YY1,YY2)
        YY1 -= np.median(O[:,0])
        YY2 -= np.median(O[:,1])

        #mark, = plt.plot(YY1[::5], YY2[::5], "o",linewidth=3)
        pl, = ax.plot(YY1+3*k, YY2, "-",linewidth=lw)
        dYY1 = np.diff(YY1)
        dYY2 = np.diff(YY2)
        L = np.sqrt(np.power(dYY1,2) + np.power(dYY2,2))*100.
        for j in np.arange(0,len(dYY1),5):
            ax.arrow(YY1[j]+3*k, YY2[j], dYY1[j]/L[j], dYY2[j]/L[j], head_width=0.1, head_length=0.15, fc=pl.get_color(), ec=pl.get_color())
    else:
        print ("Not included: ", i)



ax.set_aspect('equal', 'box')
ax.set_xlim(-2, 2*(nn+1)+1)
#ax.set_ylim(-1.5, 1.5)

for iii in range(2):
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

plt.show()       
