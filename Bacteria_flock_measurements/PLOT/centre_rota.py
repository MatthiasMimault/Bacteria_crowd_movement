import numpy as np
from numpy.linalg import*
import pylab as plt
nrm = norm 
def centre_rota_i(v1,v2):
    cross = np.cross(v2,v1)
    dot   = np.dot(v2,v1)
    
    # Distance between centre and M1
    x     = nrm(v1)*dot/nrm(cross) 
    
    # inward normal vector to the trajectory
    
    n     = np.cross(v1,cross)
    n     = n/nrm(n)
    
    return x*n
    
    
def centre_rota_trajectory(X,Y):
    CENTREs = []
    for i in range(len(X)-2):
        x1 = X[i]
        y1 = Y[i]
        x2 = X[i+1]
        y2 = Y[i+1]
        x3 = X[i+2]
        y3 = Y[i+2]  

        v1 = np.array([x2-x1,y2-y1,0])
        v2 = np.array([x3-x2,y3-y2,0])
        
        CENTREs.append(np.array([x1,y1,0])+centre_rota_i(v1,v2))
    return np.array(CENTREs)
        
# X=np.arange(-1,1,0.01)-1
# Y=np.sqrt(1-(X+1)**2) + 5

# O = centre_rota_trajectory(X,Y)
# plt.plot(np.mean(O[:,0]),np.mean(O[:,1]),'o')
# plt.plot(X,Y,'-')
# plt.show()
