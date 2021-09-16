# -*- coding: utf-8 -*-
import numpy as np

def quantization(subbands,bpp):    
    A = np.ones(64) #Adjustments
    A[52] = A[56] = 1.32
    A[53] = A[58] = 1.08
    A[54] = A[57] = 1.42
    A[55] = A[59] = 1.08
    
    
    Q = np.empty(64)
    variances = np.empty(60)
    
    for i in range(4):
        Y,X = subbands[i].shape
        x = int(np.floor(3*X/4))        #Width of subregion for variance
        y = int(np.floor(7*Y/16))       #Height of subregion for variance
        
        x_0 = int(np.floor(X/8))
        y_0 = int(np.floor(9*Y/32))
        variances[i] = np.var(subbands[i][y_0:y+y_0-1,x_0:x+x_0-1],ddof=1)              #ddof=1 -> unbiased estimator
    
    Q[0:4] = 1
    for i in range(4,60):
        Y,X = subbands[i].shape
        x = int(np.floor(3*X/4))        #Width of subregion for variance
        y = int(np.floor(7*Y/16))       #Height of subregion for variance
        
        x_0 = int(np.floor(X/8))
        y_0 = int(np.floor(9*Y/32))
        variances[i] = np.var(subbands[i][y_0:y+y_0-1,x_0:x+x_0-1],ddof=1)              #ddof=1 -> unbiased estimator
        
        if(variances[i] < 1.01):
            Q[i] = 0
        else:
            Q[i] = 10/(A[i] * np.log(variances[i]))
    Q[60:] = 0
    """ """
    from constants.dwt import filter_bank_path
    j=0
    K = [k for k in range(64) if Q[k]!=0]

    while(1):
        S = np.sum([1/2**(len(filter_bank_path[k])*2) for k in K])
        
        q = 1/2.5 * 2**(bpp/S-1) 
        z = [(np.sqrt(variances[k])/Q[k])**(1/2**(len(filter_bank_path[k])*2)) for k in K]
        q *= np.prod(z)**(-1/S)
        
        xi = [k for k in K if Q[k]/min(1e300,q) >= 2 * 2.25 * np.sqrt(variances[k])]
        
        if(len(xi) > 0):
            K = [k for k in K if k not in xi]
            j+=1
        else:
            break

    Q = Q/q
    
    Z = 1.2 * Q
    #print(Q)
    
    #Actual quantization
    p = [0]*64
    for i in range(64):
        if(Q[i]!=0):
            p[i] = (subbands[i]>Z[i]/2) * np.floor((subbands[i]-Z[i]/2)/Q[i] + 1) +\
                   (subbands[i]<-Z[i]/2) * np.ceil((subbands[i]+Z[i]/2)/Q[i] - 1)
        else:
            p[i] = np.zeros(subbands[i].shape, dtype=np.intc)
               
    return p, Q, Z, 0.44
               
def iquantization(p ,Q_k, Z, C):

    s_hat = [0]*64
    
    for i in range(64):
        s_hat[i] = (p[i]>0) * ((p[i]-C)*Q_k[i] + Z[i]/2) +\
                   (p[i]<0) * ((p[i]+C)*Q_k[i] - Z[i]/2)
    
    return s_hat
