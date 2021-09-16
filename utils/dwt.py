# -*- coding: utf-8 -*-
import numpy as np
def expand(x, L):
    """
    Expands a 1-D array with L-1 zeros
    
    y = (↑x)(n) = x(k/L) if k%L==0; 0 otherwise
    """
    y = np.zeros(x.size*L,dtype=np.float64)
    y[::L] = x
    return y


def expand2D(x, L, axis=0):
    """
    Expands a 2-D array with L-1 zeros on each column/row
    
    x: array_like
        Input data
    L: int
        number of zeros (plus one) between each element
    axis: int, optional
        Selects the direction where the zeros are places. If an axis grater than one is selected, an error is raised.
    
    see expand() for further informations
    """
    
    M,N = x.shape
    
    if(axis==0):
        M = M*L
    elif(axis==1):
        N = N*L
    else:
        raise NotImplementedError()
    
    y = np.zeros((M,N),dtype=x.dtype)
    
    if(axis==0):
        y[::L] = x
    elif(axis==1):
        y[:,::L] = x
    return y


def getSubbandStructure(x, syntesis_function):
    
    def apply_dwt_at_pos(s, pos):
        y = syntesis_function(s)
        if(pos==0):
            return y
        elif(pos==1):
            return y[:2][::-1] + y[2:][::-1]
        elif(pos==2):
            return y[2:] + y[:2]
        elif(pos==3):
            return y[::-1]
        else:
            raise ValueError(f"Cannot apply idwt2 at position {pos}")
            
    subbands = [0]*64
    y = syntesis_function(x)
    for i in range(4):
        y[i] = syntesis_function(y[i])
    
    
    subbands[60:] = y[3][::-1]
    subbands[56:60] = y[2][2:] + y[2][:2]
    subbands[52:56] = y[1][:2][::-1] + y[1][2:][::-1]
    subbands[51] = y[0][3]
    
    y = y[0]        #subbands 0-51
    
    for i in range(3):
        y[i] = apply_dwt_at_pos(y[i],i)
        for j in range(4):
            z = apply_dwt_at_pos(y[i][j],j)
            if(i==0 and j==0):

                subbands[:4] = syntesis_function(z[0])   #subband 0-3
                subbands[4:7] = z[1:]       #subbands 4-6
            else:
                pos = 3+16*i+4*j
                subbands[pos:pos+4] = z

    return subbands


def applyFunctionAtSubbandStructure(s, f):
    subbands = s[:]
        
    layer2 = []
    for i in range(3):
        layer3 = []
        for j in range(4):
            if(i==j==0):
                subbands[3] = f(subbands[:4])
            pos = 3 + 4*j + 16*i
            
            layer3 += [f(subbands[pos:pos+4])]
        layer2 += [f(layer3)]
    
    subbands[48:51] = layer2
        
    layer1 = []
    for i in range(4):
        pos = 48 + 4*i
        layer1 += [f(subbands[pos:pos+4])]

    y = f(layer1)    
    return y

def visualizeSubbandStructure(subbands):
    def merge4(x):
        x1,x2,x3,x4 = x
        return np.concatenate([np.concatenate([x1,x2],axis=1)
                               ,np.concatenate([x3,x4],axis=1)],axis=0)
    return applyFunctionAtSubbandStructure(subbands, merge4)


def getImageFromSubbandStructure(s):
    from SWT import idwt2
    
    def apply_idwt_at_pos(s, pos):
        if(pos==0):
            return idwt2(s)
        elif(pos==1):
            return idwt2(s[:2][::-1] + s[2:4][::-1])
        elif(pos==2):
            return idwt2(s[2:4] + s[:2])
        elif(pos==3):
            return idwt2(s[:4][::-1])
        else:
            raise ValueError(f"Cannot apply idwt2 at position {pos}")
    
    subbands = s[:]
        
    layer2 = []
    for i in range(3):
        layer3 = []
        for j in range(4):
            if(i==j==0):
                subbands[3] = idwt2(subbands[:4])
            pos = 3 + 4*j + 16*i

            layer3 += [apply_idwt_at_pos(subbands[pos:pos+4], j)]

        layer2 += [apply_idwt_at_pos(layer3, i)]
            
    
    subbands[48:51] = layer2
        
    layer1 = []
    for i in range(4):
        pos = 48 + 4*i
        layer1 += [apply_idwt_at_pos(subbands[pos:pos+4], i)]

    y = idwt2(layer1)    
    return y
    