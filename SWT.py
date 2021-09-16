# -*- coding: utf-8 -*-
import numpy as np

def dwt(x):
    from constants.dwt import h_0,h_1
    
    
    N = len(x)
    y = np.concatenate((x,x[1:-1][::-1]))    

    Y = np.fft.fft(y)
    H0 = np.fft.fft(h_0,2*N-2)
    H1 = np.fft.fft(h_1,2*N-2)

    t = np.asarray([np.exp(4 * 1j * 2 * np.pi * k/(2*N-2)) for k in range(2*N-2)])
    
    L = np.real(np.fft.ifft(Y * H0 * t))[::2]
    H = np.real(np.fft.ifft(Y * H1 * t))[::2]

    L = L [:N//2]
    H = H [:N//2]
    

    return L,H
    
def dwt2(x):
    from constants.dwt import h_0,h_1
    M,N = x.shape
    
    #
    #  Columns
    #

    y = np.concatenate((x,x[1:-1][::-1]),axis=0)
    
    Y = np.fft.fft(y,axis=0)
    H0 = np.fft.fft(h_0,2*M-2)
    H1 = np.fft.fft(h_1,2*M-2)
    
    t = np.asarray([np.exp(4 * 1j * 2 * np.pi * k/(2*M-2)) for k in range(2*M-2)])
    
    l = np.real(np.fft.ifft(np.einsum("ij,i->ij", Y, H0 * t),axis=0))[::2]
    h = np.real(np.fft.ifft(np.einsum("ij,i->ij", Y, H1 * t),axis=0))[::2]
    
    if(M%2==0):
        l = l[:M//2]
        h = h[:M//2]
    else:
        l = l[:(M+1)//2]
        h = h[:(M-1)//2]
    
    #
    #  Rows
    #

    
    l = np.concatenate((l,l[:,1:-1][:,::-1]),axis=1)
    h = np.concatenate((h,h[:,1:-1][:,::-1]),axis=1)

    
    L = np.fft.fft(l,axis=1)
    H = np.fft.fft(h,axis=1)
    
    H0 = np.fft.fft(h_0, 2*N-2)
    H1 = np.fft.fft(h_1, 2*N-2)
    
    t = np.asarray([np.exp(4 * 1j * 2 * np.pi * k/(2*N-2)) for k in range(2*N-2)])

    LL = np.einsum("ij,j->ij", L, H0 * t)
    LH = np.einsum("ij,j->ij", L, H1 * t)
    HL = np.einsum("ij,j->ij", H, H0 * t)
    HH = np.einsum("ij,j->ij", H, H1 * t)
    
    ll = np.real(np.fft.ifft(LL,axis=1))[:,::2]
    hl = np.real(np.fft.ifft(HL,axis=1))[:,::2]
    lh = np.real(np.fft.ifft(LH,axis=1))[:,::2]
    hh = np.real(np.fft.ifft(HH,axis=1))[:,::2]
    
    if(N%2==0):
        ll = ll[:,:N//2]
        hl = hl[:,:N//2]
        lh = lh[:,:N//2]
        hh = hh[:,:N//2]
    else:
        ll = ll[:,:(N+1)//2]
        hl = hl[:,:(N+1)//2]
        lh = lh[:,:(N-1)//2]
        hh = hh[:,:(N-1)//2]

    return [ll,lh,hl,hh]
    

    
    
def idwt(x):
    def get_analisys_filters(h_0,h_1):
        asgn = [{0: -1, 1: 1}[k % 2] for k in range(len(h_1))]
        f_0 = np.concatenate(([0],asgn * h_1))
        
        asgn = [{0: 1, 1: -1}[k % 2] for k in range(len(h_0))]
        f_1 = np.concatenate(([0],asgn * h_0))
        
        return f_0, f_1
        
    from constants.dwt import h_0,h_1
    from utils.dwt import expand
    
    l,h = x
    
    N = len(l)
    
    l = np.concatenate((l,l[1:][::-1]))

    h = np.concatenate((h,h[:-1][::-1]))
    
    l = expand(l,2)
    h = expand(h,2)
    
    L = np.fft.fft(l)
    H = np.fft.fft(h)

    f_0,f_1 = get_analisys_filters(h_0, h_1)

    F0 = np.fft.fft(f_0,4*N-2)
    F1 = np.fft.fft(f_1,4*N-2)
    
    t = np.asarray([np.exp(4 * 1j * 2 * np.pi * k/(4*N-2)) for k in range(4*N-2)])
    
    out = np.real(np.fft.ifft(t*(L * F0 + H * F1)))[:2*N]
    
    return out
   
def idwt2(x, synthesis_lowpass = None, syntesis_hipass = None):         #Wavelet WSS
    from utils.dwt import expand2D
    
    def get_analisis_filters(h_0,h_1):
        asgn = [{0: -1, 1: 1}[k % 2] for k in range(len(h_1))]
        f_0 = np.concatenate(([0],asgn * h_1))
        
        asgn = [{0: 1, 1: -1}[k % 2] for k in range(len(h_0))]
        f_1 = np.concatenate(([0],asgn * h_0))
        
        return f_0,f_1
    
    if ((isinstance(x,tuple) or isinstance(x,list)) and len(x) == 4):
        ll,lh,hl,hh = x
    else:
        raise ValueError(f"Unable to apply idwt2 at type: {type(x)}")
    
    M_low,N_low = ll.shape
    M_high,N_high = hh.shape
    
    if(synthesis_lowpass is None or syntesis_hipass is None):
        from constants.dwt import h_0,h_1
        synthesis_lowpass = h_0
        syntesis_hipass = h_1
    
    f_0, f_1 = get_analisis_filters(synthesis_lowpass, syntesis_hipass)
    
    #
    #   Rows
    #
    
    if(N_low == N_high):
        N = 2*N_low-1
        ll = np.concatenate((ll,ll[:,1:][:,::-1]),axis=1)
        hl = np.concatenate((hl,hl[:,1:][:,::-1]),axis=1)
        lh = np.concatenate((lh,lh[:,:-1][:,::-1]),axis=1)
        hh = np.concatenate((hh,hh[:,:-1][:,::-1]),axis=1)
    else:
        N = 2*N_low-2
        ll = np.concatenate((ll,ll[:,1:-1][:,::-1]),axis=1)
        hl = np.concatenate((hl,hl[:,1:-1][:,::-1]),axis=1)
        lh = np.concatenate((lh,lh[:,::-1]),axis=1)
        hh = np.concatenate((hh,hh[:,::-1]),axis=1)
        
    ll = expand2D(ll, 2, axis=1)
    hl = expand2D(hl, 2, axis=1)
    lh = expand2D(lh, 2, axis=1)
    hh = expand2D(hh, 2, axis=1)
    
    LL = np.fft.fft(ll,axis=1)
    HL = np.fft.fft(hl,axis=1)
    LH = np.fft.fft(lh,axis=1)
    HH = np.fft.fft(hh,axis=1)
    
    F0 = np.fft.fft(f_0, 2*N)
    F1 = np.fft.fft(f_1, 2*N)
    
    t = np.asarray([np.exp(4 * 1j * 2 * np.pi * k/(2*N)) for k in range(2*N)])    
    
    L = np.einsum("ij,j->ij", LL, F0 * t) + np.einsum("ij,j->ij", LH, F1 * t)
    H = np.einsum("ij,j->ij", HL, F0 * t) + np.einsum("ij,j->ij", HH, F1 * t)

    l = np.real(np.fft.ifft(L,axis=1))[:,:(N+1)]
    h = np.real(np.fft.ifft(H,axis=1))[:,:(N+1)]
    
    #
    #   Columns
    #
    
    if(M_low == M_high):
        M = 2*M_low-1
        l = np.concatenate((l,l[1:][::-1]),axis=0)
        h = np.concatenate((h,h[:-1][::-1]),axis=0)
    else:
        M = 2*M_low-2
        l = np.concatenate((l,l[1:-1][::-1]),axis=0)
        h = np.concatenate((h,h[::-1]),axis=0)
    l = expand2D(l, 2, axis=0)
    h = expand2D(h, 2, axis=0)
    
    L = np.fft.fft(l,axis=0)
    H = np.fft.fft(h,axis=0)
    
    
    F0 = np.fft.fft(f_0, 2*M)
    F1 = np.fft.fft(f_1, 2*M)
    
    t = np.asarray([np.exp(4 * 1j * 2 * np.pi * k/(2*M)) for k in range(2*M)])
    
    Y = np.einsum("ij,i->ij", L, F0 * t) + np.einsum("ij,i->ij", H, F1 * t)
    
    x = np.real(np.fft.ifft(Y,axis=0))[:(M+1)]
    return x

