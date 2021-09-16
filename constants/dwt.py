# -*- coding: utf-8 -*-

import numpy as _np

A = _np.cbrt((-14 * _np.sqrt(15)+ 63 ) / (1080* _np.sqrt(15)))
B = _np.cbrt((-14 * _np.sqrt(15)- 63 ) / (1080* _np.sqrt(15)))

x_1 = A + B - 1/6
x_2 = complex(-(A+B)/2-1/6,_np.sqrt(3)*(A-B)/2)


h_0 = _np.array([-5 * _np.sqrt(2) * x_1/64,                                                 #-4
                -5 * _np.sqrt(2) * x_1 * _np.real(x_2)/8,                                   #-3
                -5 * _np.sqrt(2) * x_1 * ( 4*abs(x_2)**2 + 4*_np.real(x_2)-1)/16,           #-2
                -5 * _np.sqrt(2) * x_1 * ( 8*abs(x_2)**2 - _np.real(x_2))/8,                #-1
                -5 * _np.sqrt(2) * x_1 * (48*abs(x_2)**2 - 16*_np.real(x_2)+3)/32,          #0
                -5 * _np.sqrt(2) * x_1 * ( 8*abs(x_2)**2 - _np.real(x_2))/8,                #1
                -5 * _np.sqrt(2) * x_1 * ( 4*abs(x_2)**2 + 4*_np.real(x_2)-1)/16,           #2
                -5 * _np.sqrt(2) * x_1 * _np.real(x_2)/8,                                   #3
                -5 * _np.sqrt(2) * x_1/64])                                                 #4


h_1 = _np.array([-_np.sqrt(2)/(64*x_1),                     #-4
                _np.sqrt(2)* (2*x_1 + 1)/(32*x_1),          #-3
                -_np.sqrt(2)*(16*x_1 - 1)/(64*x_1),         #-2
                _np.sqrt(2)* (6*x_1 - 1)/(16*x_1),          #-1
                -_np.sqrt(2)*(16*x_1 - 1)/(64*x_1),         #0
                _np.sqrt(2)* (2*x_1 + 1)/(32*x_1),          #1
                -_np.sqrt(2)/(64*x_1)])                     #2


                      #-3,-2,-1, 0, 1, 2, 3, 4, 5
f_0 = h_1 * _np.array([-1, 1,-1, 1,-1, 1,-1,])          # h_1(n-1) * (-1)^n
f_1 = h_0 * _np.array([ 1,-1, 1,-1, 1,-1, 1,-1, 1,])    # h_0(n-1) * (-1)^(n-1)
"""
#
#   HSS filter bank
#


h_1 =  [ 0.03226944131446922,   #-4
        -0.05261415011924844,   #-3
        -0.18870142780632693,   #-2
         0.60328894481393847,   #-1
        -0.60328894481393847,   # 0
         0.18870142780632693,   # 1
         0.05261415011924844,   # 2
        -0.03226944131446922    # 3
        ]

h_0 = [   0.07565691101399093,  #-4
         -0.12335584105275092,  #-3
         -0.09789296778409587,  #-2
          0.85269867900940344,  #-1
          0.85269867900940344,  # 0
         -0.09789296778409587,  # 1
         -0.12335584105275092,  # 2
          0.07565691101399093   # 3
          ]

                      #-4,-3,-2,-1, 0, 1, 2, 3,
f_0 = h_1 * _np.array([ 1,-1, 1,-1, 1,-1, 1,-1,])    # h_1(n-1) * (-1)^n
f_1 = h_0 * _np.array([-1, 1,-1, 1,-1, 1,-1, 1,])    # h_0(n-1) * (-1)^(n-1)

"""

filter_bank_path = [
    [0b00,0b00,0b00,0b00,0b00,],
    [0b00,0b00,0b00,0b00,0b10,],
    [0b00,0b00,0b00,0b00,0b01,],
    [0b00,0b00,0b00,0b00,0b11,],
    [0b00,0b00,0b00,0b10],
    [0b00,0b00,0b00,0b01],
    [0b00,0b00,0b00,0b11],
    [0b00,0b00,0b10,0b10],
    [0b00,0b00,0b10,0b00],
    [0b00,0b00,0b10,0b11],
    [0b00,0b00,0b10,0b01],
    [0b00,0b00,0b01,0b01],
    [0b00,0b00,0b01,0b11],
    [0b00,0b00,0b01,0b00],
    [0b00,0b00,0b01,0b10],
    [0b00,0b00,0b11,0b11],
    [0b00,0b00,0b11,0b01],
    [0b00,0b00,0b11,0b10],
    [0b00,0b00,0b11,0b00],
    [0b00,0b10,0b10,0b00],
    [0b00,0b10,0b10,0b10],
    [0b00,0b10,0b10,0b01],
    [0b00,0b10,0b10,0b11],
    [0b00,0b10,0b00,0b10],
    [0b00,0b10,0b00,0b00],
    [0b00,0b10,0b00,0b11],
    [0b00,0b10,0b00,0b01],
    [0b00,0b10,0b11,0b01],
    [0b00,0b10,0b11,0b11],
    [0b00,0b10,0b11,0b00],
    [0b00,0b10,0b11,0b10],
    [0b00,0b10,0b01,0b11],
    [0b00,0b10,0b01,0b00],
    [0b00,0b10,0b01,0b10],
    [0b00,0b10,0b01,0b00],
    [0b00,0b01,0b01,0b00],
    [0b00,0b01,0b01,0b10],
    [0b00,0b01,0b01,0b01],
    [0b00,0b01,0b01,0b11],
    [0b00,0b01,0b11,0b10],
    [0b00,0b01,0b11,0b00],
    [0b00,0b01,0b11,0b11],
    [0b00,0b01,0b11,0b01],
    [0b00,0b01,0b00,0b01],
    [0b00,0b01,0b00,0b11],
    [0b00,0b01,0b00,0b00],
    [0b00,0b01,0b00,0b10],
    [0b00,0b01,0b10,0b11],
    [0b00,0b01,0b10,0b01],
    [0b00,0b01,0b10,0b10],
    [0b00,0b01,0b10,0b00],
    [0b00,0b11],
    [0b10,0b10],
    [0b10,0b00],
    [0b10,0b11],
    [0b10,0b01],
    [0b01,0b01],
    [0b01,0b11],
    [0b01,0b00],
    [0b01,0b10],
    [0b11,0b11],
    [0b11,0b01],
    [0b11,0b10],
    [0b11,0b00],
    ]


__all__ = ["h_0", "h_1", "f_0", "f_1", "filter_bank_path"]