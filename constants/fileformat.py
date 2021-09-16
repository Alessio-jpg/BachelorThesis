# -*- coding: utf-8 -*-

SOI = b'\xFF\xA0'    #Start of Image
EOI = b'\xFF\xA1'    #End of Image
SOF = b'\xFF\xA2'    #Start of Frame
SOB = b'\xFF\xA3'    #Start of Block
DTT = b'\xFF\xA4'    #Define Transform Table
DQT = b'\xFF\xA5'    #Define Quantization Table
DHT = b'\xFF\xA6'    #Define Huffman Table(s)
DRI = b'\xFF\xA7'    #Define Restart Interval
RST0 = b'\xFF\xB0'   #Restart Modulo 0
RST1 = b'\xFF\xB1'   #Restart Modulo 1
RST2 = b'\xFF\xB2'   #Restart Modulo 2
RST3 = b'\xFF\xB3'   #Restart Modulo 3
RST4 = b'\xFF\xB4'   #Restart Modulo 4
RST5 = b'\xFF\xB5'   #Restart Modulo 5
RST6 = b'\xFF\xB6'   #Restart Modulo 6
RST7 = b'\xFF\xB7'   #Restart Modulo 7
COM = b'\xFF\xA8'    #Comment