# -*- coding: utf-8 -*-
import numpy as np
from Huffman import HuffmanCoding
from bitstring import BitStream, pack
from constants.fileformat import *

class WSQFileWriter:
    def __init__(self):
        self.buffer = BitStream()
        
    def tofile(self, filename):
        if(self.SHOW_DETAILS):
            print(f"Writing file format to: \"{filename}\"")
        self.buffer.tofile(open(filename,"wb"))
        if(self.SHOW_DETAILS):
            print(f"Done")
        return self
        
        
    def __write_signed_float(self, value, size):
        if(size<1):
            raise ValueError("Size of signed float should be greater than zero")

        sign = int(value < 0)   #1 if negative, 0 if positive
        
        self.buffer.append(pack("uint:8",sign))
        
        self.__write_unsigned_float(abs(value), size)
        
    def __write_unsigned_float(self, value, size):
        if(size<1):
            raise ValueError("Size of unsigned float should be greater than zero")
            
        exponent = 0
        mantissa = value

        if(mantissa == 0):        
            self.buffer.append(pack("int:8, uint:"+str(8*size), 0, 0))   
            return
        
        max_int = 2**(8*size) - 1
        
        while(mantissa < max_int):
            mantissa *= 10
            exponent += 1
            
        while(mantissa > max_int):
            mantissa /= 10
            exponent -= 1
        
        self.buffer.append(pack("int:8, uint:"+str(8*size), exponent, int(mantissa)))
        return
            
        
    
    def encode(self, image_path, show_details=False):
        import skimage.io as io
        from utils.dwt import getSubbandStructure, visualizeSubbandStructure
        from utils.quantization import quantization
        from SWT import dwt2
        
        self.SHOW_DETAILS = show_details
        
        img = np.float64(io.imread(image_path))

        if(self.SHOW_DETAILS):
            print(f"Image size is: {img.shape}")

        M = np.mean(img)
        R = np.max([np.max(img)-M,M-np.min(img)])/128
        
        img = (img-M)/R
        if(self.SHOW_DETAILS):
            print("Applying SWT")
        subbands = getSubbandStructure(img, dwt2)

        if(self.SHOW_DETAILS):
            print("Quantizing image, target bitrate: 0.75 bpp")
        p, Q, Z, C = quantization(subbands, 0.75)       #0.1709
        
        blocks = []
        
        blocks.append(np.concatenate([p_i.flatten() for p_i in p[:19]]).tolist())
        blocks.append(np.concatenate([p_i.flatten() for p_i in p[19:52]]).tolist())
        blocks.append(np.concatenate([p_i.flatten() for p_i in p[52:64]]).tolist())
        
        huff_encoder = HuffmanCoding()
        
        
        huffdata = []
        hufftables = []
        for i in range(3):
            a,b = huff_encoder.encode_data(blocks[i],self.SHOW_DETAILS)
            huffdata.append(a)
            hufftables.append(b)
        
        if(self.SHOW_DETAILS):
            print("Generating file format")
           
        self.buffer.append(SOI) #Start of File
        self.buffer.append(SOF) #Start of Frame
        
        self.buffer.append(pack("uint:16", 17))
        
        self.buffer.append(pack("uint:8",0))     #Black calibration - (useless)
        self.buffer.append(pack("uint:8",0))     #White calibration - (useless)
        
        self.buffer.append(pack("uint:16",img.shape[0]))     #Number of lines
        self.buffer.append(pack("uint:16",img.shape[1]))     #Number of samples per line
        
        self.__write_unsigned_float(M, 2)
        self.__write_unsigned_float(R, 2)
        
        self.buffer.append(pack("uint:8",0))        #Encoder alghorithm identifier
        self.buffer.append(pack("uint:16",0))       #Encoder software implementation identifier
        
        
        #
        #   Place all optional blocks
        #
        #DTT
        self.buffer.append(DTT)
        self.buffer.append(pack("uint:16", 58)) #size
        
        from constants.dwt import h_0, h_1 #Hard-coded transformation filters coefficients
        self.buffer.append(pack("uint:8", len(h_0)))
        self.buffer.append(pack("uint:8", len(h_1)))
        
        h_0 = h_0[4:]
        h_1 = h_1[3:]
        
        for elem in h_0:
            self.__write_signed_float(elem, 4)
            
        for elem in h_1:
            self.__write_signed_float(elem, 4)
        
        #DQT
        self.buffer.append(DQT)
        self.buffer.append(pack("uint:16", 389)) #size
        
        self.__write_unsigned_float(C, 2)
        
        for i in range(64):
            self.__write_unsigned_float(Q[i], 2)
            self.__write_unsigned_float(Z[i], 2)
        
        #DHT
        identifier=0
        for HT in hufftables:
            self.buffer.append(DHT)    
            self.buffer.append(pack("uint:16", 19 + np.sum(HT.bits))) #size
            
            self.buffer.append(pack("uint:8", identifier)) #Table Identifier
            
            for elem in HT.bits:
                self.buffer.append(pack("uint:8", elem))    
            for row in HT.huffval:
                for elem in row:
                    self.buffer.append(pack("uint:8", elem))    
                    
            identifier+=1
            
        identifier=0            
        for elem in huffdata:
            self.buffer.append(SOB)
            self.buffer.append(pack("uint:16",3))
            self.buffer.append(pack("uint:8",identifier))
            identifier+=1
            
            for byte in elem.bytes:
                if(byte == 255):
                    self.buffer.append(pack("uint:8, uint:8",255,00))
                else:
                    self.buffer.append(pack("uint:8",byte))
            
            
            
        #EOI
        self.buffer.append(EOI)
        return self
        

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d",action='store_true',dest='details',help='Show details')
    parser.add_argument("input", help="path to the input image")
    parser.add_argument("output", help="path to the output image")
    args = parser.parse_args()
    
    uncompressed_file = args.input
    WSQcompressed_file = args.output
    showDetails = args.details
    
    WSQFileWriter().encode(uncompressed_file,showDetails).tofile(WSQcompressed_file)
    
if __name__ == "__main__":
    main()
        