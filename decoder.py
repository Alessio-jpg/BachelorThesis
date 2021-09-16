# -*- coding: utf-8 -*-
import numpy as np
from Huffman import HuffmanCoding, HuffTable
from collections import namedtuple
from constants.fileformat import *

class WSQFileReader:
    def __init__(self,):
        self.h_0 = None     #Lowpass Filter
        self.h_1 = None     #Highpass Filter
        self.TT = None      #Transformation Table
        self.Q = None       #Quantization Table
        self.HuffTables = dict()        #Huffman Tables
        
        
    
    def decode(self, filename, show_details=False):
        self.SHOW_DETAILS = show_details
        
        
        self.__file = open(filename,"rb")
        self.__file.seek(0)
        
        #
        #   Start of Image
        #
        if(self.__file.read(2) != SOI):
            raise ValueError("Image is not in the proper format")
        
        self.__checkOptBlocks()
        
        #
        #   Start of Frame
        #
        if(self.__file.read(2) != SOF):
            raise ValueError("Image is not in the proper format")
        
        size = self.__read_uint(2)-2
        
        A = self.__read_uint(1)     #Black calibration - (useless)
        B = self.__read_uint(1)     #White calibration - (useless)
        size-=2
        
        Y = self.__read_uint(2)     #Number of lines
        X = self.__read_uint(2)     #Number of samples per line
        size-=4
        
        self.Y = Y
        self.X = X
        
        M = self.__read_unsigned_float(2)
        R = self.__read_unsigned_float(2)
        size-=6
        
        Ev = self.__read_uint(1)    #Encoder alghorithm identifier
        Sf = self.__read_uint(2)    #Encoder software implementation identifier
        size-=3
        
        assert(size==0)
        
        self.M = M
        self.R = R
        
        self.__checkOptBlocks()
        
        #
        #   Start of Block
        #
        huff_decoder = HuffmanCoding()
        blocks = []
        
        
        header = self.__file.read(2)
        while(header == SOB):
            size = self.__read_uint(2)-2
            Td = self.__read_uint(1)
            
            HT = self.HuffTables.get(Td, None)
            
            if(HT is None):
                raise ValueError(f"Huffman table #{Td} is not present")
                
            size-=1
            
            assert(size==0)
            

            hufftext = self.__read_huffman_bytes()
            

            blocks.append(huff_decoder.decode_text(hufftext, HT.bits, HT.huffval, self.SHOW_DETAILS))

            self.__checkOptBlocks()
            header = self.__file.read(2)
        
        #
        #   End of Image
        #
        
        if(header != EOI):
            raise ValueError("Image is not in the proper format")
        
        self.blocks = blocks
        self.X = X
        self.Y = Y
        
        #
        # Actual decoding
        #
        from utils.dwt import getImageFromSubbandStructure
        subbands = [0]*64    
        prev = 0
        from constants.dwt import filter_bank_path
        for i in range(64):
            M,N = self.Y,self.X
            for elem in filter_bank_path[i]:
                lowM = -1 if elem & 0b1  else 1
                lowN = -1 if elem & 0b10 else 1
                
                M = (M + (M%2!=0) * lowM)/2
                N = (N + (N%2!=0) * lowN)/2
                
            M,N = int(M), int(N)
            
            if(self.Q[i] != 0):
                size = M*N
                subbands[i] = np.asarray((self.blocks[0] + self.blocks[1] + self.blocks[2])[prev: prev + size]).reshape((M,N))
                prev +=size
            else:
                subbands[i] = np.zeros((M,N))

    
        from utils.quantization import iquantization
        
        subbands = iquantization(subbands, self.Q, self.Z, self.C)
        
        self.img = np.uint8((getImageFromSubbandStructure(subbands) * self.R) + self.M)
        
        return self
    
    def tofile(self, filename):
        
        if(self.__file.tell() == 0):
            raise Exception("File not decoded!")
        
        if(self.SHOW_DETAILS):
            print(f"Saving image \"{filename}\" . . .")
        import skimage.io as io
        io.imsave(filename,self.img)
        if(self.SHOW_DETAILS):
            print("Done")
        return
        
    def __read_int(self, size):
        if(size<1):
            raise ValueError("Size of uint should be greater than zero")
        return int.from_bytes(self.__file.read(size), byteorder="big", signed=True)
    
    def __read_uint(self, size):
        if(size<1):
            raise ValueError("Size of uint should be greater than zero")
        return int.from_bytes(self.__file.read(size), byteorder="big")
    
    def __read_signed_float(self, size):
        if(size<1):
            raise ValueError("Size of float should be greater than zero")
        sign = (1 if self.__read_uint(1)==0 else -1)    #1 Byte for sign
    
        return sign * self.__read_unsigned_float(size)
    
    def __read_unsigned_float(self, size):        #a.k.a. positive float
        if(size<1):
            raise ValueError("Size of unsigned float should be greater than zero")
        scale = self.__read_int(1)                     #1 Byte for exponent
        coeff = self.__read_uint(size)                 #size Bytes for mantissa
        
        return coeff * 10 ** (-scale)
    
    def __decodeComment(self):
        size = self.__read_uint(2) - 2
        if(self.SHOW_DETAILS):
            print(self.__file.read(size).decode("utf-8"))
        else:
            self.__file.seek(size,1)
        return
    
    def __decodeTransformationTable(self):
        size = self.__read_uint(2) - 2
        L0 = self.__read_uint(1)
        L1 = self.__read_uint(1)
        size-=2
        h_0 = np.zeros(L0)
        h_1 = np.empty(L1)
        
        if(L0%2==0):
            coeffs = L0//2
        else:
            coeffs = (L0+1)//2
        
        for i in range(coeffs):
            h_0[coeffs+i-1] = self.__read_signed_float(4)
            size-=6
        h_0[:coeffs-1] = h_0[coeffs:][::-1]                         #Replicate the filter by symmetry
        
        if(L1%2==0):
            coeffs = L1//2
        else:
            coeffs = (L1+1)//2
        
        for i in range(coeffs):
            h_1[coeffs+i-1] = self.__read_signed_float(4)
            size-=6
        h_1[:coeffs-1] = h_1[coeffs:][::-1]                         #Replicate the filter by symmetry    
        
        assert(size==0)
        
        self.h_0 = h_0
        self.h_1 = h_1
        return 
    
    def __decodeQuantizationTable(self):
        size = self.__read_uint(2) - 2
        
        C = self.__read_unsigned_float(2)
        size-=3
        
        Q = np.empty(64)
        Z = np.empty(64)
        for i in range(64):
            Q[i] = self.__read_unsigned_float(2)
            Z[i] = self.__read_unsigned_float(2)

            size-=6
        
        assert(size==0)
        
        self.Q = Q
        self.Z = Z
        self.C = C
        return 
    
    def __decodeHuffmanTable(self):
        size = self.__read_uint(2) - 2
        while(size != 0):
            Th = self.__read_uint(1)   #Table identifier
            self.Th = Th
            size-=1
            L = []
            for i in range(16):
                L.append( self.__read_uint(1))
                size-=1
            #
            #   Symbol-lenght assignment parameters
            #
            V = [0]*16
            for i in range(16):
                V[i] = np.array([self.__read_uint(1) for _ in range(L[i])],dtype=np.uint8)
                size-=L[i]
            
            self.HuffTables[Th] = HuffTable(L, V)
        
        assert(size==0)
        return
    
    def __checkOptBlocks(self):
        while(1):
            header = self.__file.read(2)
            if(header == COM):
                if(self.SHOW_DETAILS):
                    print("Decoding COM block")
                self.__decodeComment()
            elif(header == DTT):
                if(self.SHOW_DETAILS):
                    print("Decoding DTT block")
                self.__decodeTransformationTable()
            elif(header == DQT):
                if(self.SHOW_DETAILS):
                    print("Decoding DQT block")
                self.__decodeQuantizationTable()
            elif(header == DHT):
                if(self.SHOW_DETAILS):
                    print("Decoding DHT block")
                self.__decodeHuffmanTable()
            else:
                self.__file.seek(-2,1)
                return
    
    def __read_huffman_bytes(self):
        hufftext = b''
        while(1):
            byte = self.__file.read(1)
            if(byte == b'\xFF'):
                nextbyte = self.__file.read(1)
                if(nextbyte != b'\x00'):
                    self.__file.seek(-2,1)
                    return hufftext
            hufftext+=byte
        

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d",action='store_true',dest='details',help='Show details')
    parser.add_argument("input", help="path to the input image")
    parser.add_argument("output", help="path to the output image")
    args = parser.parse_args()
    
    WSQcompressed_file = args.input
    WSQuncompressed_file = args.output
    showDetails = args.details
    
    WSQFileReader().decode(WSQcompressed_file,showDetails).tofile(WSQuncompressed_file)

    
if __name__ == "__main__":
    main()
        