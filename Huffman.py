# -*- coding: utf-8 -*-
import numpy as np
from collections import namedtuple


HuffTable = namedtuple("HuffTable", "bits huffval")


class HuffTree(object):
    def __init__(self ):
        self.__root = None
        return
    
    def generate_from_frequencies(self, freq):
        elems = freq[:]
        while(len(elems) > 1):
            elems = sorted(elems, key=lambda x: x[1], reverse=True)
            
            l1 = elems[-1]      #Pop least frequent elems
            l2 = elems[-2]
            elems = elems[:-2]
            
            elems += [(self.HuffNode(l1[0],l2[0]), l1[1] + l2[1])]
            
        self.__root = elems[0][0]
        
        return
    
    def get_HuffTable(self):
        if(self.__root is None):
            raise ValueError("HuffTree not initialized")
        bits = [0] * 16
        huffval = [list() for _ in range(16)]
        
        def explore_hufftree(node, depth = 0):
            if(not isinstance(node, self.HuffNode)):
                bits[depth-1] += 1
                huffval[depth-1].append(node)
                return
                
            (left, right) = node.children()
            
            explore_hufftree(left, depth+1)
            explore_hufftree(right, depth+1)
            return
        
        explore_hufftree(self.__root)
        
        return HuffTable(bits, huffval)
    
    class HuffNode(object):
        def __init__(self, left, right):
            self.left = left
            self.right = right
            """
            distL = 0
            distR = 0
            
            if(isinstance(left, self.__class__)):
                distL = left.max_dist_from_leaf + 1 
            if(isinstance(right, self.__class__)):
                distR = right.max_dist_from_leaf + 1 
                
            self.max_dist_from_leaf = max(distL,distR)
            """
            return
              
        def children(self):
            return self.left,self.right
        
        def __repr__(self):
            return f"{self.left}, {self.right}"
                    
            

class HuffmanCoding:
            
    def __init__(self):
        pass


    """ Encoding """
    def encode_data(self, block, verbose = False):
        freq = dict()    
    
        zero_run = 0
        for n in block:
            if(n>74):
                if(n<256):
                    freq[101] = freq.get(101,0) + 1
                else:
                    freq[103] = freq.get(103,0) + 1
            elif(n<-73):
                if(n>-256):
                    freq[102] = freq.get(102,0) + 1
                else:
                    freq[104] = freq.get(104,0) + 1
            
            if(n == 0):
                zero_run +=1
            if(n != 0):
                if(zero_run != 0):
                    if(zero_run > 100):
                        if(zero_run<256):
                            freq[105] = freq.get(105,0) + 1
                        else:
                            freq[106] = freq.get(106,0) + 1
                    else:
                        freq[zero_run] = freq.get(zero_run,0) + 1
                
                zero_run = 0     
                
                if(n>74):
                    if(n<256):
                        freq[101] = freq.get(101,0) + 1
                    else:
                        freq[103] = freq.get(103,0) + 1
                elif(n<-73):
                    if(n>-256):
                        freq[102] = freq.get(102,0) + 1
                    else:
                        freq[104] = freq.get(104,0) + 1
                else:
                    freq[n+180] = freq.get(n+180,0) + 1
            
                
        if (zero_run != 0):
            if(zero_run > 100):
                if(zero_run<256):
                    freq[105] = freq.get(105,0) + 1
                else:
                    freq[106] = freq.get(106,0) + 1
            else:
                freq[zero_run] = freq.get(zero_run,0) + 1
            
           
        freq = sorted(freq.items(), key=lambda x: x[1], reverse=True)
        
        tree = HuffTree()
        
        tree.generate_from_frequencies(freq)
        
        bits, huffval = tree.get_HuffTable()
        
        huffsize = self.generate_size_table(bits, huffval)
        huffcode = self.generate_code_table(huffsize)
        ehufco, ehufsi = self.order_codes(huffval, huffcode, huffsize)
        
        huffdict = dict()
        for i in range(len(ehufco)):
            if(ehufco[i] is not None):
                huffdict[i] = format(ehufco[i],"0"+str(ehufsi[i])+"b")
        
        from bitstring import BitStream, pack
        coded_data = BitStream()
        
        zero_run = 0
        for n in block:
            if(n == 0):
                zero_run += 1
            else:
                if(zero_run != 0):
                    if(zero_run < 101):
                        coded_data.append("0b" + huffdict[zero_run])
                    elif(zero_run< 256):
                        coded_data.append("0b" + huffdict[105])
                        coded_data.append(pack("uint:8", zero_run))
                    else:
                        coded_data.append("0b" + huffdict[106])
                        coded_data.append(pack("uint:16", zero_run))
                
                zero_run=0
                if(n>74):
                    if(n<256):
                        coded_data.append("0b" + huffdict[101])
                        coded_data.append(pack("uint:8",n))
                    else:
                        coded_data.append("0b" + huffdict[103])
                        coded_data.append(pack("uint:16",n))
                elif(n<-73):
                    if(n>-256):
                        coded_data.append("0b" + huffdict[102])
                        coded_data.append(pack("uint:8",-1 * n))
                    else:
                        coded_data.append("0b" + huffdict[104])
                        coded_data.append(pack("uint:16",-1 * n))
                else:
                    coded_data.append("0b" + huffdict[n+180])

        #Check if a zero_run needs to be coded at the end
        if(zero_run != 0):
            if(zero_run < 101):
                coded_data.append("0b" + huffdict[zero_run])
            elif(zero_run< 256):
                coded_data.append("0b" + huffdict[105])
                coded_data.append(pack("uint:8", zero_run))
            else:
                while(zero_run>0):
                    coded_data.append("0b" + huffdict[106])
                    coded_data.append(pack("uint:16", min(zero_run, 2**16-1)))
                    zero_run -= 2**16-1
        while(coded_data.len%8 != 0):
            coded_data.append("0b1")
        
        if(verbose):
            print(f"Encoded {coded_data.len} bits")
        #raise Exception()
        return coded_data, HuffTable(bits, huffval)
    
    

            
            
    """
    def get_bits_from_hufftree(self, node, bitstring=''):        
        if(bitstring == ''):
            self.bits = [0]*16
            self.huffval = [list() for _ in range(16)]
        
        if(not isinstance(node, HuffNode)):
            self.bits[len(bitstring)-1] += 1
            self.huffval[len(bitstring)-1].append(node)
            return {node: bitstring}

        (left, right) = node.children()
        
        l_bit = "1"
        r_bit = "0"
        
        hufftree = dict()

        if(isinstance(right, HuffNode) and right.max_dist_from_leaf == node.max_dist_from_leaf -1):
            l_bit, r_bit = r_bit, l_bit
            hufftree.update(self.get_bits_from_hufftree(left,bitstring + l_bit))
            hufftree.update(self.get_bits_from_hufftree(right,bitstring + r_bit))
        else:
            hufftree.update(self.get_bits_from_hufftree(right,bitstring + r_bit))
            hufftree.update(self.get_bits_from_hufftree(left,bitstring + l_bit))
        
        return hufftree
    """ 
        
    """ Decoding"""
    def decode_text(self, text, bits, huffval, verbose = True):
        huffsize = self.generate_size_table(bits, huffval)
        huffcode = self.generate_code_table(huffsize)
        ehufco, ehufsi = self.order_codes(huffval, huffcode, huffsize)
        
        reverse_dict = dict()
        for i in range(len(ehufco)):
            if(ehufco[i] is not None):
                reverse_dict[format(ehufco[i],"0"+str(ehufsi[i])+"b")] = i
        
        
        code = ""
        decompressed_block = []
        

        from bitstring import BitStream
        data = BitStream(text)
        #data.replace("0xff00","0xff")
        
        data.pos = 0
        
        if (verbose):
            print("Decoding ",len(data)," bits . . . ")
        
        while(data.pos != data.len):
            if(data.pos == data.len-32):
                _p = data.pos
                data.pos = _p
            bit = data.read(1).bin
            code += bit
            if(code in reverse_dict):
                decoded_code = reverse_dict[code]
                code = ""
                if(decoded_code<1 or decoded_code == 180):
                    raise ValueError()
                elif(decoded_code < 101):
                    elem = [0] * decoded_code
                elif (decoded_code <= 104):
                    sign = -1 if (decoded_code)%2==0 else 1
                    if(decoded_code <= 102):
                        elem = [sign * data.read("uint:8")]
                    else:
                        elem = [sign * data.read("uint:16")]
                elif(decoded_code == 105):
                    elem = [0] * data.read("uint:8")
                elif(decoded_code == 106):
                    elem = [0] * data.read("uint:16")
                else:
                    elem = [decoded_code-180]
                    
                decompressed_block += elem
                
        assert(code == len(code) * "1")
        return decompressed_block
        
    def generate_size_table(self, bits, huffval):

        huffsize = []
        for i in range(len(bits)):
            huffsize += [i+1] * bits[i]
    
        return huffsize

    def generate_code_table(self, huffsize):
        code = k = 0
        prev = huffsize[0]
        
        size = len(huffsize)
        huffcode = [0]*size
        while(k<size):
            if(huffsize[k] == prev):
                huffcode[k] = code
                code+=1
                k+=1
            else:
                code = code << 1
                prev += 1
                
        return huffcode
    
    def order_codes(self, huffval, huffcode, huffsize):
        k=0
        huffval_f = np.concatenate([np.asarray(huffval[i]).flatten() for i in range(len(huffval))])
        
        ehufco = [None]*256        #CODE
        ehufsi = [None]*256        #CODE_SIZE (bits)
        
        while(k<len(huffsize)):
            i = int(huffval_f[k])
            ehufco[i] = huffcode[k]
            ehufsi[i] = huffsize[k]
            
            k+=1

        return ehufco, ehufsi
        