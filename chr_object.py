import random as rand
import math
import time
from read_rec_maps import *
import time
import os

#   This modified Chr object will have for different functions, one for each process
class Chr:
    def __init__(self, chr_size_M ):
        self.chr_size_M=chr_size_M
        self.rec_pos_M=[]
        self.generate_rec_pos()
        self.rec_pos_bp=[]

    def inverse_CDF(self,x, lbda):
        return -(math.log(1-x)/1)
    def generate_rec_pos(self):
        #   Create a list containing the recombintation positions in cM
        rec_pos=0
        iterator=0
        while True:
            if rec_pos==0:
                rec_pos=self.inverse_CDF(rand.uniform(0,1), self.chr_size_M)
            else:
                rec_pos=self.rec_pos_M[(iterator-1)]+self.inverse_CDF(rand.uniform(0,1), self.chr_size_M)
            if rec_pos > self.chr_size_M:
                break
            self.rec_pos_M.append(rec_pos)
            iterator+=1

    def convert_cM_to_bP_binary_search(self, M_list, bP_list):
            #   Convert the list of recombnation sites in M into a list of recombination sites in bp
            self.rec_pos_bp=[]
            for rec_pos in self.rec_pos_M:

                # if element not comprised between the maximum and minimum value of the Recombination Map's M_list
                if rec_pos < M_list[0] or rec_pos > M_list[-1]:
                    self.rec_pos_bp.append("None")
                    break

                low = 0  
                high = len(M_list) - 1  
                mid = 0  

                while low <= high:  
                    # for get integer result   
                    mid = (high + low) // 2  


                    # Check if n is present at mid   
                    if M_list[mid] < rec_pos:  
                        low = mid + 1  

                    # If n is greater, compare to the right of mid   
                    elif M_list[mid] > rec_pos:  
                        high = mid - 1  

                    # If n is smaller, compared to the left of mid  
                    else:  
                        self.rec_pos_bp.append(bP_list[mid])
                        break

                # element was not present in the list, return the elements between which it is comprised
                top_value_bp, bottom_value_bp = bP_list[low], bP_list[high]
                top_value_M, bottom_value_M = M_list[low], M_list[high]
                converted_rec_pos=((rec_pos-bottom_value_M)/(top_value_M-bottom_value_M))*(top_value_bp-bottom_value_bp)+bottom_value_bp
                self.rec_pos_bp.append(converted_rec_pos)
