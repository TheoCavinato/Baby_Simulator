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
        return -(math.log(1-x)/lbda)
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
        #   Convert the list of recombnation sites in cM into a list of recombination sites in bp
        self.start_time=time.time()
        self.rec_pos_bp=[]
        for rec_pos in self.rec_pos_M:
            l=0
            r=len(M_list)-1
            if rec_pos < M_list[0]: # if the recombination site in M is inferior to the first recombination site in M of the recombination map
                self.rec_pos_bp.append("NA")
                break
            if rec_pos>M_list[-2]:  # if the recombination site in M is between the last element and th before last element of the M recombination map
                top_value_bp, bottom_value_bp=bP_list[-1], bP_list[-2]
                top_value_M, bottom_value_M=M_list[-1], M_list[-2]
                converted_rec_pos=((rec_pos-bottom_value_M)/(top_value_M-bottom_value_M))*(top_value_bp-bottom_value_bp)+bottom_value_bp
                self.rec_pos_bp.append(converted_rec_pos)
                break
            while True:
                mid=(l+r)//2
                top_value_M=M_list[mid]
                bottom_value_M=M_list[mid-1]
                if mid < 0:
                    self.rec_pos_bp.append("NA")
                    break
                if rec_pos > bottom_value_M and rec_pos < top_value_M:
                    top_value_bp=bP_list[mid]
                    bottom_value_bp=bP_list[mid-1]
                    converted_rec_pos=((rec_pos-bottom_value_M)/(top_value_M-bottom_value_M))*(top_value_bp-bottom_value_bp)+bottom_value_bp
                    self.rec_pos_bp.append(converted_rec_pos)
                    break
                elif rec_pos < bottom_value_M:
                    r=mid-1
                elif rec_pos > top_value_M:
                    l=mid+1