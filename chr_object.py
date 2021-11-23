import random
import sys
import math
import time
from read_rec_maps import *
import time
import os

#   This modified Chr object will have for different functions, one for each process
class Chr:
    def __init__(self, chr_size_M, rng):
        self.chr_size_M=chr_size_M
        self.rng = rng
        self.rec_pos_M=[]
        self.generate_rec_pos()
        self.rec_pos_bp=[]

    def inverse_CDF(self,x):
        return -(math.log(1-x)/1)

    def generate_rec_pos(self):
        #   Create a list containing the recombintation positions in cM
        rec_pos=0
        iterator=0
        while True:
            if rec_pos==0:
                rec_pos=self.inverse_CDF(self.rng.uniform(0,1))
            else:
                rec_pos=self.rec_pos_M[(iterator-1)]+self.inverse_CDF(self.rng.uniform(0,1))
            if rec_pos > self.chr_size_M:
                break
            self.rec_pos_M.append(rec_pos)
            iterator+=1

    def convert_cM_to_bp(self, list_bp, list_M):
        rec_pos, pos = 0, 0
        if len(self.rec_pos_M): # if recombinating sites created, use them
            while pos < len(list_M)-1:
                while self.rec_pos_M[rec_pos] > list_M[pos] and self.rec_pos_M[rec_pos] < list_M[pos+1]:
                    conversion=self.linear_conversion(self.rec_pos_M[rec_pos], list_M[pos], list_M[pos+1], list_bp[pos], list_bp[pos+1])
                    self.rec_pos_bp.append(conversion)
                    rec_pos+=1
                    #	If we reach the end of the recombinating positons generated
                    if rec_pos == len(self.rec_pos_M):
                        break
                
                #	If we reach the end of the recombinating positons generated
                if rec_pos == len(self.rec_pos_M):
                    break
                pos+=1

    def linear_conversion(self, X, A_1, A_2, B_1, B_2):
        return ((X-A_1)/(A_2-A_1))*(B_2-B_1)+B_1
