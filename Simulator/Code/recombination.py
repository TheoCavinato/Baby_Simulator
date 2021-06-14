import random as rand
import math

#   Define the recombination sites: là c pour un chr seulement
class Chr:
    def __init__(self, chr_size_M ):
        self.chr_size_M=chr_size_M
        self.rec_pos_M=[]
        self.generate_rec_pos()
        self.rec_pos_bp=None

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

    def  find_M_range_in_Map(self, x, M_list):
        if x > max(M_list) or x < min(M_list):
            return "number out of range"
        
        if x in M_list:
            return x

        while True:
            iterator=round(len(M_list)/2)
            bottom_index=(iterator-1)
            top_index=(iterator)
            if bottom_index < 0:
                bottom_index=0

            if x >= M_list[bottom_index] and x <= M_list[top_index]:
                break
            elif x < M_list[bottom_index]:
                M_list=M_list[0:(bottom_index+1)]
            elif x > M_list[top_index]:
                M_list=M_list[(top_index-1):]
        return [M_list[bottom_index], M_list[top_index]]

    def convert_cM_to_bP(self, x, M_range, M_list, bP_list):
        if isinstance(M_range, list):
            max_M=max(M_range)
            min_M=min(M_range)
            max_bP=bP_list[M_list.index(max_M)]
            min_bP=bP_list[M_list.index(min_M)]
            fraction=(x-min_M)/(max_M-min_M)
            return (min_bP+fraction*(max_bP-min_bP))
        else:
            return bP_list[M_list.index(x)]

    #   Transform the cM rcombinating positions into bp and stock the result in rec_pos_bP!
    def calculate_rec_pos_bp(self, rec_map_cM, rec_map_bp):
        M_ranges = [self.find_M_range_in_Map(rec_pos,rec_map_cM) for rec_pos in self.rec_pos_M]
        self.rec_pos_bp=[self.convert_cM_to_bP(rec_pos, M_ranges[n], rec_map_cM, rec_map_bp) for n, rec_pos in enumerate(self.rec_pos_M)]

