from read_rec_maps import *
from recombination import *
import os

#   Read and store the recombination maps
map_files=["../Recombination_Maps/"+file for file in os.listdir("../Recombination_Maps/") if file[0:3]=='chr'] #list the recombination map available
rec_maps=read_all_read_maps(map_files) #read the recombination maps and store it in a dictionnary
#rec_maps -> {'1':[[cM_chr1][bp_chr1]], '2':[[cM_chr2][bp_chr2]], ....}

#   Create a Chr objet for each chromosome, which will contain all interesting information
chromosomes={}
for chromosome in rec_maps:

    chromosome_size_M=max(rec_maps[chromosome][0])/100
    new_chromosome=Chr(chromosome_size_M)
    curr_chr_rec_map_M=[cM/100 for cM in rec_maps[chromosome][0]] #convet cM to M
    curr_chr_rec_map_bp=rec_maps[chromosome][1]
    new_chromosome.calculate_rec_pos_bp(curr_chr_rec_map_M, curr_chr_rec_map_bp) #convert the recombinating sites into bp
    chromosomes[chromosome]=new_chromosome

#chromosomes -> {'1': chr1, '2': chr2, ...}

#   Example of how to retrieve every recominatin sites
print([chromosomes[key].rec_pos_bp for key in chromosomes])