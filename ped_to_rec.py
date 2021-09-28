from read_rec_maps import *
from chr_object import *
import argparse, os, time

#########################
#   USER PARAMETERS     #	
#########################

parser = argparse.ArgumentParser(description = 'To do')
parser.add_argument('-rec', '--rec_path', help = 'recombination folder', required=True)
parser.add_argument('-ped', '--ped_path', help = 'ped (pedigree) file containing the parents of the simulated individuals', required=True)
parser.add_argument('-o', '--output_path', help = 'output file', required=False)
args = parser.parse_args()

#########################
#	FUNCTIONS	#
#########################

def write_rec_file(output_path, sorted_ped_path, rec_maps):
	sorted_ped=open(sorted_ped_path, 'r')
	output_file=open(output_path,'w')

	#retrieve chromosome sizes about the rec_maps:
	all_chr_sizes_M = {chr:max(rec_maps[chr][0]) for chr in rec_maps.keys()}
	
	for line in sorted_ped:
		new_rec_line=line.split()
		parent_1_ID, parent_2_ID = new_rec_line[2:4]
		if parent_1_ID!='NA' and parent_2_ID!='NA': # if we know the parents of the children
			for sex_parent in range(2): #one recombination list per parent
				parental_rec_sites=[]

				# sort chromosomes by assending order
				chromosomes = [int(chr) for chr in rec_maps.keys()]
				chromosomes.sort()
				
				for chromosome in chromosomes:
					
					chromosome=str(chromosome)
					chromosome_size_M=all_chr_sizes_M[chromosome]
					curr_chr_rec_map_M=rec_maps[chromosome][0]
					curr_chr_rec_map_bp=rec_maps[chromosome][1]
					new_chromosome=Chr(chromosome_size_M)
					new_chromosome.convert_cM_to_bp(curr_chr_rec_map_bp, curr_chr_rec_map_M)
					parental_rec_sites.append(new_chromosome.rec_pos_bp)

				new_rec_line.append(str(parental_rec_sites))

			output_file.write("\t".join(new_rec_line))
			output_file.write("\n")

#########################
#	ACTUAL CODE	#
#########################

rec_maps_folder = args.rec_path
map_files=[rec_maps_folder+"/"+file for file in os.listdir(rec_maps_folder)] #list the recombination map available
recombination_maps=read_all_read_maps(map_files) #read the recombination maps and store it in a dictionnary
recombination_maps.pop('X')
recombination_maps[str(1)]

write_rec_file(args.output_path, args.ped_path, recombination_maps)