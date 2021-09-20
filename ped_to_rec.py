from read_rec_maps import *
from chr_object import *
import argparse, os

#########################
#   USER PARAMETERS     #	
#########################

parser = argparse.ArgumentParser(description = 'To do')
parser.add_argument('-vf', '--variant_file', help = 'variant file (vcf/bcf)', required=True)
parser.add_argument('-rec_m', '--rec_path_male', help = 'male recombination map', required=True)
parser.add_argument('-rec_f', '--rec_path_female', help = 'female recombination map', required=True)
parser.add_argument('-ped', '--ped_path', help = 'ped (pedigree) file containing the parents of the simulated individuals', required=True)
parser.add_argument('-o', '--output_path', help = 'output file', required=False)
args = parser.parse_args()

#########################
#	FUNCTIONS	#
#########################

def get_inds_sex(sorted_ped_path):
	sorted_ped=open(sorted_ped_path, 'r')
	female_sex_list, male_sex_list=set(),set()

	for line in sorted_ped:
		splt_line=line.split()
		child_ID, sex= splt_line[1], int(splt_line[4])
		if sex==2:
			female_sex_list.add(child_ID)
		elif sex==1:
			male_sex_list.add(child_ID)

	return (female_sex_list, male_sex_list) 

def write_rec_file(output_path, sorted_ped_path, female_rec_map, male_rec_map):
	sorted_ped=open(sorted_ped_path, 'r')
	output_file=open(output_path,'w')

	female_list, male_list = get_inds_sex(sorted_ped_path)
	if len(female_list & male_list):
		raise Exception("Some individuals have both sex in the pedigree file")

	#ensure boh rec_map have the same chromosome
	male_chr, female_chr = set(female_rec_map.keys()), set(male_rec_map.keys())
	if male_chr != female_chr:
		common_chr = male_chr & female_chr
		female_rec_map= {chr: female_rec_map[chr] for chr in common_chr}
		male_rec_map= {chr: male_rec_map[chr] for chr in common_chr}

	#retrieve information about the rec_maps:
	all_chr_sizes = [{chr:max(female_rec_map[chr][0]) for chr in female_rec_map.keys()}, {chr:max(male_rec_map[chr][0]) for chr in male_rec_map.keys()}]
	
	for line in sorted_ped:
		new_rec_line=line.split()
		parent_1_ID, parent_2_ID = new_rec_line[2:4]
		if parent_1_ID!='NA' and parent_2_ID!='NA': # if we know the parents of the children
			for parent in [parent_1_ID, parent_2_ID]: #one recombination list per parent
				parental_rec_sites=[]

				if parent in female_list:
					rec_map=female_rec_map
					chr_sizes_M=all_chr_sizes[0]
				elif parent in male_list:
					rec_map=male_rec_map
					chr_sizes_M=all_chr_sizes[1]
				else:
					raise Exception("The sex is lacking for one of the parent, please make sure you integrated all the parents with their corresponding sex in the ped file")

				chromosomes=[int(chr.split('r')[-1]) for chr in list(rec_map.keys())] # need to sort them or they end up unsorted in the rec lists
				chromosomes.sort()
				chromosomes=["chr"+str(chr) for chr in chromosomes]
				
				for chromosome in chromosomes:
					chromosome=str(chromosome)
					chromosome_size_M=chr_sizes_M[chromosome]
					curr_chr_rec_map_M=rec_map[chromosome][0]
					curr_chr_rec_map_bp=rec_map[chromosome][1]
					new_chromosome=Chr(chromosome_size_M)
					new_chromosome.convert_cM_to_bP_binary_search(curr_chr_rec_map_M, curr_chr_rec_map_bp)
					parental_rec_sites.append(new_chromosome.rec_pos_bp)
				new_rec_line.append(str(parental_rec_sites))
			output_file.write("\t".join(new_rec_line))
			output_file.write("\n")

#########################
#	ACTUAL CODE	#
#########################

#read the corresponding recombination map
male_rec_map=read_rec_map(args.rec_path_male)
female_rec_map=read_rec_map(args.rec_path_female)
write_rec_file(args.output_path, args.ped_path, female_rec_map, male_rec_map)