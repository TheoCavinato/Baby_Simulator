from read_rec_maps import *
from chr_object import *
import argparse, os

"""
python3 laboratory.py --rp ../Recombination_Maps -sp ../Intermediary_files/sorted_ind.ped > test2.rec
"""
#########################
#   USER PARAMETERS     #	
#########################

parser = argparse.ArgumentParser(description = 'To do')
parser.add_argument('-rp', '--rec_path', help = 'folder with the recombination maps', required=True)
parser.add_argument('-sp', '--sorted_ped_path', help = 'ped (pedigree) file containing the parents of the simulated individuals', required=True)
parser.add_argument('-m', '--m_suffix', help = 'suffix of male recombination maps', required=False)
parser.add_argument('-f', '--f_suffix', help = 'suffix of female recombination maps', required=False)
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

def write_rec_file(output_path, sorted_ped_path, *rec_maps):
	sorted_ped=open(sorted_ped_path, 'r')
	output_file=open(output_path,'w')

	#verify if the sex as an importance
	if len(rec_maps)>1:
		print("You are using different rec file for each sex")
		female_list, male_list = get_inds_sex(sorted_ped_path)
		#make sure both sex have the same keys in their recombination maps
		female_rec_map, male_rec_map={chr:rec_maps[0][chr] for chr in rec_maps[0]} ,{chr:rec_maps[1][chr] for chr in rec_maps[1]} 
		female_rec_map.pop("chrX")

	#retrieve information about the rec_maps:
	all_chr_sizes=[]
	for rec_map in rec_maps:
		all_chr_sizes.append({chr:max(rec_map[chr][0]) for chr in rec_map.keys()})
	
	for line in sorted_ped:
		new_rec_line=line.split()
		parent_1_ID, parent_2_ID = new_rec_line[2:4]
		if parent_1_ID!='NA' and parent_2_ID!='NA':
			for parent in [parent_1_ID, parent_2_ID]: #one recombination list per parent
				parental_rec_sites=[]
				if len(rec_maps)>1:
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
				else:
					rec_map=rec_maps[0]
					chr_sizes_M=all_chr_sizes[0]
					chromosomes=[int(chr.split('r')[-1]) for chr in list(rec_map.keys())] # need to sort them or they end up unsorted in the rec lists
					chromosomes.sort()

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

if args.m_suffix and args.f_suffix:
	#list the recombination maps available
	male_map_files=[args.rec_path+file for file in os.listdir(args.rec_path) if file[0:len(args.m_suffix)]==args.m_suffix]
	female_map_files=[args.rec_path+file for file in os.listdir(args.rec_path) if file[0:len(args.f_suffix)]==args.f_suffix]
	#read the corresponding recombination map
	male_rec_maps=read_rec_map(male_map_files)
	female_rec_maps=read_rec_map(female_map_files)
	write_rec_file(args.output_path, args.sorted_ped_path, female_rec_maps, male_rec_maps)
elif args.m_suffix or args.f_suffix:
	raise Exception("Lacking suffix for one of the sex: please use -m and -f arguments together, but not one alone") 
else:
	#list the recombination maps available
	map_files=[args.rec_path+file for file in os.listdir(args.rec_path)] 
	#read the corresponding recombination map
	rec_maps=read_rec_map(map_files) 
	write_rec_file(args.output_path, args.sorted_ped_path, rec_maps)