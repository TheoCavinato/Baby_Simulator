from pysam import VariantFile
import ast, time, argparse, random

#########################
#	USER ARGUMENTS	#
#########################
parser=argparse.ArgumentParser()
parser.add_argument('-rec', '--rec_file_path', required=True)
parser.add_argument('-vcf', '--vcf_file_path', required=True)
args=parser.parse_args()

#########################################
#	Read the recombination file	#
#########################################

rec_file=open(args.rec_file_path,'r')

child_parents, child_rec_sites, all_childs = {}, [], []
for line in rec_file:
	splt_line=line.split("\t")
	child_id, dad_id, mom_id, = splt_line[1:4]

	dad_rec, mom_rec=splt_line[-2:]
	all_childs.append(child_id)

	child_parents[child_id]=(dad_id, mom_id)
	child_rec_sites.append((ast.literal_eval(dad_rec), ast.literal_eval(mom_rec)))

#je test des trucs
child_rec_sites_reversed=[[] for i in range(len(child_rec_sites[0][0]))]
for chr in range(len(child_rec_sites[0][0])):
	for ind_rec_site in child_rec_sites:
		child_rec_sites_reversed[chr].append((ind_rec_site[0][chr], ind_rec_site[1][chr]))

#################################
#	Read the variant file	#
#################################

def haplo_chooser(meiosis,toss_result):
	haplo_chooser=toss_result
	if len(meiosis):
		for rec_site in meiosis:
			if rec_site > position:
				break
			haplo_chooser+=1
	else:
		return 0

	if haplo_chooser%2:
		return 0
	else:
		return 1


bcf_input= VariantFile(args.vcf_file_path)

output_samples=list(bcf_input.header.samples)
original_alleles_length=len(output_samples)
output_samples.extend(all_childs)

original_header=str(bcf_input.header)[:-1]
original_header+='\t'+'\t'.join(all_childs)

print(original_header)

toss_a_coin=[random.randint(0,1) for _ in range(len(all_childs))]
parent_position=[(output_samples.index(child_parents[child][0]), output_samples.index(child_parents[child][1])) for child in all_childs]

for record in bcf_input:
	start_time=time.time()
	chromosome = int(record.chrom)-1
	position = record.pos
	alleles=[al['GT'] for al in record.samples.values()]

	#-------Take the allele resulting from the meiosis for each parent-------#

	for parent_rec_sites, pos_parent, toss in zip(child_rec_sites_reversed[chromosome], parent_position, toss_a_coin):
		out_haplo=tuple(str(alleles[parent][haplo_chooser(rec_sites, toss)]) for rec_sites, parent in zip(parent_rec_sites, pos_parent))
		alleles.append(out_haplo)
			
	#-------Add alleles to the record-------#

	output_record=str(record)[:-1]
	output_record+='\t'+'\t'.join(['|'.join(all) for all in alleles[original_alleles_length:]])
	print(output_record)
