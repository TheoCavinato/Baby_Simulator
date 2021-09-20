import argparse

class Child():
	def __init__(self, family, id, father, mother, sex, phenotype):
		self.family, self.id, self.father, self.mother, self.sex, self.phenotype= family, id, father, mother, sex, phenotype
	def __str__(self):
		return ' '.join([self.family,self.id,self.father,self.mother,self.sex,self.phenotype])
	    

#	Read the PED file	#
parser=argparse.ArgumentParser()
parser.add_argument('-i', '--input_ped', help="Original unsorted .ped file", required=True)
parser.add_argument('-o', '--output_path', help="Path to the output sorted .ped", required=True)
args=parser.parse_args()

ped_file_stream=open(args.input_ped,'r')
output_ped=open(args.output_path,'w')

children_seen=set()
children_list=[]
for line in ped_file_stream:
	line_splt=line.split()
	new_child=Child(*line_splt)
	if new_child.mother=='NA' or new_child.father=='NA':
		output_ped.write(line)
	elif new_child.id in children_seen:
		print("WARNING: Duplicated CHILDID. Only the first one has been taken in account")
	else:
		children_seen.add(new_child.id)
		children_list.append(new_child)

#	Retrieve the parents already present in the vcf file
father_seen, mother_seen=set(),set()
for child in children_list:
	if child.father not in children_seen and child.father not in father_seen:
		father_seen.add(child.father)
	if child.mother not in children_seen and child.mother not in mother_seen:
		mother_seen.add(child.mother)
	if child.father in mother_seen or child.mother in father_seen:
		print("WARNING: of the samples present in the vcf file has both sex. His id is", child.id)

parent_seen=father_seen.union(mother_seen)

#	Sort using a queue process
while len(children_list)>0:
	curr_child=children_list.pop(0)
	if curr_child.mother in parent_seen and curr_child.father in parent_seen:
		parent_seen.add(curr_child.id)
		output_ped.write(str(curr_child))
		output_ped.write("\n")
	else:
		children_list.append(curr_child)