# Baby_Simulator #
Based on recombination maps, a pedigree and a vcf/bcf, <b>Baby_Simulator</b> simulate genomes.

# How to #
Here we imagine you have a .ped "example_ped.ped" and a bcf file "example_chr22.bcf" as the one present in Example/ \
`$ python3 sort_ped.py -i Example/example_ped.ped -o Example/sorted_example_ped.ped` \
`$ python3 ped_to_rec.py -rec Recombination_Maps/genetic_maps.b37/ -ped Example/sorted_example_ped.ped -o Example/example_rec.rec `
You will obtain a ".rec" file (here named "example_rec.rec") containing a list of recombinating sites in bp for each chromosome for each parent of each individual. \
Finally, to create a variant file (vcf/bcf) containing the parent individuals + the simulated children, do: \
`$ python3 write_vcf.py -rec example_rec.rec -vcf Example/example_ch22.bcf`
