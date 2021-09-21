# Baby_Simulator #
Based on recombination maps, a pedigree and a vcf/bcf, <b>Baby_Simulator</b> simulate genomes.

# How to #
Here we imagine you have a .ped "example_ped.ped" and a bcf file "example_chr22.bcf" as the one present in Example/ \
`$ python3 sort_ped.py -i Example/example_ped.ped -o Example/sorted_example_ped.ped` \
`$ python3 ped_to_rec.py -vf Example/example_chr22.bcf -rec_m Recombination_Maps/male.gmap -rec_f Recombination_Maps/female.gmap -ped Example/sorted_example_ped.ped -o example_rec.rec`