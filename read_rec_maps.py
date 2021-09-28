def read_one_rec_map(file_path):
    rec_map_file=open(file_path,'r')

    header=rec_map_file.readline().strip().split('\t')
    col_cM=header.index("cM")
    col_pB=header.index("pos")
    col_chr=header.index("chr")

    pos_M=[]
    pos_pB=[]

    rec_map_info={}

    first_line=True
    while True:
        line=rec_map_file.readline()
        line_splt=line.strip().split("\t")
        if not line:
            rec_map_info[str(chr)]=[pos_M,pos_pB]
            break

        pos_M.append(float(line_splt[col_cM])/100)
        pos_pB.append(float(line_splt[col_pB]))
        chr=line_splt[col_chr]
    
    return rec_map_info

def read_all_read_maps(map_files):
    rec_maps={}
    concatenated_maps=[]
    for rec_map in map_files:
        rec_maps.update(read_one_rec_map(rec_map))
    return(rec_maps)