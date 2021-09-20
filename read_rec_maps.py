def read_rec_map(map_files):
    rec_maps={}

    for rec_map_path in map_files:
        rec_map_file=open(rec_map_path,'r')
        header=rec_map_file.readline().strip().split('\t')
        col_cM, col_pB, col_chr=header.index("cM"), header.index("pos"), header.index("chr")

        while True:
            line=rec_map_file.readline()
            line_splt=line.strip().split("\t")
            if not line:
                break
            
            value_cM, value_pB, value_chr = line_splt[col_cM], line_splt[col_pB], str(line_splt[col_chr])

            try:
                rec_maps.setdefault(value_chr, [[],[]])[0].append(float(value_cM)/100)
                rec_maps[value_chr][1].append(float(value_pB))
        
            except ValueError: 
                print("Warning: " , line_splt[col_cM] ," at position ", line_splt[col_pB]) 
    return rec_maps
        