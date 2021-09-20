def read_rec_map(rec_map_path):
    rec_maps={}

    #   Open the file and read the header
    rec_map_file=open(rec_map_path,'r')
    header=rec_map_file.readline().strip().split('\t')
    col_cM, col_pB, col_chr=header.index("cM"), header.index("pos"), header.index("chr")

    #   Read the file line by line
    while True:
        line=rec_map_file.readline()
        line_splt=line.strip().split("\t")

        #   Stop reading if we reach the end of the file
        if not line:
            break

        #   Retrieve values
        value_cM, value_pB, value_chr = line_splt[col_cM], line_splt[col_pB], str(line_splt[col_chr])

        try:
            rec_maps.setdefault(value_chr, [[],[]])[0].append(float(value_cM)/100+last_cM)
            rec_maps[value_chr][1].append(float(value_pB))
            last_cM = rec_maps[value_chr][0][-1]
    
        except ValueError: 
            print("Warning: " , line_splt[col_cM] ," at position ", line_splt[col_pB]) 
            last_cM = 0
        

    return rec_maps