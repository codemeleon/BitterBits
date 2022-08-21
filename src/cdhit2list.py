lst_dict = {}
with open(cdhit_out) as fin:
    tlst = []
    tvar = ""
    i = False
    for line in fin:
        if line[0] == ">":
            if i:
                lst_dict[tvar] = tlst
                tlst = []
                tvar = ""
            i = True
            continue
        data = line.split()
        sequenceid = data[2][1:-3]
        tlst.append(sequenceid)
        if line[-2] == "*":
            tvar = sequenceid

lst_dict[tvar] = tlst
