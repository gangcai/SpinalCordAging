outfile = open("DEG_Summary.tsv","w")
i = 0
summary = {}
for line in open("SCAge_DEGs_WilCox_Sig.tsv"):
    i += 1
    if i == 1:
        continue
    items = line.rstrip().split("\t")
    fc = items[4]
    cname = items[2]


    if not cname in summary:
        summary[cname] = [0,0]
    fc = float(fc)
    if fc > 0:
        summary[cname][0] += 1
    else:
        summary[cname][1] += 1

outfile.write("cellType\tDEGType\tDEGNum\n")
for cname in summary.keys():
    up_num, down_num = summary[cname]
    outfile.write(cname + "\tUP\t"+str(up_num)+"\n")
    outfile.write(cname + "\tDOWN\t" + str(down_num) +  "\n")

outfile.close()
