#!/usr/bin/env python
import csv, sys,re

idir='source/JGI'
odir="genomes"

namefile=sys.argv[1]
with open(namefile) as csvfile:
    rdr = csv.reader(csvfile,delimiter=",",quotechar='"')
    for row in rdr:
        pref=row[2]
        if row[1] == "name":
            continue
        name_long = re.sub(r'\"','',row[1])
        name_long = re.sub(r' ','_',name_long)
        print("pigz -dc %s/DNA/%s.nt.fasta.gz > %s/DNA/%s.nt.fasta"%(idir,pref,odir,name_long))
