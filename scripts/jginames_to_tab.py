#!/usr/bin/env python3
import csv,re,sys

reader = csv.reader(sys.stdin,delimiter=",")
writer = csv.writer(sys.stdout,delimiter="\t",lineterminator='\n')
for line in reader:
    name=re.sub(" ","_",line[1])
    name=re.sub(";","",name)
    name=re.sub("\r","",name)
    writer.writerow([line[2],name])
