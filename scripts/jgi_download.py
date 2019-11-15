#!/usr/bin/env python3

import os, sys, logging, csv, re

import xml.etree.ElementTree as ET
DEBUG=1


# SOME FUNCTIONS TO HELP
# GET ALL THE FILES
def get_files(folder,dtype,species):
    for file in folder.findall('file'):
        name = file.get('label')
        url  = file.get('url')
        if name in local_names:
            name = local_names[name]

        if name not in species:
            species[name] = {dtype: [url,name]}
        elif dtype not in species[name]:
            species[name][dtype] = [url,name]
        else:
            print("warning - updating %s %s with url (%s) when it was previously %s"
                  % (name,dtype,url,species[name][dtype][0]))
            species[name][dtype] = [url, name]
    return species

#GET ALL THE TARFILES
def get_tarfiles(folder,species):
    for file in folder.findall('file'):
        name = file.get('label')
        url  = file.get('url')
        filename = file.get('filename')
        dtype = ""
        if name in local_names:
            name = local_names[name]

        m = re.search(r'(\S+).tar.gz$',filename)
        if ( m ):
            dtype = 'tar'
            # print("%s match for %s outpref = %s"%(dtype,filename,name))
        else:
            continue

        if name not in species:
            species[name] = {dtype: [url,name]}
        elif dtype not in species[name]:
            species[name][dtype] = [url,name]
        else:
            print("warning - updating %s %s with url (%s) when it was previously %s"
                  % (name,dtype,url,species[name][dtype][0]))
            species[name][dtype] = [url, name]
    return species

#GET THE ANNOTATION FILES
def get_annotation_files(folder,species):
    for file in folder.findall('file'):
        name = file.get('label')
        url  = file.get('url')
        filename = file.get('filename')
        dtype = ""
        if name in local_names:
            name = local_names[name]
        m = re.search(r'([^_]+)_(\S+).(gff|tsv|domtblout)$',filename)
        if ( m ):
            dtype = m.group(3)
            name = name + "_" + m.group(2)
            # print("%s match for %s outpref = %s"%(dtype,filename,name))
        elif re.search(r'_(prodigal|genemark)',filename):
            continue
        elif re.search(r'_contigs\.fna',filename):
            dtype = "DNA"
        elif re.search(r'_genes\.fna',filename):
            dtype = "CDS"
        elif re.search(r'_proteins\.faa',filename):
            dtype = "pep"
        else:
            #print("Cannot match filename %s to a category"%(filename))
            continue
#        print("'%s' '%s'"%(filename,name))

        if name not in species:
            species[name] = {dtype: [url,name]}
        elif dtype not in species[name]:
            species[name][dtype] = [url,name]
        else:
            print("warning - updating %s %s with url (%s) when it was previously %s"
                  % (name,dtype,url,species[name][dtype][0]))
            species[name][dtype] = [url, name]
    return species

def traverse(xmlroot,species):

    for topfolder in xmlroot.findall('folder'):
        foldername = topfolder.get('name')
        if (  foldername == 'Sequence' or foldername == 'Assembly' or
              foldername == 'Alignment'):
            for filesfolder in topfolder.findall("folder"):
                print("names are ",filesfolder.get("name"))
                if filesfolder.get('name') == "QC and Genome Assembly":
                    if foldername == 'Assembly':
                        species = get_files(filesfolder,foldername+"-dna",species)
                    for subfolder in filesfolder.findall("folder"):
                        species = get_files(subfolder,foldername,species)
                elif filesfolder.get('name') == "IMG Data":
                    if foldername == 'Assembly':
                        species = get_annotation_files(filesfolder,species)
        elif foldername == "Annotation":
            for filesfolder in topfolder.findall("folder"):
                if filesfolder.get('name') == "IMG Data":
                    species = get_annotation_files(filesfolder,species)
        elif foldername == "Other":
            for filesfolder in topfolder.findall("folder"):
                if filesfolder.get('name') == "IMG Data":
                    species = get_tarfiles(filesfolder,species)
    return species

#INIT SOME VARIABLES
hosturl='https://genome.jgi.doe.gov'

mapext = {'Sequence'   : [ 'fastq', 'corr.fastq.gz'],
          'Assembly'  : [ 'info', 'coverage.txt' ],
          'Assembly-dna'   : [ 'DNA' , 'fasta'],
          'domtblout' : ['annotation', 'domtblout'],
#          'blasttab' : ['annotation', 'blasttab'],
          'gff' : ['annotation', "gff"],
          'tsv' : ['annotation', "tsv"],
          'DNA'   : [ 'DNA' , 'contig.fasta'],
          'pep'   : [ 'pep' , 'aa.fasta'],
          'CDS'   : [ 'CDS' , 'cds.fasta'],
          'Alignment'  : [ 'SAM', 'sam.gz' ],
          'tar'   : [ 'Other', 'tar.gz' ],

      }

outdir="source/JGI"

if not os.path.exists(os.path.dirname(outdir)):
    os.mkdir(os.path.dirname(outdir))

# make the top level folder for storing the data
if not os.path.exists(outdir):
    os.mkdir(outdir)
# make the folders for storing the data - process mapext dictionary
for t in mapext.values():
    if not os.path.exists(os.path.join(outdir,t[0])):
        os.mkdir(os.path.join(outdir,t[0]))


# READ IN SOME SKIP DIRECTIVES
skipme='lib/skip_data_jgi.csv'

skip_data = {}
if os.path.exists(skipme):
    with open(skipme,'r') as skips:
        skipread = csv.reader(skips,delimiter=",")
        skipheader = next(skipread)
        for row in skipread:
            if len(row):
                skip_data[row[0]] = 1

local_names = {}

# HARCODED THE EXPECTED NAME OF THE INPUT NAMES
localnamesfile = 'lib/zygo_names.csv'

if os.path.exists(localnamesfile):
    with open(localnamesfile,'r') as newnm:
        nmread = csv.reader(newnm,delimiter="\t")
        header = next(nmread)
        for row in nmread:
            local_names[row[0]] = row[1]

names = ['mucoromycota', 'zoopagomycota']
if len(sys.argv) > 1: # if a command line argument is provided
    names = sys.argv[1:]

for base in names:
    xmlfile  = "lib/%s.xml" % (base)

    if not os.path.exists(xmlfile):
        print("expecting %s - did you run scripts/init_jgi_download.sh?"%(xmlfile))
        exit()

    species = {}
    tree = ET.parse(xmlfile)
    root = tree.getroot()
    for topfolder in root.findall('folder'):
        foldername = topfolder.get('name')
        #print("top",foldername)
        for subfolder in topfolder.findall('folder'):
            #print("\tsub1",subfolder.get('name'))
            if subfolder.get('name') == "Files":
                species = traverse(subfolder,species)

    # main code
    with open("lib/%s.csv"%(base),"w") as jgiout:
        with open("lib/%s_jgi_download.sh"%(base),"w") as dwnload:
            jgicsv = csv.writer(jgiout,delimiter=",",lineterminator="\n")
            jgicsv.writerow(['Label','read_URL','asm_URL'])
            for sp in sorted(species.keys()):
                row = ['',sp]
                for t in mapext.keys():
                    if t in species[sp]:
                        if len(row[0]) == 0:
                            row[0] = species[sp][t][1]
                        row.append(hosturl + species[sp][t][0])
                        outfile="%s.%s"%(os.path.join(outdir,mapext[t][0],
                                                      row[0]),mapext[t][1])
                        if ( not os.path.exists(outfile) and
                        not os.path.exists(outfile+".gz") ):
                            dwnload.write("curl -o %s '%s' -b cookies\n"
                                          %(outfile,row[-1]))
                    else:
                        row.append("NO_%s_URL"%(t))
