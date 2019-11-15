#!/usr/bin/env python3

import os, sys, logging, csv, re

import xml.etree.ElementTree as ET
DEBUG=1

long2shortasm = { "Mitochondrial assembly": 'mito',
                  "Assembled scaffolds (unmasked)": 'nuclear' }

def get_assemblies (fileroot):
    for dir in fileroot.findall("folder"):
#        print("->",dir.tag,dir.attrib)
        dirname = dir.get('name')
        if ( dirname == "Mitochondrial assembly" or
             dirname == "Assembled scaffolds (unmasked)"):
            asmtype = long2shortasm[dirname]

            for file in dir:
                # print(file.tag, file.attrib)
                url = file.get('url')
                filename = file.get('filename')
                label = file.get('label')
                name = ascii(label)
#                if filename.startswith('Aciri1_iso'):
#                    name=name + " iso"
                thisasmtype = asmtype
                if( "MitoAssembly" in filename or
                    'mito.scaffolds' in filename or
                    'MitoScaffolds' in filename):
                    thisasmtype = 'mito'
                elif "PlasmidAssembly" in filename:
                    thisasmtype = 'plasmid'

                if name not in species:
                    species[name] = {thisasmtype: [url,filename]}
#                    print("storing %s for %s file=%s"%(name,thisasmtype,filename))
                elif thisasmtype not in species[name]:
                    species[name][thisasmtype] = [url,filename]
#                    print("storing %s for %s file=%s"%(name,thisasmtype,filename))
                else:
                    print("skipping this data as already have a type %s stored for %s - trying to save %s" %(thisasmtype,name,filename))
                    print(species[name])

# we probably should defined a priority here based on age of file
# primary alleles over secondary, etc
# for now relying on the XML file order and how it is parsed as to who
# wins
def get_annotations(fileroot):
    for dir in fileroot.findall("folder"):
        print("->",dir.tag,dir.attrib)
        dirname = dir.get('name')

        for file in dir:
            if file.tag == "file":
                url = file.get('url')
                filename = file.get('filename')
                if ("secondary_alleles" in filename or
                    'Secondary_Alleles' in filename):
                    continue

                label = file.get('label')
                name = ascii(label)
#                if filename.startswith('Aciri1_iso'):
#                    name=name + " iso"

                dtype = ''

                if filename.endswith("gff3.gz") or filename.endswith("gff.gz"):
                    dtype = 'gff'
                elif dirname == "CDS":
                    dtype = 'CDS'

                if dtype:
                    if name not in species:
                        species[name] = {dtype: [url,filename]}
#                        print("storing %s for %s file=%s"%(name,dtype,filename))
                    elif dtype not in species[name]:
                        species[name][dtype] = [url,filename]
#                        print("storing %s for %s file=%s"%(name,dtype,filename))
                    else:
                        # replace gff3 with gff when available
                        if (filename.endswith(".gff3.gz") and
                            species[name][dtype][1].endswith(".gff.gz")):
                            species[name][dtype] = [url,filename]
                        else:
                            print("skipping this data as already have a type %s stored for %s - trying to save %s" %(dtype,name,filename))
                            print(species[name])



skipme='lib/skip_data_jgi.csv'

skip_data = {}
if os.path.exists(skipme):
    with open(skipme,'r') as skips:
        skipread = csv.reader(skips,delimiter=",")
        skipheader = next(skipread)
        for row in skipread:
            if len(row):
                skip_data[row[0]] = 1

outdir="source/JGI"
if not os.path.exists(os.path.dirname(outdir)):
    os.mkdir(os.path.dirname(outdir))
if not os.path.exists(outdir):
    os.mkdir(outdir)
for t in ['CDS','GFF','DNA','MT','pep']:
    if not os.path.exists(os.path.join(outdir,t)):
        os.mkdir(os.path.join(outdir,t))


hosturl='https://genome.jgi.doe.gov'

# HARCODED THE EXPECTED NAME OF THE INPUT NAMES
localnamesfile = 'lib/zygo_names.csv'
#local_names = {}
#if os.path.exists(localnamesfile):
#    with open(localnamesfile,'r') as newnm:
#        nmread = csv.reader(newnm,delimiter="\t")
#        header = next(nmread)
#        for row in nmread:
#            local_names[row[0]] = row[1]

species = {}

names = ['mucoromycota', 'zoopagomycota']
if len(sys.argv) > 1: # if a command line argument is provided
    names = sys.argv[1:]

for base in names:
    xmlfile  = "lib/%s.xml" % (base)

    tree = ET.parse(xmlfile)
    root = tree.getroot()

    # main code
    for topfolder in root.findall('folder'):
        if topfolder.get('name') == 'Unknown':
            for filesfolder in topfolder.findall("folder"):
                if filesfolder.get('name') == "Files":
                    for files in filesfolder.findall("folder"):
                        print("\t",files.tag,files.attrib)
                        if files.get('name') == "Assembly":
                            get_assemblies(files)
                        elif files.get('name') == "Annotation":
                            for annot in files.findall('folder'):
                                if annot.get('name') == 'Filtered Models ("best")':
                                    get_annotations(annot)

mapext = {'nuclear': [ 'DNA', 'nt.fasta.gz'],
          'mito'   : [ 'MT' , 'mt.fasta.gz'],
          'gff'    : [ 'GFF', 'gff3.gz' ],
          'CDS'    : [ 'CDS', 'cds.fasta.gz'] }

for base in names:
    xmlfile  = "lib/%s.xml" % (base)
    # main code
    with open("lib/%s.csv"%(base),"w") as jgiout:
        with open("lib/%s_jgi_download.sh"%(base),"w") as dwnload:
            jgicsv = csv.writer(jgiout,delimiter=",",lineterminator="\n")
            jgicsv.writerow(['Prefix','Species','Full Name','DNA_URL','MITO_URL',
                             'GFF_URL','CDS_URL'])
            for sp in sorted(species.keys()):
                s_sp = re.sub(r'\s*v\d+\.\d+\s*','',sp)
                s_sp = re.sub(r'\s*\(Environmental single-cell\)\s*','',s_sp)
                s_sp = re.sub(r'\s+',' ',s_sp)
                s_sp = re.sub(r'\'','',s_sp)
                print(s_sp)
                m = re.split(r'\s+',s_sp,maxsplit=3)
                print(m)
                row = [''," ".join(m[0:2]),sp]
                for t in ['nuclear','mito','gff','CDS']:
                    if t in species[sp]:
                        if len(row[0]) == 0:
                            row[0] = species[sp][t][0].split('/')[2]
                        row.append(hosturl + species[sp][t][0])
                        outfile="%s.%s"%(os.path.join(outdir,mapext[t][0],
                                                      row[0]),mapext[t][1])
                        if not os.path.exists(outfile):
                            dwnload.write("curl -o %s '%s' -b cookies\n"
                                          %(outfile,row[-1]))
                    else:
                        row.append("NO_%s_URL"%(t))
                if row[0] not in skip_data:
                    jgicsv.writerow(row)
