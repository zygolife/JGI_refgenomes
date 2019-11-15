#!/usr/bin/bash
#SBATCH -p short --out logs/rename_fasta.log -n 2 --mem 2gb
perl -i -p -e 's/>jgi\|(\S+)\|(\d+)\|/>$1|$1_$2 /' source/JGI/pep/*.aa.fasta
perl -i -p -e 's/>jgi\|(\S+)\|(\d+)\|/>$1|$1_$2 /' source/JGI/CDS/*.cds.fasta
