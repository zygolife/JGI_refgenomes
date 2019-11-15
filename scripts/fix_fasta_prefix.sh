#!/usr/bin/bash
#SBATCH -p short --out logs/rename_fasta.log

pushd source/JGI/pep
for file in *.aa.fasta
do
	perl -i.bak -p -e 's/^>([^_\|]+)_/>$1|$1_/' $file 
done
popd

pushd source/JGI/CDS
for file in *.cds.fasta
do
	perl -i.bak -p -e 's/^>([^_\|]+)_/>$1|$1_/' $file
done
popd
