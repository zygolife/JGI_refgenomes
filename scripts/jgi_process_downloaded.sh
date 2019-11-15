#!/usr/bin/bash
#SBATCH --nodes 1 --ntasks 48 -p stajichlab --mem 32g --out logs/jgi_process.log
module load parallel
module unload perl
translate=$(which bp_translate_seq.pl)
if [ -z $translate ]; then
	translate=$(realpath scripts/bp_translate_seq.pl)
fi
pushd source/JGI/CDS
echo $translate
parallel -j 48 'gzip -dc {} > {.}' ::: *.cds.fasta.gz
parallel -j 48 "cat {} | $translate --out  ../pep/{= s:\.cds\.fasta:.aa.fasta: =}" ::: *.cds.fasta
popd
