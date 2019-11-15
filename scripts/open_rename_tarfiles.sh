#!/usr/bin/bash
#SBATCH -p short -N 1 -n 4 --mem 8gb --out logs/unzip_tarfiles.log

module load parallel

pushd source/JGI/Other
for tar in *.tar.gz
do
	label=$(basename $tar .tar.gz)
	names=$(pigz -dc $tar | tar xvf -)
	for name in $names
	do
		folder=$(dirname $name)
		break
	done
	#label=$(grep taxon_name $folder/*.config | perl -p -e '$_ = (split(/\s+/,$_))[-1]."\n"')
	echo "$label,$folder"
	rsync -a $folder/$folder.a.faa  ../pep/$label.aa.fasta
	rsync -a $folder/$folder.a.fna  ../CDS/$label.cds.fasta
	for file in $(ls $folder | grep -v "a\.f[an]a" )
	do
		b=$(echo $file | perl -p -e 's/^(\d+)\.//')
		echo "$folder/$file ../annotation/$label.$b"
		rsync -a $folder/$file ../annotation/$label.$b
	done
done
