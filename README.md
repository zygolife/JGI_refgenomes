# Downloading Zygo Genomes from JGI

```
cp scripts/init_jgi_download.sh.template init_jgi_download.sh
# edit this script to add your username and password
bash scripts/init_jgi_download.sh # will download .xml files
python scripts/jgi_download.py # to setup download steps
# now run the download steps
bash lib/mucoromycota_jgi_download.sh
bash lib/zoopagomycota_jgi_download.sh
# check these as they assume a large number of CPUS
sbatch -p short scripts/jgi_process_downloaded.sh
sbatch -p short scripts/fix_fasta_prefix.sh
```

# Downloading Zygo Genomes from  non-JGI sources

