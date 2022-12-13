!/bin/bash

wget -O datasets https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets
wget -O dataformat https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/dataformat

chmod +x dataformat
chmod +x datasets

./datasets download virus genome taxon sars-cov-2 --complete-only
unzip -p ncbi_dataset.zip ncbi_dataset/data/genomic.fna | gzip -c > all_genomes.fa.gz