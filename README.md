# CQ_chrY
T2T Genome Assembly of the Y Chromosome in the Chinese Quartet

01. Genome assemble
   
export PATH=~/miniconda3/envs/hifiasm0.25.0/bin:$PATH
hifiasm -o CQ_chrY -t114 --ont LCL7.all.pass.fastq.gz
