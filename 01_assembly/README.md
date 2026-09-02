# ONT-UL only assembly
export PATH=~/miniconda3/envs/hifiasm0.25.0/bin:$PATH
hifiasm -o L7 -t114 --ont /data/DATA/ChineseQuartet/RAWDATA/ONTUL_100kb/LCL7/LCL7.all.pass.fastq.gz
# HiFi ONT-UL HiC assembly
hifiasm -o L7 -t96 --dual-scaf --telo-m CCCTAA --h1 /data/DATA/ChineseQuartet/RAWDATA/HiC/PLCL7/PLCL7-1-1/R1.fq.gz,/data/DATA/ChineseQuartet/RAWDATA/HiC/PLCL7/PLCL7-1-2/R1.fq.gz --h2 /data/DATA/ChineseQuartet/RAWDATA/HiC/PLCL7/PLCL7-1-1/R2.fq.gz,/data/DATA/ChineseQuartet/RAWDATA/HiC/PLCL7/PLCL7-1-2/R_2.fq.gz --ul /data/DATA/ChineseQuartet/RAWDATA/ONTUL_100kb/LCL7/PLCL7.100kb.pass.fastq.gz /data/DATA/ChineseQuartet/RAWDATA/HiFi/fastqs/LCL7/*.fq.gz
# Gap filling
python fill_gap.py # Use utg000737l.fa to fill gaps in CQ_chrY_v1.2_1gap.fasta.
