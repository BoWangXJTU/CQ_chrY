# ONT-UL only assembly
export PATH=~/miniconda3/envs/hifiasm0.25.0/bin:$PATH
hifiasm -o L7_ont -t114 --ont /data/DATA/ChineseQuartet/RAWDATA/ONTUL_100kb/LCL7/LCL7.all.pass.fastq.gz
# HiFi ONT-UL HiC assembly
hifiasm -o L7_hifi_ont_hic -t96 --dual-scaf --telo-m CCCTAA --h1 /data/DATA/ChineseQuartet/RAWDATA/HiC/PLCL7/PLCL7-1-1/R1.fq.gz,/data/DATA/ChineseQuartet/RAWDATA/HiC/PLCL7/PLCL7-1-2/R1.fq.gz --h2 /data/DATA/ChineseQuartet/RAWDATA/HiC/PLCL7/PLCL7-1-1/R2.fq.gz,/data/DATA/ChineseQuartet/RAWDATA/HiC/PLCL7/PLCL7-1-2/R_2.fq.gz --ul /data/DATA/ChineseQuartet/RAWDATA/ONTUL_100kb/LCL7/PLCL7.100kb.pass.fastq.gz /data/DATA/ChineseQuartet/RAWDATA/HiFi/fastqs/LCL7/*.fq.gz
# Gap filling
python fill_gap.py # Use utg000737l.fa to fill gaps in CQ_chrY_v1.2_1gap.fasta.
# Construct the non-Y scaffold genome sequence
nohup ragtag.py scaffold -t 48 -o hap1_rag -r ~/chm13-v2.0/chm13v2.0.fa L7_hifi_ont_hic.haplotype1.fasta > hap1.log 2>&1 &
nohup ragtag.py scaffold -t 48 -o hap2_rag -r ~/chm13-v2.0/chm13v2.0.fa L7_hifi_ont_hic.haplotype2.fasta > hap2.log 2>&1 &
cat hap1_ragtag.scaffold.fasta hap2_ragtag.scaffold.fasta | awk 'BEGIN{RS=">"; FS="\n"} !/chrY_RagTag/{printf ">%s", $0}' > non-Y_scaffold_assembly.fasta
# Build the complete genome assembly including both non-Y scaffolds and the curated Y chromosome
cat non-Y_scaffold_assembly.fasta CQ_chrY_v1.3.fasta > Y_scaffold_assembly.fasta # Y_scaffold_assembly.fasta will be used for subsequent read alignment and Merqury analysis.

