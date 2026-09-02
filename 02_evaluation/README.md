# CRAQ analysis
samtools view -@ 24 -b Y_scaffold_assembly_hifi_winnowmap_sort.bam CQ_chrY_v1.3 > CQ_chrY_v1.3_hifi_winnowmap.bam

samtools index -@ 8 CQ_chrY_v1.3_hifi_winnowmap.bam

samtools view -@ 24 -b Y_scaffold_assembly_ont_winnowmap_sort.bam CQ_chrY_v1.3 > CQ_chrY_v1.3_ont_winnowmap.bam

samtools index -@ 8 CQ_chrY_v1.3_ont_winnowmap.bam

samtools merge -@ 24 hifi_ont_winnowmap_merge.bam  CQ_chrY_v1.3_hifi_winnowmap.bam CQ_chrY_v1.3_ont_winnowmap.bam

samtools index -@ 8 hifi_ont_winnowmap_merge.bam

export PATH=~/miniconda3/envs/CRAQ/bin:$PATH

perl /data/home/wangbo/software/CRAQ-v1.10/bin/craq -g CQ_chrY_v1.3.fasta -sms hifi_ont_winnowmap_merge.bam -t 48 -D hifi_ont
# Merqury analysis
~/data/miniconda3/envs/merqury/bin/meryl k=21 count /data/DATA/ChineseQuartet/RAWDATA/ILM_PCR-free/finished/LCL7_R*.gz output illuminak21.meryl

bash ~/data/miniconda3/envs/merqury/bin/merqury.sh illuminak21.meryl Y_scaffold_assembly.fasta merqury_res_k21

