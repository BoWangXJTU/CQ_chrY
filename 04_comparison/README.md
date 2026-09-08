# Collinearity comparison
minimap2 -ax asm5 -t48 --eqx chm13v2.0_chrY_reid.fa CQ_chrY_v1.3.fasta > CQ_vs_chm13.sam

minimap2 -ax asm5 -t48 --eqx chm13v2.0_chrY_reid.fa CN1.v1.0.1_chrY.fasta > CN1_vs_chm13.sam

minimap2 -ax asm5 -t48 --eqx chm13v2.0_chrY_reid.fa YAOv2.0_chrY.fasta > YAO_vs_chm13.sam

~/data/miniconda3/envs/syri/bin/syri -c CQ_vs_chm13.sam -r chm13v2.0_chrY_reid.fa -q CQ_chrY_v1.3.fasta -k -F S --prefix CQ_

~/data/miniconda3/envs/syri/bin/syri -c CN1_vs_chm13.sam -r chm13v2.0_chrY_reid.fa -q CN1.v1.0.1_chrY.fasta -k -F S --prefix CN1_

~/data/miniconda3/envs/syri/bin/syri -c YAO_vs_chm13.sam -r chm13v2.0_chrY_reid.fa -q YAOv2.0_chrY.fasta  -k -F S --prefix YAO_

~/data/miniconda3/envs/plotsr/bin/plotsr --sr CQ_syri.out --sr CN1_syri.out --sr YAO_syri.out --genomes genomes.txt -o chrY_plotsr_output.pdf -S 0.6 -W 24 -H 4 -f 8 --chr chrY

