# Regional annotation
We use the annotation pipeline developed by the T2T consortium: [https://github.com/arangrhie/T2T-chrY/blob/main/src/annotate.sh]

#Centromere annotation (CenMAP)
snakemake -c 48 -p --workflow-profile none --configfile config-chrYv1.3.yaml --show-failed-logs --conda-cleanup-pkgs cache

bash hmmer-run.sh CQchrYv1.3/ AS-HORs-hmmer3.3.2-120124.hmm 48 (HOR annotation)

# Gene annotation
liftoff -g ref_annotation_coding.gff -sc 0.8 -copies -cds -o CQ_v1.3_assembly_chrY_liftoff_0.8_cds.gff -u unmapped_features_0.8_cds.txt -p 48 -f feature_types.txt CQ_chrY_v1.3.fasta chm13v2.0_chrY.fa
