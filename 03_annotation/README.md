# Regional annotation
We used the Y chromosome region annotation pipeline developed by the T2T-Y consortium (https://github.com/arangrhie/T2T-chrY/blob/main/src/annotate.sh
) to annotate the sequence classes of CQ_chrY v1.3. After annotation, we performed the following validation steps:

1. Self dot-plot alignment​ of CQ_chrY v1.3 to verify internal consistency of the annotated regions.

2. Dot-plot alignment against CHM13v2Y​ to compare the annotated regions between the two assemblies.

3. Back-mapping of our annotated intervals to CHM13v2Y​ using minimap2, followed by manual inspection and boundary refinement to ensure annotation accuracy.

This multi-step validation process ensured that the regional annotation boundaries were correctly placed and consistent with the reference.


#  Centromere annotation (CenMAP)
snakemake -c 48 -p --workflow-profile none --configfile config-chrYv1.3.yaml --show-failed-logs --conda-cleanup-pkgs cache

bash hmmer-run.sh CQchrYv1.3/ AS-HORs-hmmer3.3.2-120124.hmm 48 (HOR annotation)

# Gene annotation
liftoff -g ref_annotation_coding.gff -sc 0.8 -copies -cds -o CQ_v1.3_assembly_chrY_liftoff_0.8_cds.gff -u unmapped_features_0.8_cds.txt -p 48 -f feature_types.txt CQ_chrY_v1.3.fasta chm13v2.0_chrY.fa
