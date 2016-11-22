LOBREF=/home/mag50/hs37d5_v3.0.3/lobstr_v3.0.2_hg19_ref_nochr.bed
TMPLOC=/n/scratch2/gymrek
REFDB=/groups/reich/melissa/dbase/GRch37/hs37d5.fa
PROPFILE=/groups/reich/melissa/mutation_rate_project/autosomal_byclass/explore_classes/lobSTR_ref_GRCh37_properties.bed
TORFILE=/groups/reich/melissa/mutation_rate_project/analysis_round1/autosomal_estimates/constraint2/features/replication_timing/Koren_etal_TableS2_hg19.bed
INTERGENIC=/groups/reich/melissa/mutation_rate_project/autosomal_byclass/explore_classes/lobSTR_ref_GRCh37_intergenic.bed

BASEDIR=/groups/reich/melissa/mutation_rate_project/analysis_round1/
PREDFILE=${BASEDIR}/constraint/autosomal_perlocus_estimates.bed
#DATADIR=/oasis/projects/nsf/ddp268/mgymrek/mutea-autosomal/sgdp_asdt_vcf

NUMSIM=10 # For making null, how many samples to draw
BATCHSIZE=100 # Number of loci to include in each batch for making the null. Actual Mutea batch size is $NUMSIM*$BATCHSIZE
