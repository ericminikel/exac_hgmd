#!/bin/bash

# set up some variables
gatkjar=/humgen/gsa-hpprojects/GATK/bin/current/GenomeAnalysisTK.jar
b37ref=/humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta


# HGMD Annotations
# ----------------
# The HGMD SNP mutations are from version 2012.4 (i.e. HGMD IDs that are prefixed with 
# CM and CS). The HGMD insertion and deletion mutions are from version 2012.1 (i.e. 
# HGMD IDs that are prefixed with CI and CD). Lastly, the HGMD gene annotations are 
# from version 2012.1.

# Table 4. HGMD annotations.
# HGMD Annotation Description
# HGMD_MUT Overlaps HGMD mutation and annotated with HGMD ID
# HGMD_SITE Overlaps HGMD mutation site but allele does not match
# HGMD_GENE Overlaps Gene with HGMD mutation

# find one example of this annotation
zcat $annovcf | grep -v ^# | head -5647 | grepk HGMD
# example: HGMD_SITE=CM128668;HGMD_MUT=CM128668
# note that the above example is from 1 949739  .   G   T,A
# i.e. there are two alt alleles but the HGMD fields are not sub-divided per allele
# therefore it is ambiguous which allele matches the HGMD mutation

# and if you try this:
zcat $annovcf | grep -v ^# | head -5650 | grepk HGMD
# there are two entries that just say HGMD_GENE

mkdir -p jobfiles

# convert whole VCF to table
# first, try to figure out full list of fields in INFO
zcat $annovcf | grep -v ^# | head -10000 | cut -f8 | perl -pe 's/=.*?(;|\n)/\n/g' | sed 's/;/\n/g' | sed 's/|/\n/g' | sort | uniq > annovcf_info_fields.txt
# then construct a GATK VariantsToTable command with all of them
echo -e "java -Xmx8g -jar $gatkjar \
    -R $b37ref \
    -T VariantsToTable \
    -V $annovcf \
    --variant_index_type LINEAR \
    --allowMissingData \
    --splitMultiAllelic \
    --showFiltered \
    -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER \\" > variantstotable.bash
cat annovcf_info_fields.txt | awk '{print "-F "$1" \\"}' >> variantstotable.bash
echo -e "-o annovcf.table" >> variantstotable.bash
bsub -P $RANDOM -q bhour -W 04:00 -o jobfiles/variantstotable.out -e jobfiles/variantstotable.err bash variantstotable.bash
# the above job took ~20 mins

wc -l annovcf.table
# 11993923 annovcf.table

# HGMD fields are cols 30-32 in resulting table
HGMD_GENE_colno=`cat annovcf.table | head -1 | tr '\t' '\n' | grep -n HGMD_GENE | sed 's/:.*//'`
HGMD_MUT_colno=`cat annovcf.table | head -1 | tr '\t' '\n' | grep -n HGMD_MUT | sed 's/:.*//'`
HGMD_SITE_colno=`cat annovcf.table | head -1 | tr '\t' '\n' | grep -n HGMD_SITE | sed 's/:.*//'`

# create a file of just the HGMD-relevant sites
# only keep rows where >0 of these cols are non-blank
cat annovcf.table | awk -F"\t" -v col1=$HGMD_GENE_colno \
    -v col2=$HGMD_MUT_colno \
    -v col3=$HGMD_SITE_colno \
    '$col1 != "NA" || $col2 != "NA" || $col3 != "NA" {print $0}' > annovcf.hgmdonly.table
# the above command takes ~5 mins to run

# create a separate list of multi-allelic sites
echo -e "CHROM\tPOS" > multiallelicsites.txt
cat annovcf.table | cut -f1,2 | uniq -d >> multiallelicsites.txt
# that also takes ~5 mins

# number of lines per file as of 2014-07-04:
wc -l annovcf.table
# 11993923 annovcf.table
wc -l annovcf.hgmdonly.table
# 2730324 annovcf.hgmdonly.table
wc -l multiallelicsites.txt
# 1004223 multiallelicsites.txt

# analysis plan:
# - add to the HGMD file a column indicating whether multi-allelic
# - sort on AC
# - check if multi-allelics are enriched
# - check if kinases are enriched

# questions for group
# - do we have the HGMD database?
# - do we know whether the annotator correctly handles non-minimal representations?

# try to get tables small enough to manage with R
# - the dbSNP id is per site not variant, so that's useless. 
cat annovcf.hgmdonly.table | head -1 | tr '\t' '\n' | \
    egrep -n -w "CHROM|POS|REF|ALT|QUAL|FILTER|AC|AF|HGMD_GENE|HGMD_MUT|HGMD_SITE|HGNC_GENE|LOF|LOF_FLAG" | \
    sed 's/:.*//' > desired_column_indices.txt
awk_col_list=`cat desired_column_indices.txt | sed 's/^/$/' | tr '\n' ',' | sed 's/,$//'`
echo "awk '{print "$awk_col_list"}'" > awk_command.bash
head annovcf.hgmdonly.table | bash awk_command.bash

/humgen/atgu1/fs03/lek/HGMD/Original/2012.4/hgmd-2012.4-alleles.zip








