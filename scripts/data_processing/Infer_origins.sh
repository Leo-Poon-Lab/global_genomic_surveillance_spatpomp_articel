# gisaid process
## input: metadata.tsv cross_check_table_ihme_input_completed.xlsx 
## output: df_meta_clean.csv BA1.csv 
Rscript  r1.gisaid_process.r 
## select clean sequences
awk -F "," '{print $1}' data/gisaid_data/df_meta_clean.csv > results/model_data/origin_infer/clean_ID1
seqkit grep -n -f clean_ID1 data/gisaid_data/sequences_rename.fasta > results/model_data/origin_infer/clean_ID1.fasta
seqkit grep -n -f  BA1.csv clean_ID1.fasta |seqkit rmdup > results/model_data/origin_infer/BA1_rmdup.fasta 
nextclade run --input-dataset  data/sars-cov-2  --output-all=Delta/Delta_rmdup.fasta  --max-band-area  5000000000  --output-selection=fasta,csv,insertions,errors
mkdir results/model_data/origin_infer/BA1/data
Rscript r2.gisaid_process.r
awk -F "\t" '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' results/model_data/origin_infer/BA1/vn_clean_gisaid.tsv|uniq > results/model_data/origin_infer/BA1/metadata.tsv

## subsampler install subsampler
setwd("results/model_data/origin_infer/BA1/data/")
awk '{print $1}'  results/model_data/origin_infer/BA1/metadata.tsv |seqkit grep -f - results/model_data/origin_infer/BA1_rmdup.fasta  > sequences.fasta
python scripts/subsampler/get_daily_matrix_global.py --download yes --end-date 2022-03-15 --start-date 2021-09-15 
python3 scripts/subsampler/get_genome_matrix.py --metadata metadata.tsv --index-column  code --extra-columns country --date-column  date  --output  outputs/genome_matrix_days.tsv
python3 scripts/subsampler/aggregator.py --input outputs/genome_matrix_days.tsv --unit week --format integer --output outputs/matrix_genomes_unit.tsv
python3 scripts/subsampler/aggregator.py --input time_series_covid19_global_reformatted.tsv --unit week --format integer --start-date 2021-09-15 --output outputs/matrix_cases_unit.tsv
python3 scripts/subsampler/correct_bias.py --baseline 0.0001 --genome-matrix outputs/matrix_genomes_unit.tsv --case-matrix outputs/matrix_cases_unit.tsv --index-column code --output1 outputs/weekly_sampling_proportions.tsv --output2 outputs/weekly_sampling_bias.tsv --output3 outputs/matrix_genomes_unit_corrected.tsv
python3 scripts/subsampler/subsampler_timeseries.py  --metadata metadata.tsv   --genome-matrix outputs/matrix_genomes_unit_corrected.tsv --max-missing 99    --seed 2021  --index-column strain --geo-column code  --date-column date --time-unit week  --weekasdate no --start-date 2021-09-15 --end-date 2022-03-15   --sampled-sequences outputs/selected_sequences.txt --sampled-metadata outputs/selected_metadata.tsv --report outputs/sampling_stats.txt --keep ~/Project/04.SARS-COV-2/01.global_epi/00.data/BA2/data/config/keep.txt
seqkit grep -n -f outputs/selected_sequences.txt sequences.fasta  > sequences_select.fasta
sed -i 's/ //g' sequences_select.fasta 
## remove outliers
nextalign  run  -r data/inferred_origins/nextalign/reference.fasta -o sequences_select_align.fasta  sequences_select.fasta
awk -F "\t" '{print $1"\t"$3}' outputs/selected_metadata.tsv |sed 's/ //g' >  BA1_tempest_date.tsv
#remove outliers by tempest
seqkit grep -n -v -f data/tempest_remove sequences_select_align.fasta > sequences_select2_align.fasta
seqkit subseq -r 266:29674 sequences_select2_align.fasta > sequences_select2_align_codon.fasta
fasttree -nt -gtr  -boot 100 sequences_select2_align_codon.fasta  > BA1_step2.tree
treetime --tree BA1_step2.tree  --dates outputs/selected_metadata.tsv --aln sequences_select2_align_codon.fasta --outdir timetree --clock-rate 0.0008  --date-column date  --name-column strain  --confidence  --clock-std-dev 0.0004 --stochastic-resolve

# reference tree & annotation
treetime mugration --tree timetree/timetree.nexus  --states  outputs/selected_metadata.tsv --attribute code --confidence --outdir mugration --name-column strain
faToVcf -maskSites=data/inferred_origins/problematic_sites_sarsCov2_codon.vcf <(cat data/inferred_origins/wuhan01_codon.fasta sequences_select2_align_codon.fasta )  usher/sequences_select2_align_codon.vcf
head -7 mugration/annotated_tree.nexus |tail -1 |sed 's/^.\{12\}//g'>  mugration/annotated_tree.tree 
usher -t mugration/annotated_tree.tree  -v usher/sequences_select2_align_codon.vcf -o usher/BA2_usher_step1.pb -d usher
treetime mugration --tree usher/final-tree.nh --states outputs/selected_metadata.tsv --attribute code --confidence --outdir usher/mugration --name-column strain 

# add samples to usher tree
awk '{print $2}' |sort |uniq results/model_data/origin_infer/BA1/metadata.tsv > results/model_data/origin_infer/BA1/time
mkdir newsamples
setwd("results/model_data/origin_infer/BA1/data/newsamples/")

while read a; do
        awk -F "\t" '($2=="'$a'"){print $1}' results/model_data/origin_infer/BA1/metadata.tsv  |seqkit grep  -f - results/model_data/origin_infer/BA1/data/sequences.fasta > ${a}.fasta ;
       nextalign  run  -r data/inferred_origins/nextalign/reference.fasta -o  ${a}_align.fasta  ${a}.fasta        ;
       seqkit subseq -r 266:29674  ${a}_align.fasta >  ${a}_align_codon.fasta;
       cat data/inferred_origins/wuhan01_codon.fasta ${a}_align_codon.fasta > ${a}_align_codon2.fasta
       faToVcf -maskSites=data/inferred_origins/problematic_sites_sarsCov2_codon.vcf  ${a}_align_codon2.fasta ${a}_align_codon.vcf ;
      usher -i  results/model_data/origin_infer/BA1/usher/BA1_usher_step1.pb -v  ${a}_align_codon.vcf -p -d ${a} ;
    awk '($4=="y"){print $1"\t"$2}' ${a}/parsimony-scores.tsv |awk   '!a[$1]++'  |sort -t "," -k1   > ${a}/parsimony-scores2.tsv;
      rm ${a}*fasta ${a}*vcf ${a}*fai ${a}/parsimony-scores.tsv;

      done < results/model_data/origin_infer/BA1/time

# infer 1st 
Rscript r3.importe_auto.r ####
cat */parsimony-scores2_import_result.csv > import_result_step1.csv
Rscript r4.importe_less.r  import_result_step1.csv > import_result_step2.csv
awk '($11=="Imported"){print $1}' import_result_step2.csv > imported_ID
# infer 2nd 
seqkit grep -n -f imported_ID  results/model_data/origin_infer/BA1/data/sequences.fasta > imported.fasta
nextalign run -r data/inferred_origins/nextalign/reference.fasta -o imported_align.fasta imported.fasta
seqkit subseq -r 266:29674 imported_align.fasta > imported_align_codon.fasta
cat data/inferred_origins/wuhan01_codon.fasta imported_align_codon.fasta > imported_align_codon2.fasta
faToVcf -maskSites=data/inferred_origins/problematic_sites_sarsCov2_codon.vcf imported_align_codon2.fasta imported_align_codon2.vcf
usher -i results/model_data/origin_infer/BA1/usher/BA1_usher_step1.pb -v imported_align_codon2.vcf  -d usher2
cat results/model_data/origin_infer/BA1/data/outputs/selected_metadata.tsv  > imported_metadata.tsv
awk '($10=="Imported"){print $1}' all_imported.tsv |grep -f - results/model_data/origin_infer/BA1/metadata.tsv >> imported_metadata.tsv
treetime mugration --tree usher2/final-tree.nh --states outputs/selected_metadata.tsv --attribute code --confidence --outdir usher2/mugration2 --name-column strain
Rscript r5.importe_final.r
