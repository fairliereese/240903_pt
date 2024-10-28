conda activate /gpfs/home/bsc/bsc083001/miniconda3/envs/lr-kallisto
/gpfs/home/bsc/bsc083001/miniconda3/envs/lr-kallisto/bin/kallisto bus \
    -t 32  \
    --long  \
    --threshold 0.8 \
    -x bulk \
    -i ../../ref/HG002_maternal/HG002_maternal_k-63.idx \
    -o ../../data/personal_kallisto/HG002_maternal/GM18631_1/ \
    /gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/02_ONT_preprocessing/data/q7/15_CH3_GM18631_preprocessed_Q7.fastq.gz

shell:
 """
 conda activate /gpfs/projects/bsc83/utils/conda_envs/SQANTI3-5.2.1
 python /gpfs/projects/bsc83/Projects/pantranscriptome/fairlie/240903_pt/scripts//sqanti3_protein.py \
        ../../data/poder_protein/transcripts.transcript_exons_only.gtf \
        ../../data/poder_protein/transcripts.cds_renamed_exon.gtf \
        ../../data/poder_protein/transcripts_best_orf.tsv \
        ../../data/poder_protein/transcripts_gencode.transcript_exons_only.gtf \
        ../../data/poder_protein/transcripts_gencode.cds_renamed_exon.gtf \
        d ./                     -p ../../data/poder_protein/transcripts
"""

"""
python /gpfs/projects/bsc83/Projects/pantranscriptome/fairlie/240903_pt/scripts//create_protein_overview_table.py \
 --best_orf_path ../../data/poder_protein/transcripts_best_orf.tsv \
    --sqanti_protein_path ../../data/poder_protein/transcripts.sqanti_protein_classification.tsv \
    --orf_completeness_path ../../data/poder_protein/transcripts_ORF_completeness.tsv \
    --output_name ../../data/poder_protein/transcripts_protein_annotation.tsv \
    --gtf_original_path ../../data/transcripts_filt_with_gene.gtf \
    --gtf_predicted_path ../../data/poder_protein/transcripts_protein.gtf \
    --protein_fasta_path ../../data/poder_protein/transcripts_protein.fa \
    --blastp_path ../../data/poder_protein/transcripts_blastp.out 2> wha_happen.txt
"""
