conda activate /gpfs/home/bsc/bsc083001/miniconda3/envs/lr-kallisto
/gpfs/home/bsc/bsc083001/miniconda3/envs/lr-kallisto/bin/kallisto bus \
    -t 32  \
    --long  \
    --threshold 0.8 \
    -x bulk \
    -i ../../ref/HG002_maternal/HG002_maternal_k-63.idx \
    -o ../../data/personal_kallisto/HG002_maternal/GM18631_1/ \
    /gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/02_ONT_preprocessing/data/q7/15_CH3_GM18631_preprocessed_Q7.fastq.gz
