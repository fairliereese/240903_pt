def get_df_val(df, col1, col_dict):
    temp = df.copy(deep=True)

    for key, item in col_dict.items():
        temp = temp.loc[temp[key] == item]

    val = temp[col1].unique()
    assert len(val) == 1
    return val[0]


rule kallisto_build_ind:
    input:
        fa = config['ref']['fa'],
        gtf = config['lr']['novel_gene']['gtf_added_gene_entries']
    params:
        kallisto_path = '/gpfs/home/bsc/bsc083001/miniconda3/envs/lr-kallisto/bin/kallisto'
    output:
        ind = config['ref']['kallisto']['ind'],
        fa = config['ref']['kallisto']['t_fa'],
        t2g = config['ref']['kallisto']['t2g'],
    resources:
        threads = 8,
        nodes = 4
    conda:
        'base'
    shell:
        """
        conda activate /gpfs/home/bsc/bsc083001/miniconda3/envs/lr-kallisto
        kb ref \
            --kallisto {params.kallisto_path} \
            -i {output.ind} \
            -k 63 \
            -f1 {output.fa} \
            -g {output.t2g} \
            {input.fa} \
            {input.gtf}
        """

rule kallisto_pseudoalign:
    input:
        ind = config['ref']['kallisto']['ind'],
        fq = lambda wc: expand(config['lr']['fastq'],
                    lab_sample=get_df_val(df,
                            'lab_rep',
                            {'tech_rep': wc.sample})),
    params:
        kallisto_path = '/gpfs/home/bsc/bsc083001/miniconda3/envs/lr-kallisto/bin/kallisto',
        odir = config['lr']['kallisto']['pseudoalign']
    output:
        pseudoalign = directory(config['lr']['kallisto']['pseudoalign'])
    resources:
        threads = 32,
        nodes = 8
    conda:
        'base'
    shell:
        """
        conda activate /gpfs/home/bsc/bsc083001/miniconda3/envs/lr-kallisto
        {params.kallisto_path} bus \
            -t 32 \
            --long \
            --threshold 0.8 \
            -x bulk \
            -i {input.ind} \
            -o {params.odir} \
            {input.fq}
        """
# conda activate /gpfs/home/bsc/bsc083001/miniconda3/envs/bustools
# ${path_to_bustools} sort -t 32 ${output}/output.bus \
#  -o ${output}/sorted.bus; \
#  ${path_to_bustools} count ${output}/sorted.bus \
#  -t ${output}/transcripts.txt \
#  -e ${output}/matrix.ec \
#  -o ${output}/count --cm -m \
#  -g ${ref}.t2g; \
# ${path_to_lr_kallisto} quant-tcc -t 32 \
# 	--long -P ONT -f ${output}/flens.txt \
# 	${output}/count.mtx -i ${ref}_k-63.idx \
# 	-e ${output}/count.ec.txt \
# 	-o ${output};
