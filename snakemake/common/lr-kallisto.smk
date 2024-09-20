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
        odir = config['lr']['kallisto']['odir']
    output:
        bus = config['lr']['kallisto']['bus'],
        flens = config['lr']['kallisto']['flens'],
        transcripts = config['lr']['kallisto']['transcripts'],
        matrix = config['lr']['kallisto']['matrix'],
    resources:
        threads = 32,
        nodes = 8
    conda:
        'base'
    shell:
        """
        conda activate /gpfs/home/bsc/bsc083001/miniconda3/envs/lr-kallisto
        {params.kallisto_path} bus \
            -t {resources.threads} \
            --long \
            --threshold 0.8 \
            -x bulk \
            -i {input.ind} \
            -o {params.odir} \
            {input.fq}
        """

rule bustools_sort:
    input:
        bus = config['lr']['kallisto']['bus']
    params:
        bustools_path = '/gpfs/home/bsc/bsc083001/miniconda3/envs/bustools/bin/bustools'
    output:
        bus = config['lr']['kallisto']['bus_sort']
    resources:
        threads = 32,
        nodes = 4
    conda:
        'base'
    shell:
        """
        conda activate /gpfs/home/bsc/bsc083001/miniconda3/envs/bustools
        {params.bustools_path} sort \
            -t {resources.threads} \
            {input.bus} \
            -o {output.bus}
        """

rule bustools_count:
    input:
        bus = config['lr']['kallisto']['bus_sort'],
        transcripts = config['lr']['kallisto']['transcripts'],
        matrix = config['lr']['kallisto']['matrix'],
        t2g = config['ref']['kallisto']['t2g']
    params:
        bustools_path = '/gpfs/home/bsc/bsc083001/miniconda3/envs/bustools/bin/bustools',
        count_pref = config['lr']['kallisto']['count_pref']
    output:
        mtx = config['lr']['kallisto']['count_mtx'],
        ec = config['lr']['kallisto']['count_ec']
    resources:
        threads = 32,
        nodes = 4
    conda:
        'base'
    shell:
        """
        conda activate /gpfs/home/bsc/bsc083001/miniconda3/envs/bustools
        {params.bustools_path} count \
            {input.bus} \
             -t {input.transcripts} \
             -e {input.matrix} \
             -o {params.count_pref} \
            --cm \
            -m \
            -g {input.t2g}
        """

rule lr_kallisto:
    input:
        flens = config['lr']['kallisto']['flens'],
        mtx = config['lr']['kallisto']['count_mtx'],
        ind = config['ref']['kallisto']['ind'],
        ec = config['lr']['kallisto']['count_ec']
    resources:
        threads = 32,
        nodes = 8
    conda:
        'base'
    params:
        kallisto_path = '/gpfs/home/bsc/bsc083001/miniconda3/envs/lr-kallisto/bin/kallisto',,
        odir = config['lr']['kallisto']['quant']['odir']
    output:
        quant = config['lr']['kallisto']['matrix']
    shell:
        """
        conda activate /gpfs/home/bsc/bsc083001/miniconda3/envs/lr-kallisto
        {params.kallisto_path} quant-tcc \
            -t {resources.threads} \
            --long \
            -P ONT \
            -f {input.flens} \
            {input.mtx} \
            -i {input.ind} \
            -e {input.ec} \
            -o {params.odir}
        """


# ${path_to_lr_kallisto} quant-tcc -t 32 \
# 	--long -P ONT -f ${output}/flens.txt \
# 	${output}/count.mtx -i ${ref}_k-63.idx \
# 	-e ${output}/count.ec.txt \
# 	-o ${output};
