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

rule kallisto_t2g_8col:
    input:
        t2g = config['ref']['kallisto']['t2g']
    output:
        t2g = config['ref']['kallisto']['t2g_8col']
    resources:
        nodes = 1,
        threads = 1
    shell:
        """
        awk '{{print $1, $2, $3, $1, $4, $5, $6, $7}}' {input.t2g} > {output.t2g}
        """

rule kallisto_get_t2t:
    input:
        t2g = config['ref']['kallisto']['t2g_8col']
    output:
        t2t = config['ref']['kallisto']['t2t']
    resources:
        nodes = 1,
        threads = 1
    shell:
        """
        awk '{{print $1"\\t"$1"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7"\\t"$8}}' {input.t2g} > {output.t2t}
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

rule bustools_count_uniq:
    input:
        bus = config['lr']['kallisto']['bus_sort'],
        transcripts = config['lr']['kallisto']['transcripts'],
        matrix = config['lr']['kallisto']['matrix'],
        t2t = config['ref']['kallisto']['t2g_8col'] # TODO
    params:
        bustools_path = '/gpfs/home/bsc/bsc083001/miniconda3/envs/bustools/bin/bustools',
        count_pref = config['lr']['kallisto']['count_pref_uniq']
    output:
        mtx = config['lr']['kallisto']['count_mtx_uniq'],
        t = config['lr']['kallisto']['count_transcripts_uniq']
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
            --genecounts  \
            -g {input.t2t}
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
        kallisto_path = '/gpfs/home/bsc/bsc083001/miniconda3/envs/lr-kallisto/bin/kallisto',
        odir = config['lr']['kallisto']['quant']['odir']
    output:
        quant = config['lr']['kallisto']['quant']['matrix'],
        tpm = config['lr']['kallisto']['quant']['matrix_tpm'],
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

rule fmt_mtx_transcripts:
    resources:
        nodes = 2,
        threads = 1
    run:
        import scipy
        import numpy as np
        count = scipy.io.mmread(input.mtx)
        labels = pd.read_csv(input.ts, header=None, sep='\t')
        kallisto_df = pd.DataFrame(count.todense().T, columns=[params.col])
        kallisto_df['transcript_id'] = [labels.values[i][0] for i in range(np.shape(labels.values)[0])]
        kallisto_df.to_csv(output.tsv, sep="\t", columns=['transcript_id',params.col], header=1, index=0)

use rule fmt_mtx_transcripts as fmt_mtx_transcripts_counts with:
    input:
        mtx = config['lr']['kallisto']['quant']['matrix'],
        ts = config['lr']['kallisto']['transcripts']
    params:
        col = 'counts'
    output:
        tsv = config['lr']['kallisto']['quant']['matrix_tsv']

# use rule fmt_mtx_transcripts as fmt_mtx_transcripts_tpm with:
#     input:
#         mtx = config['lr']['kallisto']['quant']['matrix_tpm'],
#         ts = config['lr']['kallisto']['transcripts']
#     params:
#         col = 'counts'
#     output:
#         tsv = config['lr']['kallisto']['quant']['matrix_tpm_tsv']

use rule fmt_mtx_transcripts as fmt_mtx_transcripts_uniq with:
    input:
        mtx = config['lr']['kallisto']['count_mtx_uniq'],
        ts = config['lr']['kallisto']['count_transcripts_uniq']
    params:
        col = 'counts'
    output:
        tsv = config['lr']['kallisto']['matrix_tsv_uniq']
