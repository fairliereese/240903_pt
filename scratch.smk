rule all:
    input:
        config['lr']['personal_kallisto']['quant']['merge_matrix_tpm_tsv'],
        config['lr']['personal_kallisto']['quant']['merge_matrix_tsv']

use rule kallisto_build_ind as kallisto_ind with:
    input:
        fa = config['ref']['alt']['personal']['fa_gz'],
        # gtf = config['lr']['gtf_filt_with_genes'] # TODO fabiens GTF
    output:
        ind = config['ref']['personal']['kallisto']['ind'],
        fa = config['ref']['personal']['kallisto']['t_fa'],
        t2g = config['ref']['personal']['kallisto']['t2g'],


use rule kallisto_pseudoalign as pseudoalign with:
    input:
        ind = config['ref']['personal']['kallisto']['ind'],
        fq = lambda wc: expand(config['lr']['fastq'],
                    lab_sample=get_df_val(df,
                            'lab_rep',
                            {'tech_rep': wc.sample})),
    params:
        kallisto_path = '/gpfs/home/bsc/bsc083001/miniconda3/envs/lr-kallisto/bin/kallisto',
        odir = config['lr']['personal_kallisto']['odir']
    output:
        bus = config['lr']['personal_kallisto']['bus'],
        flens = config['lr']['personal_kallisto']['flens'],
        transcripts = config['lr']['personal_kallisto']['transcripts'],
        matrix = config['lr']['personal_kallisto']['matrix']

use rule bustools_sort as bt_sort with:
    input:
        bus = config['lr']['personal_kallisto']['bus']
    output:
        bus = config['lr']['personal_kallisto']['bus_sort']

use rule bustools_count as bt_count with:
    input:
        bus = config['lr']['personal_kallisto']['bus_sort'],
        transcripts = config['lr']['personal_kallisto']['transcripts'],
        matrix = config['lr']['personal_kallisto']['matrix'],
        t2g = config['ref']['personal']['kallisto']['t2g']
    params:
        bustools_path = '/gpfs/home/bsc/bsc083001/miniconda3/envs/bustools/bin/bustools',
        count_pref = config['lr']['personal_kallisto']['count_pref']
    output:
        mtx = config['lr']['personal_kallisto']['count_mtx'],
        ec = config['lr']['personal_kallisto']['count_ec']

use rule lr_kallisto as run_lr_kallisto with:
    input:
        flens = config['lr']['personal_kallisto']['flens'],
        mtx = config['lr']['personal_kallisto']['count_mtx'],
        ind = config['ref']['personal']['kallisto']['ind'],
        ec = config['lr']['personal_kallisto']['count_ec']
    params:
        kallisto_path = '/gpfs/home/bsc/bsc083001/miniconda3/envs/lr-kallisto/bin/kallisto',
        odir = config['lr']['personal_kallisto']['quant']['odir']
    output:
        quant = config['lr']['personal_kallisto']['quant']['matrix'],
        tpm = config['lr']['personal_kallisto']['quant']['matrix_tpm'],

use rule fmt_mtx_transcripts as get_counts_mtx with:
    input:
        mtx = config['lr']['personal_kallisto']['quant']['matrix'],
        ts = config['lr']['personal_kallisto']['transcripts']
    params:
        col = 'counts'
    output:
        tsv = config['lr']['personal_kallisto']['quant']['matrix_tsv']

use rule fmt_mtx_transcripts as get_tpm_mtx with:
    input:
        mtx = config['lr']['personal_kallisto']['quant']['matrix_tpm'],
        ts = config['lr']['personal_kallisto']['transcripts']
    params:
        col = 'counts'
    output:
        tsv = config['lr']['personal_kallisto']['quant']['matrix_tpm_tsv']

rule merge_matrices:
    resources:
        nodes = 3,
        threads = 1
    run:
        merge_df = pd.DataFrame()
        samples = list(params.samples)
        i = 0
        for t, s in zip(list(input.tsvs), samples):
            temp = pd.read_csv(t, sep='\t')
            temp.rename({'counts': s}, axis=1, inplace=True)
            if i == 0:
                merge_df = temp.copy(deep=True)
            else:
                merge_df = merge_df.merge(temp, how='outer',
                            on='transcript_id')
            i+=1
        merge_df.to_csv(output.tsv, sep='\t', index=False)

use rule merge_matrices as merge_matrices_tpm with:
    input:
        tsvs = expand(config['lr']['personal_kallisto']['quant']['matrix_tpm_tsv'],
                      sample=df.tech_rep.tolist())
    params:
        samples = df.tech_rep.tolist()
    output:
        tsv = config['lr']['personal_kallisto']['quant']['merge_matrix_tpm_tsv']

use rule merge_matrices as merge_matrices_counts with:
    input:
        tsvs = expand(config['lr']['personal_kallisto']['quant']['matrix_tsv'],
                      sample=df.tech_rep.tolist())
    params:
        samples = df.tech_rep.tolist()
    output:
        tsv = config['lr']['personal_kallisto']['quant']['merge_matrix_tsv']
