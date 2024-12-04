def get_df_val(df, col1, col_dict):
    temp = df.copy(deep=True)

    for key, item in col_dict.items():
        temp = temp.loc[temp[key] == item]

    val = temp[col1].unique()
    assert len(val) == 1
    return val[0]


rule short_kallisto_build_ind:
    params:
        kallisto_path = '/gpfs/home/bsc/bsc083001/miniconda3/envs/kallisto/bin/kallisto',
    resources:
        threads = 8,
        nodes = 2
    conda:
        'base'
    shell:
        """
        conda activate /gpfs/home/bsc/bsc083001/miniconda3/envs/kallisto
        kb ref \
            --kallisto {params.kallisto_path} \
            -i {output.ind} \
            -k 31 \
            -f1 {output.fa} \
            -g {output.t2g} \
            {input.fa} \
            {input.gtf}
        """

rule short_kb_count:
    params:
        kallisto_path = '/gpfs/home/bsc/bsc083001/miniconda3/envs/kallisto/bin/kallisto',
    resources:
        threads = 8,
        nodes = 2
    conda:
        'base'
    shell:
        """
        conda activate /gpfs/home/bsc/bsc083001/miniconda3/envs/kallisto
        kb count \
            -x BULK \
            -o {params.odir} \
            -i {input.ind} \
            -g {input.t2g} \
            --parity=paired \
            --strand=unstranded \
            --tcc \
            --matrix-to-directories \
            {input.r1_fq} {input.r2_fq}
        """
# rule short_kallisto_quant:
#     params:
#         kallisto_path = '/gpfs/home/bsc/bsc083001/miniconda3/envs/kallisto/bin/kallisto',
#     resources:
#         threads = 32,
#         nodes = 2
#     conda:
#         'base'
#     shell:
#         """
#         conda activate /gpfs/home/bsc/bsc083001/miniconda3/envs/kallisto
#         {params.kallisto_path} bus \
#             -t {resources.threads} \
#             -x bulk \
#             -i {input.ind} \
#             -g {input.t2g} \
#             -o {params.odir} \
#             --paired \
#             {input.r1_fq} {input.r2_fq} 2> temp.err
#         """
#
# rule short_bustools_count:
#     params:
#         bustools_path = '/gpfs/home/bsc/bsc083001/miniconda3/envs/bustools/bin/bustools',
#         count_pref = config['lr']['kallisto']['count_pref']
#     resources:
#         threads = 32,
#         nodes = 2
#     conda:
#         'base'
#     shell:
#         """
#         conda activate /gpfs/home/bsc/bsc083001/miniconda3/envs/bustools
#         {params.bustools_path} count \
#             {input.bus} \
#              -t {input.transcripts} \
#              -e {input.matrix} \
#              -o {params.count_pref} \
#             --cm \
#             --multimapping \
#             -g {input.t2g}
#         """
#
# rule short_kallisto_quant:
#     resources:
#         threads = 32,
#         nodes = 2
#     conda:
#         'base'
#     params:
#         kallisto_path = '/gpfs/home/bsc/bsc083001/miniconda3/envs/kallisto/bin/kallisto',
#         odir = config['lr']['kallisto']['quant']['odir']
#     shell:
#         """
#         conda activate /gpfs/home/bsc/bsc083001/miniconda3/envs/kallisto
#         {params.kallisto_path} quant-tcc \
#             -t {resources.threads} \
#             -f {input.flens} \
#             {input.mtx} \
#             -i {input.ind} \
#             -e {input.ec} \
#             -g {input.t2g} \
#             -o {params.odir}
#         """
#
#
#
# #
# # rule short_fmt_mtx_transcripts:
# #     resources:
# #         nodes = 2,
# #         threads = 1
# #     run:
# #         import scipy
# #         import numpy as np
# #         count = scipy.io.mmread(input.mtx)
# #         labels = pd.read_csv(input.ts, header=None, sep='\t')
# #         kallisto_df = pd.DataFrame(count.todense().T, columns=[params.col])
# #         kallisto_df['transcript_id'] = [labels.values[i][0] for i in range(np.shape(labels.values)[0])]
# #         kallisto_df.to_csv(output.tsv, sep="\t", columns=['transcript_id',params.col], header=1, index=0)
# #
# # # use rule short_fmt_mtx_transcripts as fmt_mtx_transcripts_counts with:
# # #     input:
# # #         mtx = config['lr']['kallisto']['quant']['matrix'],
# # #         ts = config['lr']['kallisto']['transcripts']
# # #     params:
# # #         col = 'counts'
# # #     output:
# # #         tsv = config['lr']['kallisto']['quant']['matrix_tsv']
# #
# # # use rule short_fmt_mtx_transcripts as fmt_mtx_transcripts_tpm with:
# # #     input:
# # #         mtx = config['lr']['kallisto']['quant']['matrix_tpm'],
# # #         ts = config['lr']['kallisto']['transcripts']
# # #     params:
# # #         col = 'counts'
# # #     output:
# # #         tsv = config['lr']['kallisto']['quant']['matrix_tpm_tsv']
# #
# # # use rule short_fmt_mtx_transcripts as fmt_mtx_transcripts_uniq with:
# # #     input:
# # #         mtx = config['lr']['kallisto']['count_mtx_uniq'],
# # #         ts = config['lr']['kallisto']['count_transcripts_uniq']
# # #     params:
# # #         col = 'counts'
# # #     output:
# # #         tsv = config['lr']['kallisto']['matrix_tsv_uniq']
