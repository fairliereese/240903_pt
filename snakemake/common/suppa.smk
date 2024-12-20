suppa_events = ['A3', 'A5', 'AF', 'AL', 'MX', 'RI', 'SE']

# the input format is so dumb
rule fmt_kallisto_to_suppa_ab:
    resources:
        threads = 1,
        nodes = 2
    run:
        df = pd.read_csv(input.ab, sep='\t')
        df = df.set_index('transcript_id')
        df.index.name = ''
        header = '\t'+'\t'.join(df.columns.tolist())+'\n'
        with open(output.ab, 'w') as ofile:
            ofile.write(header)
        df.to_csv(output.ab, sep='\t', header=None, mode='a')

rule suppa_generate_events:
    resources:
        nodes = 1,
        threads = 1
    conda:
        'base'
    shell:
        """
        conda activate /gpfs/projects/bsc83/utils/conda_envs/suppa2
        python /gpfs/home/bsc/bsc083001/mele_lab/bin/SUPPA/suppa.py generateEvents \
            -i {input.gtf} \
            -o {params.opref} \
            -f ioe \
            -e SE SS MX RI FL
        """

rule suppa_psi:
    resources:
        nodes = 1,
        threads = 1
    conda:
        'base'
    shell:
        """
        conda activate /gpfs/projects/bsc83/utils/conda_envs/suppa2
        python /gpfs/home/bsc/bsc083001/mele_lab/bin/SUPPA/suppa.py psiPerEvent \
            --ioe-file {input.ioe} \
            --expression-file {input.filt_ab} \
            --o {params.opref}
        """
