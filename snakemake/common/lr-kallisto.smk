rule kallisto_build_ind:
    input:
        fa = config['ref']['fa'],
        gtf = config['lr']['gtf_chr_renamed']
    params:
        kallisto_path = '/gpfs/home/bsc/bsc083001/miniconda3/envs/lr-kallisto/bin/kallisto'
    output:
        ind = config['ref']['kallisto']['ind'],
        fa = config['ref']['kallisto']['t_fa'],
        t2g = config['ref']['kallisto']['t2g'],
    conda:
        'base'
    shell:
        """
        conda activate /gpfs/home/bsc/bsc083001/miniconda3/envs/lr-kallisto
        kb ref \
            --kallisto {parmas.kallisto_path} \
            -i {output.ind} \
            -k 63 \
            -f1 {output.fa} \
            -g {output.t2g} \
            {input.fa} \
            {input.gtf}
        """
