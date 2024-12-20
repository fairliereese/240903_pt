suppa_events = ['A3', 'A5', 'AF', 'AL', 'MX', 'RI', 'SE']

rule suppa_generate_events:
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

use rule suppa_psi as suppa_get_psi with:
    input:
        ioe = lambda wc: config['lr']['suppa']['events'][wc.event],
        filt_ab = config['lr']['suppa']['fmt_ab']
    params:
        opref = config['lr']['suppa']['psi'].rsplit('.psi', maxsplit=1)[0]
    output:
        config['lr']['suppa']['psi']
