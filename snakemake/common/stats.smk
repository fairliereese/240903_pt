rule fq_count_reads:
    resources:
        nodes = 1,
        threads = 1
    shell:
        """
            n=`zcat {input.fq} | wc -l`
            n2=4
            echo $((${{n}} / ${{n2}})) > {output.txt}
        """

rule count_reads_summary:
    resources:
        nodes = 2,
        threads = 1
    run:
        df = pd.DataFrame()
        for f,d in zip(input.txts, params.samples):
            with open(f, 'r') as infile:
                for i, line in enumerate(infile):
                    if i == 0:
                        n = line.strip()

            temp = pd.DataFrame()
            temp['dataset'] = [d]
            temp['n_reads'] = [n]
            df = pd.concat([df, temp], axis=0)
            df.to_csv(output.summ, sep='\t', index=False)

rule count_bam:
    resources:
        threads = 8,
        nodes = 1
    shell:
        """
        module load samtools
        samtools view -c {input.align} > {output.txt}
        """

rule read_id_df_summary:
    resources:
        threads = 1,
        nodes = 2,
    run:
        df = pd.DataFrame()
        for f, d in zip(input.tsvs, params.samples):
            temp = pd.read_csv(f, sep='\t')
            temp['dataset'] = d
            df = pd.concat([df, temp], axis=0)
        df.to_csv(output.summ, sep='\t', index=False)

rule bam_get_mapqs:
    resources:
        threads = 1,
        nodes = 2
    shell:
        """
        module load samtools
        samtools view {input.bam} | grep -v ^@ | cut -f1,5 > {output.txt}
        """

rule bam_get_query_cov:
    resources:
        nodes = 1,
        threads = 16
    run:
        df = compute_query_coverage(input.align,
                                    resources.threads)
        df.to_csv(output.out, sep='\t', index=False)
