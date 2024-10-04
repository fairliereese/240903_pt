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

rule fq_get_read_ids:
    resources:
        nodes = 1,
        threads = 1
    shell:
        """
        grep '^@' {input.fq} | cut -d ' ' -f 1 > {output.txt}
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
        i = 0
        for f, d in zip(input.tsvs, params.samples):
            temp = pd.read_csv(f, sep='\t')
            temp['dataset'] = d
            if i == 0:
                temp.to_csv(output.summ, sep='\t', index=False)
            else:
                temp.to_csv(output.summ, sep='\t', index=False, header=None, mode='a')
            i+=1

rule bam_get_mapqs:
    resources:
        threads = 1,
        nodes = 2
    shell:
        """
        module load samtools
        echo -e 'read_id\tmapq'> {output.txt}
        samtools view {input.bam} | grep -v ^@ | cut -f1,5 >> {output.txt}
        """

rule bam_get_query_cov:
    resources:
        nodes = 1,
        threads = 16
    run:
        df = compute_query_coverage(input.align,
                                    resources.threads)
        df.to_csv(output.out, sep='\t', index=False)

##################
#### mapping summary statistics
##################

rule bool_mapq_summary:
    resources:
        nodes = 3,
        threads = 1
    run:
        files = list(input.files)
        assemblies = params.assemblies
        thresh = params.mapq_thresh

        # get the color (v important)
        sample_1 = wildcards.sample.split('_')[0]
        meta = load_meta()
        pop = meta.loc[meta['sample'] == sample_1, 'population'].values[0]
        c_dict, _ = get_population_colors()
        color = c_dict[pop]

        # get the upset plot table
        i = 0
        for f, a in zip(files, assemblies):
            temp = pd.read_csv(f, sep='\t')
            temp.rename({'mapq':a}, axis=1, inplace=True)
            assert len(temp.index) == len(temp.read_id.unique())

            if i == 0:
                df = temp.copy(deep=True)
            else:
                df = df.merge(temp, how='outer', on='read_id')
            i += 1

        # convert to binary
        df.fillna(0, inplace=True)
        df.set_index('read_id', inplace=True)
        df = df>thresh

        df.reset_index(inplace=True)
        df.set_index(assemblies, inplace=True)

        # make the upset plot
        ax_dict = upsetplot.UpSet(df, subset_size='count',
                          facecolor=color,
                          sort_by='cardinality',
                          show_counts=False,
                          show_percentages=True).plot()
        plt.suptitle(f'% of {sample} reads w/ mapq>{thresh}')
        plt.savefig(output.upset, dpi=500)

        # get the afr only reads
        afr_reads = df.copy(deep=True)
        afr_reads.reset_index(inplace=True)

        non_afr_assemblies = list(set(assemblies)-set(['afr']))

        afr_reads = afr_reads.loc[(afr_reads.afr==True)]
        for a in non_afr_assemblies:
            afr_reads = afr_reads.loc[afr_reads[a]==False]
        afr_reads = afr_reads[['read_id']]
        afr_reads.to_csv(output.afr_reads, index=False)

        # get the summary table
        df['total_reads'] = df.n_reads.sum()
        df['perc'] = (df.n_reads/df.total_reads)*100
        df['sample'] = sample
        df.to_csv(output.tsv, sep='\t', index=False)
