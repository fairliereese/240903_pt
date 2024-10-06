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
def rm_color_cats(palette, order, cats):
    if cats:
        keys = palette.keys()
        pop_list = []
        for key in keys:
            if key not in cats:
                pop_list.append(key)
        for p in pop_list:
            del palette[p]
        order = [o for o in order if o in cats]
    return palette, order

def get_population_colors(cats=None):
    palette = {'ITU': '#db72f2',
                 'PEL': '#ff3a33',
                 'HAC': '#4cb33e',
                 'AJI': '#46bff0',
                 'LWK': '#A09136',
                 'YRI': '#DFBD00',
                 'CEU': '#347eed',
                 'MPC': '#eb9d0c'}

    order = list(palette.keys())
    order.sort()

    palette, order = rm_color_cats(palette, order, cats)
    return palette, order

rule bool_mapq_summary:
    resources:
        nodes = 3,
        threads = 1
    run:
        import upsetplot
        import matplotlib.pyplot as plt
        files = list(input.files)
        assemblies = params.assemblies
        thresh = float(params.mapq_thresh)
        sample = wildcards.sample

        # get the color (v important)
        sample_1 = sample.split('_')[0]
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
        ax_dict = upsetplot.UpSet(df, subset_size='count', facecolor=color, sort_by='cardinality', show_counts=False, show_percentages=True).plot()
        # ax_dict = upsetplot.UpSet(df, subset_size='count',
        #                   facecolor=color,
        #                   sort_by='cardinality',
        #                   show_counts=False,
        #                   show_percentages=True).plot()
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
        df.reset_index(inplace=True)
        df = df.groupby(assemblies).count().reset_index().rename({'read_id':'n_reads'},axis=1)
        df['total_reads'] = df.n_reads.sum()
        df['perc'] = (df.n_reads/df.total_reads)*100
        df['sample'] = sample
        df.to_csv(output.tsv, sep='\t', index=False)

rule max_mapq_summary:
    resources:
        nodes = 3,
        threads = 1
    run:
        import upsetplot
        import matplotlib.pyplot as plt
        files = list(input.files)
        assemblies = params.assemblies
        thresh = float(params.mapq_thresh)
        sample = wildcards.sample

        # get the color (v important)
        sample_1 = sample.split('_')[0]
        meta = load_meta()
        pop = meta.loc[meta['sample'] == sample_1, 'population'].values[0]
        c_dict, _ = get_population_colors()
        color = c_dict[pop]

        # get upset plot table
        i = 0
        for f, a in zip(files, assemblies):
            temp = pd.read_csv(f, sep='\t')
            # temp.rename({'mapq':a}, axis=1, inplace=True)
            assert len(temp.index) == len(temp.read_id.unique())
            temp['assembly']=a
            if i == 0:
                df = temp.copy(deep=True)
            else:
                df = pd.concat([df, temp], axis=0)
            i += 1

        # assert min mapq
        df = df.loc[df.mapq>thresh]

        # groupby read id and mapq to find reads that map equally as well
        df = df.groupby(['read_id', 'mapq']).agg({
            'assembly': lambda x: ','.join(x)})

        import pdb; pdb.set_trace()
        # sort by mapq and dedupe by keeping max
        df.reset_index(inplace=True)
        df = df.sort_values(by='mapq', ascending=False)
        df = df.drop_duplicates(subset=['read_id'], keep='first')

        # process out using the upset plot stuff
        df['assembly'] = df.assembly.str.split(',')
        df = df.explode('assembly')
        df['val'] = True
        df = df.pivot(index='read_id', columns='assembly', values='val')
        df = df.fillna(False)

        # make the upset plot
        ax_dict = upsetplot.UpSet(df, subset_size='count', facecolor=color, sort_by='cardinality', show_counts=False, show_percentages=True).plot()
        # ax_dict = upsetplot.UpSet(df, subset_size='count',
        #                   facecolor=color,
        #                   sort_by='cardinality',
        #                   show_counts=False,
        #                   show_percentages=True).plot()
        plt.suptitle(f'% of best-mapping {sample} reads w/ mapq>{thresh}')
        plt.savefig(output.upset, dpi=500)


        # TODO

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
        df.reset_index(inplace=True)
        df = df.groupby(assemblies).count().reset_index().rename({'read_id':'n_reads'},axis=1)
        df['total_reads'] = df.n_reads.sum()
        df['perc'] = (df.n_reads/df.total_reads)*100
        df['sample'] = sample
        df.to_csv(output.tsv, sep='\t', index=False)
