import pandas as pd
from pyfaidx import Fasta
import pyranges as pr

def get_lr_encid(wc, df, kind):
    """
    Get the ENCID of the alignments file
    given the dataset and species name
    """
    dataset = wc.dataset
    temp = df.loc[(df.dataset==dataset)]
    if kind == 'filtered_alignments':
        c = 'ENCODE_alignments_id'
    elif kind == 'unfiltered_alignments':
        c = 'ENCODE_unfiltered_alignments_id'
    elif kind == 'reads':
        c = 'ENCODE_reads_id'
    return temp[c].values[0]

def parse_wgs_meta(meta):
    df = pd.read_csv(meta, sep='\t')
    df['paired_acc'] = df['Paired with'].str.split('/', expand=True)[2]
    df = df.loc[df['Paired end'] == 1]

    # get bio rep numbers -- each exp w/ same biosamp is bio rep
    temp = df[['Experiment accession', 'Biosample term name']].drop_duplicates()
    temp['biorep'] = temp.sort_values(['Biosample term name',
                             'Experiment accession'],
                              ascending=[True, True])\
                              .groupby(['Biosample term name']) \
                              .cumcount() + 1
    temp.sort_values(by=['Biosample term name', 'biorep'])
    df = df.merge(temp, how='left', on=['Biosample term name', 'Experiment accession'])

    # get tech rep numbers -- each file name w/ same exp is a tech rep
    temp = df[['Experiment accession', 'File accession']]
    temp['techrep'] = temp.sort_values(['File accession',
                             'Experiment accession'],
                              ascending=[True, True])\
                              .groupby(['Experiment accession']) \
                              .cumcount() + 1
    temp.sort_values(by=['Experiment accession', 'techrep'])
    df = df.merge(temp, how='left', on=['Experiment accession', 'File accession'])

    print(df[['Biosample term name', 'Experiment accession', 'File accession', 'biorep', 'techrep', 'paired_acc']].sort_values(by=['Biosample term name', 'biorep', 'techrep']))

    df = df[['Biosample term name', 'Experiment accession', 'File accession', 'biorep', 'techrep', 'paired_acc']]

    # sample is biosample + biorep + tech rep
    df['biosamp'] =  df['Biosample term name'].str.lower().str.replace(' ', '_')
    df['sample'] = df.biosamp+'_'+\
                   df.biorep.astype(str)+'_'+\
                   df.techrep.astype(str)

    return df

def get_ena_query(cell_line_id, seq_len, platform, lib_strat):

    fields = 'run_accession%2Cexperiment_title%2Clibrary_source%2Csample_title%2Cstudy_accession%2Csample_accession%2Cexperiment_accession%2Cscientific_name%2Cinstrument_model%2Cinstrument_platform%2Cfastq_ftp%2Clibrary_selection%2Clibrary_strategy%2Cdescription&format=tsv'

    def raise_wat_error():
        print(f'cell line: {cell_line_id}, seq_len: {seq_len}, platform: {platform}, lib_strat: {lib_strat}')
        raise ValueError('wat is this')

    # get query platform
    if seq_len == 'sr':
        if platform == 'ill':
            instrument_platform = 'illumina'
        else: raise_wat_error()
    elif seq_len == 'lr':
        if platform == 'pb':
            instrument_platform = 'pacbio_smrt'
        elif platform == 'ont':
            instrument_platform = 'oxford_nanopore'
        else: raise_wat_error()

    # get the library prep protocol
    if lib_strat == 'wgs':
        library_source = 'genomic'
        library_strategy = 'WGS'
    else: raise_wat_error()

    query = f"""https://www.ebi.ac.uk/ena/portal/api/search?query=instrument_platform%3D%22{instrument_platform}%22%20AND%20library_source%3D%22{library_source}%22%20AND%20library_strategy%3D%22{library_strategy}%22%20AND%20(library_name%20%3D%20%22*{cell_line_id}*%22%20OR%20sample_title%20%3D%20%22*{cell_line_id}*%22)&result=read_run&fields={fields}"""
    return query
    # WGS

        # https://www.ebi.ac.uk/ena/portal/api/search?query=instrument_platform%3D%22illumina%22%20AND%20library_source%3D%22genomic%22%20AND%20library_strategy%3D%22WGS%22%20AND%20(library_name%20%3D%20%22NA18906%22%20OR%20sample_title%20%3D%20%22NA18906%22)&result=read_run&fields=run_accession,experiment_title,library_source,sample_title&limit=0&download=true&format=tsv

def make_sample_transcript(fa, n_reads, tsv):
    """
    Given a set of sequences and a target set of reads,
    get a "sample.transcript" file that pbsim needs
    """
    # Load the FASTA file
    fasta = Fasta(fa)

    # Initialize lists to store chromosome names and sequences
    chromosomes = []
    sequences = []

    # Iterate through each entry in the FASTA file
    for name in fasta.keys():
        chromosomes.append(name)
        sequences.append(str(fasta[name]))

    # Create a pandas DataFrame
    df = pd.DataFrame({
        'Chromosome': chromosomes,
        'Sequence': sequences
    })

    df['morph'] = df.Chromosome.str.split('.', expand=True, n=1)[1]
    assert len(df.Sequence.unique()) == len(df.morph.unique())

    # get the number of copies per morph
    # and format it in the way that pbsim wants
    df = df[['morph', 'Chromosome', 'Sequence']].groupby(['morph', 'Sequence']).count().rename({'Chromosome':'cn'}, axis=1).reset_index()
    df.rename({'morph': 'transcript',
               'cn': 'reads'}, axis=1,
              inplace=True)
    df['antisense_reads'] = 0
    df = df[['transcript', 'reads', 'antisense_reads', 'Sequence']]
    n_copies = df.reads.sum()
    scale = n_reads/n_copies
    df['reads'] = (df.reads*scale).astype(int)

    # should I artificially make some morphs w/ really low exp?
    # like by doing a power? or something? ie prop * prop and then multply
    # by some other scaling factor?
    df.to_csv(tsv, index=False, header=None, sep='\t')

def get_odir_from_fname(fname):
    """
    For programs that want an output directory.
    Return the directory name of the output file
    """
    return fname.rsplit('/', maxsplit=1)[0]+'/'

def get_odir_and_pref_from_fname(fname, sep='_', maxsplit=1):
    odir = get_odir_from_fname(fname)
    pref = fname.rsplit('/', maxsplit=1)[-1].split(sep, maxsplit=maxsplit)[0]
    return odir+pref

def parse_config(fname):
    df = pd.read_csv(fname, sep='\t')
    df['sample'] = df.lab_rep.str.rsplit('_', n=1, expand=True)[1]

    # get tech rep id -- each lab replicate w/ the same sample is a tech rep
    temp = df.drop_duplicates()
    temp['tech_rep_num'] = temp.sort_values(['sample',
                             'lab_rep'],
                              ascending=[True, True])\
                              .groupby(['sample']) \
                              .cumcount() + 1
    temp.sort_values(by=['sample', 'tech_rep_num'])
    df = df.merge(temp, how='left', on=['sample', 'lab_rep'])

    df['tech_rep'] = df['sample']+'_'+df['tech_rep_num'].astype(str)

    return df

def fmt_list_for_cli(l, sep=','):
    return sep.join(l)

def get_gtf_ss_bed(gtf_file, extension=0):
    """
    Get SSs from input GTF file. Extend based on
    user input
    """
    # read in ref. gtf and get SJs from there
    temp = pr.read_gtf(gtf_file)

    # get the introns
    temp = temp.features.introns()
    temp = temp.df

    temp.Start = temp.Start.astype(int)
    temp.End = temp.End.astype(int)
    assert len(temp.loc[temp.Start>temp.End])==0

    # melt to 5' and 3'
    temp = temp.melt(id_vars=['Chromosome', 'Strand'],
                     value_vars=['Start', 'End'])
    temp['sj_loc'] = ''
    temp.loc[temp.variable=='Start', 'sj_loc'] = 'start'
    temp.loc[temp.variable=='End', 'sj_loc'] = 'end'

    temp.rename({'value':'Start'}, axis=1, inplace=True)

    temp['Start'] = temp.Start-2
    temp['End'] = temp.Start+1

    # verified
    temp.loc[temp.sj_loc=='start', 'Start'] = temp.loc[temp.sj_loc=='start', 'Start']+1
    temp.loc[temp.sj_loc=='start', 'End'] = temp.loc[temp.sj_loc=='start', 'End']+1


    temp.loc[temp.sj_loc=='end', 'End'] = temp.loc[temp.sj_loc=='end', 'End']+2
    temp.loc[temp.sj_loc=='end', 'Start'] = temp.loc[temp.sj_loc=='end', 'Start']+2

    # extend based on user input
    temp = pr.PyRanges(temp)
    temp = temp.extend(extension)

    return temp

def get_internal_exon_info(gtf, ref_gtf):
    """
    Get novelty and eid<->tid info for each exon
    """

    def make_exon_id(df):
        """
        Get an ID for each unique exon using
        Chromosome, Strand, Start, End
        """
        df['eid'] = df['Chromosome'].astype(str)+'_'+\
                    df['Strand'].astype(str)+'_'+\
                    df['Start'].astype(str)+'_'+\
                    df['End'].astype(str)
        return df

    # PODER
    df = pr.read_gtf(gtf).df
    df = df.loc[df.Feature=='exon']

    # remove first and last exons
    df['first_or_last_exon'] = (~df.duplicated('transcript_id', keep='first'))|\
                               (~df.duplicated('transcript_id', keep='last'))
    df = df.loc[~df.first_or_last_exon]
    df = make_exon_id(df)

    # only keep unique exons
    df = df[['Chromosome', 'Strand', 'Start', 'End', 'eid']].drop_duplicates(keep='first')

    # GENCODE
    ref_df = pr.read_gtf(ref_gtf).df
    ref_df = ref_df.loc[ref_df.Feature=='exon']

    # remove first and last exons
    ref_df['first_or_last_exon'] = (~ref_df.duplicated('transcript_id', keep='first'))|\
                               (~ref_df.duplicated('transcript_id', keep='last'))
    ref_df = ref_df.loc[~ref_df.first_or_last_exon]
    ref_df = make_exon_id(ref_df)

    # only keep unique exons
    ref_df = ref_df[['Chromosome', 'Strand', 'Start', 'End', 'eid']].drop_duplicates(keep='first')

    # do a pyranges left join w/ PODER
    df = pr.PyRanges(df)
    ref_df = pr.PyRanges(ref_df)
    df = df.join(ref_df,
                 strandedness='same',
                 how='left',
                 report_overlap=True,
                 slack=0, suffix='_gc')

    # get exon end_site_novelty
    df = df.df
    temp = df[['eid']].drop_duplicates()
    temp['novelty'] = 'Novel'
    temp.loc[temp.eid.isin(df.loc[df.Overlap>=1].eid.tolist()), 'novelty'] = "Novel 5'/3'"
    temp.loc[temp.eid.isin(df.eid_gc.tolist()), 'novelty'] = 'Known'

    # now get eid<->tid mapping
    df = pr.read_gtf(gtf).df
    df = df.loc[df.Feature=='exon']

    # remove first and last exons
    df['first_or_last_exon'] = (~df.duplicated('transcript_id', keep='first'))|\
                               (~df.duplicated('transcript_id', keep='last'))
    df = df.loc[~df.first_or_last_exon]
    df = make_exon_id(df)
    df = df[['eid', 'transcript_id']]

    df = df.merge(temp,
                  how='inner',
                  on='eid')

    return df

def make_exon_id(df):
    """
    Get an ID for each unique exon using
    Chromosome, Strand, Start, End
    """
    df['eid'] = df['Chromosome'].astype(str)+'_'+\
                df['Strand'].astype(str)+'_'+\
                df['Start'].astype(str)+'_'+\
                df['End'].astype(str)
    return df

def get_novel_part_exons(ref_gtf, gtf,
                         exon_info_file):

    # PODER
    df = pr.read_gtf(gtf).df
    df = df.loc[df.Feature=='exon']

    # remove first and last exons
    df['first_or_last_exon'] = (~df.duplicated('transcript_id', keep='first'))|\
                               (~df.duplicated('transcript_id', keep='last'))
    df = df.loc[~df.first_or_last_exon]
    df = make_exon_id(df)

    # only keep unique exons
    df = df[['Chromosome', 'Strand', 'Start', 'End', 'eid']].drop_duplicates(keep='first')

    # GENCODE
    ref_df = pr.read_gtf(ref_gtf).df
    ref_df = ref_df.loc[ref_df.Feature=='exon']

    # remove first and last exons
    ref_df['first_or_last_exon'] = (~ref_df.duplicated('transcript_id', keep='first'))|\
                               (~ref_df.duplicated('transcript_id', keep='last'))
    ref_df = ref_df.loc[~ref_df.first_or_last_exon]
    ref_df = make_exon_id(ref_df)

    # only keep unique exons
    ref_df = ref_df[['Chromosome', 'Strand', 'Start', 'End', 'eid']].drop_duplicates(keep='first')

    # get uniq genomic ranges that are in poder compared to gencode
    df = pr.PyRanges(df)
    ref_df = pr.PyRanges(ref_df)

    df = df.subtract(ref_df,
                  strandedness='same')

    df = df.df
    df[['Chromosome_og', 'Strand_og', 'Start_og', 'End_og']] = df.eid.str.split('_', expand=True)
    df.Start_og = df.Start_og.astype(int)
    df.End_og = df.End_og.astype(int)

    exon_info = pd.read_csv(exon_info_file, sep='\t')
    exon_info = exon_info.drop('transcript_id', axis=1)
    exon_info = exon_info.drop_duplicates()
    df = df.merge(exon_info, how='left', on='eid')

    assert len(df.loc[(df.novelty=='Novel')].index) == len(df.loc[(df.Start==df.Start_og)&(df.End==df.End_og)].index)
    df = df[['Chromosome', 'Strand', 'Start', 'End', 'eid', 'novelty']]
    return df
