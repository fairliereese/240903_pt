import yaml
import os
import pandas as pd
import re


def load_config(config_file=None):
    """
    Load snakemake config file as a dictionary
    """
    if not config_file:
        d = os.path.dirname(__file__)
        od = f'{d}/../snakemake/'
        config_file = f'{od}/config.yml'

    with open(config_file) as f:
        config = yaml.safe_load(f)

    return config

def proc_cfg(entry, od):
    entry = entry.replace('../../', '')
    entry = od+entry
    return entry

def load_meta():
    """
    Load metadata file from config
    """
    d = os.path.dirname(__file__)
    od = f'{d}/../'

    config = load_config()
    df = pd.read_csv(proc_cfg(config['lr']['meta'], od), sep='\t')
    return df

def set_col_order(df, col, order):
    """
    Reorders the values in a specified column of a DataFrame according to a given list.

    Parameters:
    df (pd.DataFrame): The DataFrame containing the column to reorder.
    col (str): The name of the column to reorder.
    order (list): The list defining the order of the values in the column.

    Returns:
    pd.DataFrame: A DataFrame with the reordered column.
    """
    df[col] = pd.Categorical(df[col], categories=order, ordered=True)
    df = df.sort_values(col)
    return df

def get_rdna_regions():
    """
    Get string names of different rDNA regions
    """
    return ['5_ETS', '18S', 'ITS1',
            '5.8S', 'ITS2', '28S', '3_ETS']


def get_rdna_regions_genic_status():
    """
    Return a dict mapping rdna region names to whether
    they're genic or not
    """
    return {'5_ETS': 'non-genic',
            '18S': 'genic',
            'ITS1': 'non-genic',
            '5.8S': 'genic',
            'ITS2': 'non-genic',
            '28S': 'genic',
            '3_ETS': 'non-genic'}


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

    # print(df[['Biosample term name', 'Experiment accession', 'File accession', 'biorep', 'techrep', 'paired_acc']].sort_values(by=['Biosample term name', 'biorep', 'techrep']))

    df = df[['Biosample term name', 'Experiment accession', 'File accession', 'biorep', 'techrep', 'paired_acc']]

    # sample is biosample + biorep + tech rep
    df['biosamp'] =  df['Biosample term name'].str.lower().str.replace(' ', '_')
    df['sample'] = df.biosamp+'_'+\
                   df.biorep.astype(str)+'_'+\
                   df.techrep.astype(str)

    return df

def get_mapq(fname, threads):
    import pysam
    if fname.endswith('.bam'):
        in_mode = 'rb'
    else:
        in_mode = 'r'
    input =  pysam.AlignmentFile(fname, in_mode, threads=threads)
    read_ids = []
    read_mapqs = []
    for read in input:
        read_mapqs.append(read.mapping_quality)
        read_ids.append(read.query_name)
    df = pd.DataFrame()
    df['read_id'] = read_ids
    df['mapq'] = read_mapqs
    return df


def compute_query_coverage(fname, threads):
    import pysam
    if fname.endswith('.bam'):
        in_mode = 'rb'
    else:
        in_mode = 'r'
    input =  pysam.AlignmentFile(fname, in_mode, threads=threads)

    read_covs = []
    read_lens = []
    read_ids = []
    read_chrs = []
    for read in input:
        cigar = read.cigarstring
        if cigar:
            read_len, read_cov = compute_alignment_coverage(cigar)
            read_lens.append(read_len)
            read_covs.append(read_cov)
            read_ids.append(read.query_name)
            read_chrs.append(read.reference_name)
    df = pd.DataFrame()
    df['read_id'] = read_ids
    df['query_cov'] = read_covs
    df['read_len'] = read_lens
    df['chr'] = read_chrs
    return df

def compute_alignment_coverage(CIGAR):
    """This function computes what fraction of the read is actually aligned to
    the genome by excluding hard or soft-clipped bases."""

    total_bases = 0.0
    unaligned_bases = 0.0
    ops, counts = split_cigar(CIGAR)
    for op, ct in zip(ops, counts):
        if op == "N":
            continue
        if op == "H" or op == "S":
            unaligned_bases += ct
        total_bases += ct

    return (total_bases, ((total_bases - unaligned_bases) / total_bases))

def split_cigar(cigar):
    """Takes CIGAR string from SAM and splits it into two lists:
    one with capital letters (match operators), and one with
    the number of bases that each operation applies to."""

    alignTypes = re.sub("[0-9]", " ", cigar).split()
    counts = re.sub("[=A-Z]", " ", cigar).split()
    counts = [int(i) for i in counts]

    return alignTypes, counts

def get_dupe_read_names(fname, threads, out):
    import pysam
    if fname.endswith('.bam'):
        in_mode = 'rb'
    else:
        in_mode = 'r'
    input =  pysam.AlignmentFile(fname, in_mode, threads=threads)
    read_covs = []
    read_ids = []
    read_lens = []
    for read in input:
        cigar = read.cigarstring
        if cigar:
            read_len, read_cov = compute_alignment_coverage(cigar)
        else:
            read_len = read_cov = np.nan
        read_lens.append(read_len)
        read_covs.append(read_cov)
        read_ids.append(read.query_name)
    df = pd.DataFrame()
    df['read_id'] = read_ids
    df['query_cov'] = read_covs
    df['read_len'] = read_lens

    df = df.loc[df.read_id.duplicated(keep='first')]
    df.to_csv(out, sep='\t', index=False)

def tiebreak_supp_reads(fname, threads, output):
    # record the read id and query coverage for each read.
    # tiebreak duplicates w/ the higher coverage.
    import pysam
    if fname.endswith('.bam'):
        in_mode = 'rb'
    else:
        in_mode = 'r'
    input =  pysam.AlignmentFile(fname, in_mode, threads=threads)
    read_covs = []
    read_ids = []
    read_lens = []
    for read in input:
        cigar = read.cigarstring
        if cigar:
            read_len, read_cov = compute_alignment_coverage(cigar)
        else:
            read_len = read_cov = np.nan
        read_lens.append(read_len)
        read_covs.append(read_cov)
        read_ids.append(read.query_name)
    df = pd.DataFrame()
    df['read_id'] = read_ids
    df['query_cov'] = read_covs
    df['read_len'] = read_lens

    # drop the reads that are duplicated that have the lower coverage
    df['num_ind'] = [i for i in range(len(df.index))]
    df = df.sort_values(by='query_cov', ascending=False)
    df = df.drop_duplicates(subset='read_id', keep='first')

    if fname.endswith('.bam'):
        out_mode = 'wb'
    else:
        out_mode = 'w'
    input.close()

    input =  pysam.AlignmentFile(fname, in_mode, threads=threads)
    outfile = pysam.AlignmentFile(output, out_mode, template=input)
    # outfile = pysam.AlignmentFile(output, out_mode)
    for i, read in enumerate(input):
        if i in df.num_ind.tolist():
            outfile.write(read)

def get_vcf_summary(all_vcf,
                   filt_vcf,
                   call_vcf,
                   uncall_vcf,
                   ofile,
                   **kwargs):

    # get long form table for all first,
    # then label allele status as filtered,
    # unfiltered, callable, or uncallable
    df = get_hr_vcf(read_vcf(all_vcf), **kwargs)

    # passed filter?
    temp = get_hr_vcf(read_vcf(filt_vcf), **kwargs)
    ids = temp.var_id.unique().tolist()
    df['filt_pass'] = False
    df.loc[df.var_id.isin(ids), 'filt_pass'] = True

    # callable vs. uncallable
    temp = get_hr_vcf(read_vcf(call_vcf), **kwargs)
    ids = temp.var_id.unique().tolist()
    df['callable'] = np.nan
    df.loc[df.var_id.isin(ids), 'callable'] = True

    temp = get_hr_vcf(read_vcf(uncall), **kwargs)
    df.loc[df.var_id.isin(ids), 'callable'] = False

    df.to_csv(ofile, sep='\t', index=False)

def get_seq_lens():
    return ['sr', 'lr']
def get_platforms():
    return ['ill', 'pb', 'ont']
def get_lib_strats():
    return ['wgs', 'rna']

def get_valid_seq_len_platforms(seq_lens=get_seq_lens(),
                                platforms=get_platforms()):
    d = {'ill': 'sr',
         'pb': 'lr',
         'ont': 'lr'}
    d_pairs = [(key, item) for key, item in d.items()]
    pairs = [(p,s) for p in platforms for s in seq_lens if (p,s) in d_pairs]

    df = pd.DataFrame()
    df['platform'] = [p[0] for p in pairs]
    df['seq_len'] = [p[1] for p in pairs]
    return df

def rm_sirv_ercc_gtf(ifile, ofile):
    """
    Remove SIRV and ERCC entries from a GTF
    """
    import pyranges as pr
    df = pr.read_gtf(ifile).df
    df = df.loc[~df.Chromosome.str.contains('SIRV')]
    df = df.loc[~df.Chromosome.str.contains('ERCC')]
    df = pr.PyRanges(df)
    df.to_gtf(ofile)

def write_parsed_hmmer(ifile, ofile):
    from Bio import SearchIO
    
    infile = open(ifile, 'r')
    outfile = open(ofile, 'w')

    # write header to file first
    header = ['transcript_id', 'accession', 'bias',
              'bitscore', 'description',
              'evalue', 'id']
    outfile.write('\t'.join(header)+'\n')
    for record in SearchIO.parse(infile, 'hmmscan3-domtab'):
        tid = record.id
        for hit in record.hits:
            new_line = []
            new_line.append(tid)
            new_line.append(str(hit.accession))
            new_line.append(str(hit.bias))
            new_line.append(str(hit.bitscore))
            new_line.append(str(hit.description))
            new_line.append(str(hit.evalue))
            new_line.append(str(hit.id))
            outfile.write('\t'.join(new_line)+'\n')
    infile.close()
    outfile.close()
