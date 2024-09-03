import pandas as pd
import os
import pyranges as pr

DISP_NAMES = {'pos':'Position',
              'var':'Variant',
              'allele':'Allele'}

def get_vcf_header(fname):
    """
    Get the comment lines (start with ##)
    and the column names from a VCF
    """
    comments = []
    with open(fname, 'r') as ifile:
        for line in ifile:
            if line.startswith('##'):
                comments.append(line)
                continue
            elif line.startswith('#'):
                header = line[1:-1].strip().split('\t')
                break
    return comments, header

def get_vcf_field_name_dict():
    """
    Return a dict mapping the abbreviation
    for each VCF field into a human-readable name
    """
    d = {'GT':'genotype',
         'AD':'allelic_depth',
         'DP':'read_depth',
         'GQ':'genotype_quality',
         'PL':'genotype_phred',
         'SB':'strand_bias'}
    return d

def add_vcf_info(df):
    """
    From a VCF DataFrame, add the metadata from the
    FORMAT column as human-readable columns
    """

    # copy so we can just merge results in later
    df_back = df.copy(deep=True)

    len1 = len(df.index)

    # split the : separated fields
    col = df.columns[-1]
    df['vcf_field'] = df['FORMAT'].str.split(':')
    df['vcf_field_val'] = df[col].str.split(':')

    # explode out, change to human readable,
    # drop unecessary columns, pivot to convert to columns
    ind_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT']
    keep_cols = ind_cols + ['vcf_field', 'vcf_field_val']
    df = df[keep_cols]

    df = df.explode(['vcf_field', 'vcf_field_val'])
    df['vcf_field'] = df['vcf_field'].map(get_vcf_field_name_dict())
    df = df.pivot(index=ind_cols,
             columns='vcf_field',
             values='vcf_field_val').reset_index()
    df.columns.name = ''
    len2 = len(df.index)

    # make sure that nothing got screwed up in the transformation
    assert len1 == len2

    # merge back in with the original thing
    df_back = df_back.merge(df, how='left',
                            on=ind_cols)
    len3 = len(df_back.index)
    assert len2 == len3

    return df_back

def read_vcf(fname):
    """
    Read VCF file which doesn't play nice w/ pandas
    """
    _, header = get_vcf_header(fname)
    df = pd.read_csv(fname, sep='\t',
                     header=None,
                     comment='#')
    df.columns = header

    # drop duplicates b/c we apparently have to do that?
    df = df.drop_duplicates()

    return df

def add_var_id(df):
    """
    Add a single-column identifier for each allele
    from long-form table output from `get_hr_vcf`
    """

    # add a single-column id
    df['var_id'] = df['CHROM']+'_'+\
                   df['POS'].astype(str)+'_'+\
                   df['REF']+'_'+\
                   df['alt_allele']

    return df

def add_pos_id(df):
    """
    Add a single-column identifier for each position
    from long or short table output from `get_hr_vcf`
    """

    # add a single-column id
    df['pos_id'] = df['CHROM']+'_'+\
                   df['POS'].astype(str)
    return df

def write_vcf(df, ofile, template):
    """
    Write VCF file w/ header from template
    """
    comments, header = get_vcf_header(template)
    with open(ofile, 'w') as out:
        for c in comments:
            out.write(c)
        header[0] = '#'+header[0]
        out.write('\t'.join(header)+'\n')
    df.to_csv(ofile,
              mode='a',
              sep='\t',
              index=False,
              header=False)

def call_var_type(x):
    """
    Call a variant as a snp, ins, or del
    """
    if x.alt_len_diff == 0:
        return 'snp'
    elif x.alt_len_diff > 0:
        return 'ins'
    elif x.alt_len_diff < 0:
        return 'del'
    else:
        return np.nan

def get_hr_vcf(fname,
               how='pos',
               filt_non_ref=True,
               filt_0_count=True,
               var_id=True,
               pos_id=True):
    """
    Get a human-readable VCF file where
    values are reported per allele (not position), and
    allele values can be expanded

    Parameters:
        fname (str): File name
        how (str): {'pos', 'var'}
            Default = 'pos'
        filt_non_ref (bool): Whether to remove positions with
            just <NON_REF> values
            Default = True
        filt_0_count (bool): Whether to remove 0-count (AD=0) alleles,
            which apparently sometime exist
            Only compatible w/ how=var.
            Default = True
        var_id (bool): Whether to add a single-column variant id.
            Only compatible w/ how=var
    """

    # read vcf and pull some info
    df = read_vcf(fname)

    # filter out non-ref stuff
    if filt_non_ref:
        df = df.loc[df['ALT']!='<NON_REF>']

    # make sure we have expected # of columns
    dataset = df.columns[-1]
    if len(df.columns) != 10:
        raise ValueError('Currently only supports single-sample VCFs')

    # add info for each vcf tag
    df = add_vcf_info(df)

    # get dataset name and add to table, drop original col
    df['dataset'] = dataset
    df.drop(dataset, axis=1, inplace=True)


    # convert to long form
    if how == 'var':
        inds = df.loc[df.allelic_depth.isnull()].index.tolist()
        if len(inds) > 0: print(f'Found {len(inds)} NaN entries. Removing.')
        df = df[~df.index.isin(inds)]

        # get a long-form table
        df['alt_allele'] = df['ALT'].str.split(',')
        df['allelic_depth'] = df['allelic_depth'].str.split(',')
        df['ref_allelic_depth'] = df.apply(lambda x: x.allelic_depth[0], axis=1)
        df['alt_allelic_depth'] = df.apply(lambda x: x.allelic_depth[1:], axis=1)
        keep_cols = ['CHROM', 'POS', 'ID', 'REF', 'QUAL',
                     'alt_allele', 'alt_allelic_depth',
                     'ref_allelic_depth', 'read_depth', 'dataset']
        df = df[keep_cols].explode(column=['alt_allele', 'alt_allelic_depth'])

        # type converstion
        num_cols = ['alt_allelic_depth', 'ref_allelic_depth', 'read_depth']
        for c in num_cols:
            df[c] = df[c].astype(int)
            
        num_cols = ['QUAL']
        for c in num_cols:
            df[c] = df[c].astype(float)

        # remove alt. alleles w/ 0 counts
        if filt_0_count:
            df = df.loc[df.alt_allelic_depth != 0]

        # compute % AF or RE for each allele
        df['alt_allele_freq'] = df['alt_allelic_depth']/df['read_depth']

        # remove non ref alt allele
        if filt_non_ref:
            df = df.loc[df.alt_allele != '<NON_REF>']

        # add single-col variant id
        if var_id == True:
            df = add_var_id(df)

        # add lengths and types of var
        df['ref_len'] = df['REF'].str.len()
        df['alt_len'] = df['alt_allele'].str.len()
        df['alt_len_diff'] = df['alt_len']-df['ref_len']
        df['var_type'] = df.apply(lambda x: call_var_type(x), axis=1)

    # add single-col pos id
    if pos_id == True:
        df = add_pos_id(df)

    return df

def write_hr_vcf(fname,
                 ofile,
                 **kwargs):
    """
    Write a VCF table in human-readable form
    """

    df = get_hr_vcf(fname, **kwargs)
    df.to_csv(ofile, sep='\t', index=False)

def vcf_to_bed(df,
               use_ref_len=False):
    """
    Convert VCF DataFrame to a BED DataFrame

    Parameters:
        df (pandas DataFrame): DataFrame of VCF data
        use_ref_len (bool): Calculate end coordinate
            based on length of reference allele
            instead of just adding 1
    """

    m = {'CHROM': 'Chromosome',
         'POS': 'Start',
         'ID': 'Name'}
    df.columns = [m[c] if c in m.keys() \
                  else c for c in df.columns]
    if use_ref_len:
        df['End'] = df['Start']+df['REF'].str.len()
    else:
        df['End'] = df['Start']+1

    return df

def write_vcf_header(ofile, template):
    """
    Write the header of a VCF file from a template
    """
    comments, header = get_vcf_header(template)
    with open(ofile, 'w') as out:
        for c in comments:
            out.write(c)
        header[0] = '#'+header[0]
        out.write('\t'.join(header)+'\n')

def bed_to_vcf(bed, ofile, template):
    """
    Convert a BED file (w/o end position) to a VCF file
    """
    write_vcf_header(ofile, template)

    # empty stuff
    if os.stat(bed).st_size == 0:
        open(ofile, 'w')
    else:
        df = pr.read_bed(bed, as_df=True)
        df.to_csv(ofile,
                  sep='\t',
                  mode='a',
                  index=False,
                  header=None)
