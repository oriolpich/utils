import pandas as pd
import fire
import matplotlib.pyplot as plt
import gzip
import io
import os
from pybedtools import BedTool
import matplotlib as mpl
from collections import defaultdict

def config_params(font_size=7):

    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rcParams['font.sans-serif'] = ['arial']
    plt.rcParams['font.size'] = font_size
    plt.rcParams['font.family'] = ['sans-serif']
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.cal'] = 'arial'
    plt.rcParams['mathtext.rm'] = 'arial'
    

def read_vcf(path):
    """
    Reads a vcf. From Stackoverflow.

    Args:
        path ([str]): path to vcf

    Returns:
        [df]: dataframe with correct column names
    """
    if '.gz' in os.path.basename(path):
        with gzip.open(path, 'rt') as f:
            lines = [l for l in f if not l.startswith('##')]
    else:
        with open(path, 'rt') as f:
            lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})


def chunks(l, n):
    """
    Split into even chunks

    Args:
        l (list): 
        n (int): number of chunks

    Yields:
        [list]: 
    """
    for i in range(0, len(l), n):
        yield l[i:i + n]
    
    
def snvs_order():
    """
    Get correct order for plotting

    Returns:
        list: 
    """
    
    order = []
    first = ['A', 'C', 'G', 'T']
    pyr = ['C', 'T']
    for p in pyr:
        for mut in first:
            if mut != p:
                for f in first:
                    for f2 in first:
                        comb = '{}[{}>{}]{}'.format(f, p, mut, f2)
                        order.append(comb)
    return order

def create_snv_class(df):
    """
    Create the channels, in pyrimidine reference
    Args:
        df (df): 

    Returns:
        str: Fixed channels
    """
    
    pyr = ['C', 'T']
    rev = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N':'N'}
    x = df['TRIPLET']
    if x[1] in pyr:
        out = '{}[{}>{}]{}'.format(x[0], x[1], df['ALT'], x[2])
    else:
        out = '{}[{}>{}]{}'.format(rev[x[2]], rev[x[1]], rev[df['ALT']], rev[x[0]])
    
    return out

def plot_snvs(sig, title, outpath='', fontsize = 12):
    """
    Plotting function

    Args:
        sig (list): vector, size 96 
        title (str): [description]
        outpath (str, optional): Defaults to ''.
        fontsize (int, optional): Defaults to 12.
    """

    config_params(fontsize)

    fig, axs = plt.subplots(
        nrows=2, ncols=1, figsize=(10.2, 3), gridspec_kw={'height_ratios': [1, 9]}
    )
    
    plt.title(title)
    order_plot = snvs_order()

    vals = []
    colors = []
    colors_mut = [
        '#1ebff0', '#050708', '#e62725', '#cbcacb', '#a1cf64', '#edc8c5'
    ]
    bot = -0.5
    for ix, c in enumerate(chunks(sig, 16)):
        colors.extend([colors_mut[ix] for s in c])
        axs[0].barh(1, 16, left=bot, color=colors_mut[ix])
        bot += 16
        vals.extend(c)

    axs[0].set_xlim(-1, 96)
    axs[0].spines['top'].set_visible(False)
    axs[0].spines['bottom'].set_visible(False)
    axs[0].spines['left'].set_visible(False)
    axs[0].spines['right'].set_visible(False)

    axs[0].get_yaxis().set_visible(False)
    axs[0].get_xaxis().set_visible(False)

    x = [i for i in range(len(vals))]

    axs[1].axhline(y=0.05, xmin=-1, xmax=96, lw=0.6, color='grey', alpha=0.2)
    axs[1].axhline(y=0.1, xmin=-1, xmax=96, lw=0.6, color='grey', alpha=0.2)
    axs[1].axhline(y=0.15, xmin=-1, xmax=96, lw=0.6, color='grey', alpha=0.2)

    axs[1].bar(x, vals, color=colors, width=0.8, linewidth=0, align='center')
    axs[1].set_xticks(x)
    axs[1].set_xticklabels(
        ['{}{}{}'.format(a[0], a[2], a[-1]) for a in order_plot],
        verticalalignment="center", ha='center', rotation=90, fontsize=6,
        color='grey'
    )

    plt.tight_layout()
    plt.xlim(-1, 96)

    axs[1].spines['top'].set_visible(False)
    axs[1].set_ylabel('Mutational Count')
    axs[1].spines['right'].set_visible(False)

    plt.setp([axs[1].get_xticklines(), axs[1].get_yticklines()], color='grey')

    axs[1].xaxis.set_ticks_position('none')
    for axis in ['top', 'bottom', 'left', 'right']:
        axs[1].spines[axis].set_linewidth(0.2)

    axs[1].xaxis.set_tick_params(pad=0.5)
    axs[1].yaxis.set_tick_params(pad=0.5, width=0.5)

    plt.tick_params(axis='both', which='both', bottom=False, left=False)
    
    plt.savefig(outpath + title + '.svg')
    plt.savefig(outpath + title + '.png', dpi = 600)

    plt.close()
    
def get_trinucleotide_context(df, genome):
    """
    Using bedtools to get the trinucleotide context

    Args:
        df (df): 
        genome (str): 

    Returns:
        [type]: [description]
    """
    
    # position minus two because it is a bedfile
    df['pos-2'] = df['POS'] - 2
    df['pos+2'] = df['POS'] + 1

    subs = BedTool.from_dataframe(df[['CHROM', 'pos-2', 'pos+2']])

    # this is really slow, but we want to avoid bgreference
    subs = subs.sequence(fi=genome, name = False)
    trinucleotide = []
    with open(subs.seqfn, 'rt') as infile:
        for count, line in enumerate(infile, start=1):

            # only read every two lines because the header of the fasta is also provided
            if count % 2 == 0:
                trinucleotide.append(line.rstrip().upper())
    
    return trinucleotide
           
def plot_profile(data, sample):
    """
    Get the counts of each variant class, create the vector and plot it
    Args:
        data (df): 
        sample (str): 
    """
    
    dic_counts = data['VARIANT_CLASS'].value_counts().to_dict()
    signature = []
    
    for s in snvs_order():
        signature.append(dic_counts.get(s, 0))

    plot_snvs(signature, sample)

    
def mutational_profile(f, genome):
    """
    This function will read a tsv file and plot the mutational spectra
    Args:
        f (str): A path to a tsv with CHROM, POS, REF, and ALT. A VCF is also supported. SAMPLE column is optional, if provided
                it will plot at the sample level.
        genome (str): Path to the fasta sequence to get the trinucleotide context.
    """
    
    # read the file
    if 'vcf' in os.path.basename(f):
        df = read_vcf(f)
    else:
        df = pd.read_csv(f, sep ='\t')
           
    # select whether we have SNVs or others
    df['len_alt'] = df['ALT'].str.len()

    # number of characters in ref
    df['len_ref'] = df['REF'].str.len()

    # first classification between SNV and others
    df['TYPE'] = df.apply(
        lambda x: 'SNV' if ((x['len_alt'] == 1) and (x['len_ref'] == 1)
                            and (x['ALT'] != '-') and (x['REF'] != '-')) else 'INDEL', axis=1
    )
    
    # select SNVs and get the trincleotide context and the variant type
    df = df[df['TYPE']=='SNV']
    trinucleotide = get_trinucleotide_context(df, genome)
    df['TRIPLET'] = trinucleotide
    df['VARIANT_CLASS'] = df.apply(create_snv_class, axis=1)
    
    if 'SAMPLE' in df.columns:
        for sample, data in df.groupby(by='SAMPLE'):
            plot_profile(data, sample)
    
    else:
        plot_profile(df, os.path.basename(f))
    
    dic_matrix = defaultdict(dict)
    if 'SAMPLE' in df.columns:
        for sample, data in df.groupby(by='SAMPLE'):
            dic_s = data['VARIANT_CLASS'].value_counts().to_dict()
            for s in snvs_order():
                dic_matrix[sample][s] = dic_s.get(s, 0)
    else:
        sample = 'SAMPLE'
        dic_s = df['VARIANT_CLASS'].value_counts().to_dict()
        for s in snvs_order():
            dic_matrix[sample][s] = dic_s.get(s, 0)

    matrix = pd.DataFrame(dic_matrix).fillna(0).loc[snvs_order()].astype(int)
    matrix.to_csv(f + '.matrix.tsv', sep ='\t', index = True, header = True)
if __name__ == '__main__':
    fire.Fire(mutational_profile)

