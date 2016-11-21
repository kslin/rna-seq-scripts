from string import maketrans
import sys

import numpy as np
import pandas as pd

def rev_comp(seq):
    """
    Parameters:
    ==========
    seq: string, sequence to get reverse complement of

    Returns:
    =======
    float: reverse complement of seq in all caps
    """
    seq = seq.upper().replace('T','U')
    intab = "AUCG"
    outtab = "UAGC"
    trantab = maketrans(intab, outtab)
    seq = seq[::-1]
    seq = seq.translate(trantab)
    return seq

def get_site_info(utr_file, seed):
    """
    Given a utr file and seed, find all genes with a seed match in the 3'utr
    """
    SITE = rev_comp(seed)
    SIXMER = SITE[1:]

    UTRS = pd.read_csv(utr_file, sep='\t')#.set_index('Ensembl ID')
    UTRS['UTR sequence'] = [x.replace('-','').upper().replace('T','U') for x in UTRS['UTR sequence']]
    UTRS['num 7mer-m8'] = [x.count(SITE) for x in UTRS['UTR sequence']]
    UTRS['num sites'] = [x.count(SIXMER) for x in UTRS['UTR sequence']]
    UTRS['Gene ID'] = [x.split('.')[0] for x in UTRS['Gene ID']]
    UTRS = UTRS[['Gene ID', 'num sites', 'num 7mer-m8']]
    print len(UTRS)
    UTRS = UTRS.groupby('Gene ID').first()
    print len(UTRS)
    # .set_index('Gene ID')
    # UTRS.to_csv('../../data/mir_data/{}_data.txt'.format(mir_name),sep='\t')

    return UTRS


def get_cdf(logFC_list):
    """
    Parameters:
    ==========
    logFC_list: list of floats, list of log fold-change values to plot CDF
    
    Returns:
    =======
    list of floats: bin indices
    list of floats: cdf values corresponding to the bin indices
    """
    if len(logFC_list) < 5:
        return [],[]
    num_bins = len(logFC_list)/5
    counts,bin_edges = np.histogram(logFC_list,bins=num_bins)
    counts = counts / float(sum(counts))
    return bin_edges[1:],np.cumsum(counts)


def naive_logFC(exon_file, intron_file):
    """
    columns for both files: mir_rep1, mir_rep2... mock_rep1, mock_rep2...
    indexed by transcript name

    1) Normalize for read count.
    2) Remove genes whose expression is less than the median expression level.
    3) Compute average RPM for reps
    4) Calculate log fold-change.
    """

    # read in data
    exons = pd.read_csv(exon_file, sep='\t').dropna()
    introns = pd.read_csv(intron_file, sep='\t').dropna()
    mir_cols = []
    mock_cols = []
    for col in exons.columns:
        if col.split('_')[0] == 'mir':
            mir_cols.append(col)
        else:
            mock_cols.append(col)

        # normalize for read count
        exons_normed = np.array(exons[col].astype(float)) / float(np.nanmean(exons[col]))
        med_exon = np.median(exons_normed)

        introns_normed = introns[col].astype(float) / float(np.nanmean(introns[col]))
        med_intron = np.median(introns_normed)

        # filter for low expression
        exons[col] = [x if x > med_exon else None for x in exons_normed]
        introns[col] = [x if x > med_exon else None for x in introns_normed]

    exons = exons.dropna()
    introns = introns.dropna()

    # calculate average RPM
    exons['average_mir'] = np.nanmean(exons[mir_cols], axis=1)
    exons['average_mock'] = np.nanmean(exons[mock_cols], axis=1)
    introns['average_mir'] = np.nanmean(introns[mir_cols], axis=1)
    introns['average_mock'] = np.nanmean(introns[mock_cols], axis=1)

    # calculate logFC
    exons['exon logFC'] = np.log2(exons['average_mir'] / exons['average_mock'])
    introns['intron logFC'] = np.log2(introns['average_mir'] / introns['average_mock'])

    # merge exons and introns
    merged = pd.concat([exons[['exon logFC']], introns[['intron logFC']]], axis=1, join='outer')

    return merged

def edgeR(exon_file, intron_file, exon_file_mock, intron_file_mock):
    pass

def EISA(exon_file, intron_file, exon_file_mock, intron_file_mock):
    pass

def new_EISA(exon_file, intron_file, exon_file_mock, intron_file_mock):
    pass