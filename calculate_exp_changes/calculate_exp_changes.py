import os
import sys

import numpy as np
import pandas as pd

import helpers

if __name__ == '__main__':
    GENE_FILE, COUNT_DIR, SEED, UTR_FILE, OUTFILE = sys.argv[1:]

    # read in genes and transcripts
    gene_data = pd.read_csv(GENE_FILE, sep='\t')
    gene_data = gene_data.groupby('Ensembl Transcript ID').first()
    print gene_data.head()

    # get site information for every gene for the miRNA
    site_info = helpers.get_site_info(UTR_FILE, SEED)
    print len(site_info)

    print site_info.head()

    # read in all files and aggregate into one table
    all_files = []
    headers = []
    for file in os.listdir(COUNT_DIR):
        if 'exon' in file:
            filetype = '_exon'
        else:
            filetype = '_intron'
        temp = pd.read_csv(COUNT_DIR + file, delim_whitespace=True, header=None).drop_duplicates()
        header = file.split('_')[0] + filetype
        headers.append(header)
        temp.columns = [header, 'Transcript ID']
        temp['Transcript ID'] = [x.split('.')[0] for x in temp['Transcript ID']]
        temp = temp.groupby('Transcript ID').agg(sum)
        all_files.append(temp)

    count_data = pd.concat(all_files, axis=1, join='outer')
    print count_data.head()

    # add gene information and combine same genes
    print len(count_data)
    count_data = pd.concat([count_data, gene_data], axis=1, join='inner')
    print len(count_data)
    count_data = count_data.reset_index().drop(['index'], 1).groupby('Ensembl Gene ID').agg(np.nansum)
    print len(count_data)

    # print count_data.head()
    # sys.exit()

    # add site information
    count_data = pd.concat([count_data, site_info], axis=1, join='inner')
    print len(count_data)

    print count_data.head()

    # export data as csv
    count_data.to_csv(OUTFILE, sep='\t')
