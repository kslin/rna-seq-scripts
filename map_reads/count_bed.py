import sys

import pandas as pd

if __name__ == '__main__':
    infile, outfile = sys.argv[1:]
    counts = {}
    with open(infile, 'rb') as f:
        for line in f:
            line = line.split('\t')
            transcript_name = line[-3]
            if transcript_name in counts:
                counts[transcript_name] += 1
            else:
                counts[transcript_name] = 1

    data = pd.DataFrame({'Transcript_ID': counts.keys(),
                         'count': counts.values()})

    data = data.sort_values(by='Transcript_ID')

    data.to_csv(outfile, sep='\t', index=False)