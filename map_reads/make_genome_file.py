import sys

import pandas as pd

if __name__ == '__main__':
	bed_file = sys.argv[1]

	coords = pd.read_csv(bed_file, sep='\t', header=None)
	coords['order'] = range(len(coords))
	coords = coords.groupby(0).agg({1:min, 2:max, 'order': min}).sort_values(by='order')

	coords['length'] = coords[2] - coords[1]

	print coords.drop('order', 1)

	coords.to_csv(outfile, sep='\t',)