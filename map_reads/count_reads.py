import sys

def resolve(label):
    """
    Given a result from htseq-count, resolve label
    """
    if len(label) == 0:
        return 'ID'

    # one label -> return label
    elif '__' not in label:
        pieces = label.split(':')
        return '{}:{}'.format(pieces[2], pieces[-1][:-1])

    # no feature
    elif 'no_feature' in label:
        return 'no_feature'

    # not aligned
    elif 'not_aligned' in label:
        return 'not_aligned'

    # ambiguous mapping
    elif 'ambiguous' in label:
        if ('exon' in label):

            # map to both exons and introns
            if ('intron' in label):
                ids = label.split('[')[1][:-2].split('+')
                gene_ids = list(set([x.split(':')[-1] for x in ids]))

                # if it maps to one exon and one intron of the same transcript, call intron/exon junction
                if len(ids) == 2:
                    transcript_ids = [x.split(':')[1] for x in ids]
                    if (transcript_ids[0] == transcript_ids[1]):
                        return 'intronexonjunction:{}'.format(gene_ids[0])
                
                # if it maps to exons and introns of different transcripts, same gene, call gene + ambiguous
                if len(gene_ids) == 1:
                    return 'ambiguous:{}'.format(gene_ids[0])

                # otherwise, just call ambiguous
                return 'ambiguous_intron_exon'

            # if it maps to exons of the same gene, call gene, otherwise call ambiguous
            ids = label.split('[')[1][:-2].split('+')
            ids = list(set([x.split(':')[-1] for x in ids]))
            if len(ids) != 1:
                return 'ambiguous_mult_genes'
            else:
                return 'exon:{}'.format(ids[0])

        # if it maps to introns of the same gene, call gene, otherwise call ambiguous
        elif ('intron' in label):
            ids = label.split('[')[1][:-2].split('+')
            ids = list(set([x.split(':')[-1] for x in ids]))
            if len(ids) != 1:
                return 'ambiguous_mult_genes'
            else:
                return 'intron:{}'.format(ids[0])
    else:
        return 'other'


if __name__ == '__main__':
    INFILE, OUTFILE = sys.argv[1:]
    current = ''
    count = 'count'
    count_dict = {}
    i = 0

    # read annotations, compile counts and resolve labels
    with open(INFILE, 'r') as infile:
        for line in infile:
            if line != current:
                resolved = resolve(current)
                if resolved in count_dict:
                    count_dict[resolved] += count
                else:
                    count_dict[resolved] = count
                current = line
                count = 1
            else:
                count += 1

        if current in count_dict:
            count_dict[current] += count
        else:
            count_dict[current] = count

    # write compiled counts to outfile
    with open(OUTFILE, 'w') as outfile:
        for key, count in count_dict.items():
            outfile.write('{}\t{}\n'.format(key, count))
