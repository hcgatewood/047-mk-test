import cPickle as pickle
from pprint import pprint

CCDS_ALIGNMENT_DIRECTORY = '/broad/compbio/mwolf/50kExome/Fresco/CCDSAlign'
CCDS_GENES_FILENAME = 't2d-cds-genes.txt'
SKIPPED_FILENAME = 'gene-data-skipped-genes.txt'
GENE_DATA_FILENAME = 'gene-data.p'
CODON_DATA_FILENAME = 'codon-data.p'
ALIGNMENT_FILENAME = 'alignments.p'
CCDS_EXTENSION = '.fa'

HUMAN_TAG = '>Human'
CHIMP_TAG = '>Chimp'
TAG = '>'

SYNONYMOUS = 'synonymous'
NONSYNONYMOUS = 'nonsynonymous'
NONE = 'none'
gap_chars = ['.', '-']
all_chars = ['A', 'T', 'C', 'G'] + gap_chars
codon_str_to_codon = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}


if __name__ == '__main__':
    regenerate_alignments = False
    verbose = False

    # Parse gene names into list
    gene_names_file = open(CCDS_GENES_FILENAME, 'rb')
    gene_names = []
    for line in gene_names_file:
        for word in line.split():
            gene_names.append(word)
    # pprint(gene_names)
    # print "number of genes:", len(gene_names)

    if regenerate_alignments:
        print "Regenerating alignments dictionary..."
        print "Files left to process:"
        # Generate a dictionary of {gene_name -> [human_alignment, chimp_alignment]}
        gene_to_aligns = dict()
        count = len(gene_names) + 1
        for gene_filename in gene_names:
            count -= 1
            if count % 100 == 0:
                print count

            gene = gene_filename.strip(CCDS_EXTENSION)
            # Retrieve the text of the alignment file
            try:
                alignment_txt_file = open(
                    CCDS_ALIGNMENT_DIRECTORY + '/' + gene_filename, 'rb')
            except:
                print "File " + gene_filename + " not found."
                continue
            alignment_txt = []
            for line in alignment_txt_file:
                alignment_txt.append(line.strip('\n'))
            alignment_txt_file.close()

            # Extract the alignments for just human and chimp
            all_starts = [i for i in xrange(len(alignment_txt)) if TAG in alignment_txt[i][0]]
            human_start = alignment_txt.index(HUMAN_TAG) + 1
            human_end = min([i for i in all_starts if i > human_start])
            human_lines = alignment_txt[human_start:human_end]
            human_align = ''.join(human_lines)

            chimp_start = alignment_txt.index(CHIMP_TAG) + 1
            chimp_end = min([i for i in all_starts if i > chimp_start])
            chimp_lines = alignment_txt[chimp_start:chimp_end]
            chimp_align = ''.join(chimp_lines)

            gene_to_aligns[gene] = [human_align, chimp_align]

            if verbose:
                if len(human_align)%3 != 0 or len(chimp_align)%3 != 0:
                    print "Gene " + gene + " is not playing nice."
                print "Gene:", gene
                print "Human alignment:"
                print human_align
                print ''
                print "Chimp alignment:"
                print chimp_align
                print '\n'

        print "Generated."
        print "Writing to file..."
        pickle.dump(gene_to_aligns, open(ALIGNMENT_FILENAME, 'wb'))
        print "Wrote to", ALIGNMENT_FILENAME + "."
    else:
        print "Loading gene alignments..."
        gene_to_aligns = pickle.load(open(ALIGNMENT_FILENAME))
        print "Loaded."
        print "Number of genes:", len(gene_to_aligns)
        print ""

    # For each gene, we want to walk through and add the appropriate information
    # to our table. The format is:
    #   gene_name index_of_current_base base_in_human base_in_chimp (non)synonymous
    gene_table = []
    gene_to_entries = dict()
    # Below is of the form gene -> gene_idx -> ([codon_chars], idx_in_codon)
    codon_data = dict()
    skipped_genes = []
    print "Generating MK data table..."
    print "Genes left to process:"
    count = len(gene_names) + 1
    # Iterate over each gene we're considering
    for gene_filename in gene_names:
        count -= 1
        print count

        # Strip the extension name from the gene
        gene = gene_filename.strip(CCDS_EXTENSION)
        if verbose:
            print "Gene:", gene
        if gene not in gene_to_aligns:
            skipped_genes.append(gene)
            gene_to_entries[gene] = None
            print "Gene", gene, "not found in alignments dictionary."
            continue

        human_alignment = gene_to_aligns[gene][0]
        chimp_alignment = gene_to_aligns[gene][1]

        # Return early if invariants aren't met
        if len(human_alignment) != len(chimp_alignment):
            skipped_genes.append(gene)
            gene_to_entries[gene] = None
            print "Alignments are of different length for gene", gene
            continue
        if len(human_alignment)%3 != 0 or len(chimp_alignment)%3 != 0:
            skipped_genes.append(gene)
            gene_to_entries[gene] = None
            print "Alignment lengths are not a multiple of 3 for gene", gene
            continue

        alignment_len = len(human_alignment)
        human_base_count = 0
        human_chars = []
        chimp_chars = []
        # Iterate over each codon index in the human-chimp gene alignment
        for base_i in xrange(alignment_len):
            human_char = human_alignment[base_i]
            chimp_char = chimp_alignment[base_i]

            # If both chars are gaps, just skip this index
            if human_char in gap_chars and chimp_char in gap_chars:
                continue

            human_chars.append(human_char)
            chimp_chars.append(chimp_char)

            # Return unless we have a full codon
            if len(human_chars) < 3:
                continue

            human_str = ''.join(human_chars).upper()
            chimp_str = ''.join(chimp_chars).upper()

            # Assert that each codon is either all gaps or all chars
            human_gaps_only = all([char in gap_chars for char in human_str])
            human_chars_only = all([char not in gap_chars for char in human_str])
            chimp_gaps_only = all([char in gap_chars for char in chimp_str])
            chimp_chars_only = all([char not in gap_chars for char in chimp_str])
            violates_invariant = (
                human_gaps_only is human_chars_only or
                chimp_gaps_only is chimp_chars_only)
            if violates_invariant:
                print "Gene", gene, "violates the gap invariant."
                print "Human codon:", human_str
                print "Chimp codon:", chimp_str
                print ""
                skipped_genes.append(gene)
                gene_to_entries[gene] = None
                break

            # Return early if the human chars is all gaps
            if not human_chars_only:
                human_chars = []
                chimp_chars = []
                continue

            # Return early if something's off with the codon characters
            chimp_bad = not all([base in all_chars for base in chimp_str])
            if human_str not in codon_str_to_codon or chimp_bad:
                print "Gene", gene, "has a weird codon."
                print "Human codon:", human_str
                print "Chimp codon:", chimp_str
                print "Chimp's fault?", chimp_bad
                print ""
                skipped_genes.append(gene)
                gene_to_entries[gene] = None
                break

            # Determine whether the base synonymity
            chimp_codon = codon_str_to_codon[chimp_str] if chimp_chars_only else 'None'
            human_codon = codon_str_to_codon[human_str]
            same_codon = human_codon == chimp_codon
            synonymity = []
            for i in xrange(3):
                if human_str[i] != chimp_str[i]:
                    if chimp_str[i] in gap_chars:
                        is_synonymous = NONE
                    else:
                        is_synonymous = SYNONYMOUS if same_codon else NONSYNONYMOUS
                else:
                    is_synonymous = NONE
                synonymity.append(is_synonymous)

            # Populate the table entry vals
            # Reaching this point assures the human str is a codon
            for i in xrange(3):
                human_base_count += 1
                # Update codon data dict
                if gene not in codon_data:
                    codon_data[gene] = dict()
                codon_data[gene][human_base_count] = (human_chars[:], i)
                row = [
                    gene,
                    human_base_count,
                    human_str[i],
                    chimp_str[i],
                    synonymity[i]]
                gene_table.append(row)
                if gene in gene_to_entries:
                    gene_to_entries[gene].append(row[1:])
                else:
                    gene_to_entries[gene] = [row[1:]]
                if verbose:
                    pprint(row)
            if verbose: print ""

            # Reset the chars lists
            human_chars = []
            chimp_chars = []
        if verbose:
            print "Finished gene", gene + "."
            print "Generated value:"
            pprint(gene_to_entries[gene])

    print ""
    print "Writing gene table entries to file..."
    table_data_file = open(GENE_DATA_FILENAME, 'w+b')
    pickle.dump(gene_to_entries, table_data_file)
    table_data_file.close()
    print "Wrote to file", GENE_DATA_FILENAME + "."

    print "Writing skipped genes to file..."
    skipped_genes_file = open(SKIPPED_FILENAME, 'w+b')
    for gene in skipped_genes:
        skipped_genes_file.write(gene + '\n')
    skipped_genes_file.close()
    print "Wrote to file", SKIPPED_FILENAME + "."

    print "Writing codon data to file..."
    codon_data_file = open(CODON_DATA_FILENAME, 'w+b')
    pickle.dump(codon_data, codon_data_file)
    codon_data_file.close()
    print "Wrote to file", CODON_DATA_FILENAME + "."

    num_genes_with_data = len([
        gene for gene in gene_to_entries if gene_to_entries[gene] is not None])
    print "Generated MK data for", num_genes_with_data, "genes."
