import sys
import cPickle as pickle
import gene_data

GENE_FILENAME = 'gene-data.p'
SNP_LOC_FILENAME = 'gene-to-snp-locations.p'
FREQUENCIES_FILENAME = 'population-frequencies.p'
CODON_DATA_FILENAME = 'codon-data.p'
MINOR_ALLELES_FILENAME = 'snp-minor-alleles.p'
MK_RESULTS_FILENAME_PREFIX = 'mk-results-'

# Gene data indices
POSITION_IDX = 0
BASE_IN_HUMAN_IDX = 1
BASE_IN_CHIMP_IDX = 2
SYNONYMITY_IDX = 3

# Codon data indices
CODON_STR_IDX = 0
POSITION_IN_CODON_IDX = 1

FREQUENCY_THRESHOLD = 0.01
# Allow the frequency threshold to be passed as a command line argument
try:
    print "Chosen threshold:", sys.argv[1]
    FREQUENCY_THRESHOLD = float(sys.argv[1])
except:
    print "Using default threshold of", FREQUENCY_THRESHOLD


def is_polymorphic_idx_synonymous(gene, gene_idx, minor_allele, codon_data):
    '''Returns true if replacing the polymorphic index in the gene with the
    SNP's minor allele is a synonymous polymorphism; otherwise False.
    '''
    assert minor_allele in gene_data.all_chars

    original_codon_chars = codon_data[gene][gene_idx][CODON_STR_IDX]
    polymorphic_idx = codon_data[gene][gene_idx][POSITION_IN_CODON_IDX]
    polymorphic_codon_chars = original_codon_chars[:]
    polymorphic_codon_chars[polymorphic_idx] = minor_allele

    original_codon = ''.join(original_codon_chars).upper()
    polymorphic_codon = ''.join(polymorphic_codon_chars).upper()
    is_synonymous = (
        gene_data.codon_str_to_codon[original_codon] ==
        gene_data.codon_str_to_codon[polymorphic_codon])
    return is_synonymous

# Load gene data
print "Loading gene data..."
gene_file = open(GENE_FILENAME, 'rb')
gene_to_data = pickle.load(gene_file)
gene_file.close()

# Load snp data
print "Loading SNP data..."
snp_loc_file = open(SNP_LOC_FILENAME, 'rb')
# Below is of the form gene_name -> {snp_locations}
gene_to_snp_locations = pickle.load(snp_loc_file)
snp_loc_file.close()
frequencies_file = open(FREQUENCIES_FILENAME, 'rb')
# Below is of the form gene_name -> snp_index -> (population frequencies)
population_frequencies = pickle.load(frequencies_file)
frequencies_file.close()
# The genes considered for the chosen phenotype
all_genes = gene_to_snp_locations.keys()
included_genes = [
    gene for gene in gene_to_snp_locations.keys()
    if gene_to_data[gene] is not None]

# Load codon data
print "Loading codon data..."
codon_file = open(CODON_DATA_FILENAME, 'rb')
# Below is of the form gene_name -> gene_idx -> (codon, position_in_codon), where
# gene_idx is 1-indexed, position_in_codon is 0-indexed
codon_data = pickle.load(codon_file)
codon_file.close()

# Load minor allele data
print "Loading minor allele data..."
allele_file = open(MINOR_ALLELES_FILENAME, 'rb')
# Below of the format gene_name -> gene_idx -> minor_allele, where
# gene_idx is 1-indexed, position_in_codon is 0-indexed
minor_alleles = pickle.load(allele_file)
allele_file.close()

# Generate the total number of populations
snp_idx_to_freqs = population_frequencies[population_frequencies.keys()[0]]
freq_tuple = snp_idx_to_freqs[snp_idx_to_freqs.keys()[0]]
num_populations = len(freq_tuple)

print "Generating MK values..."
print "Remaining genes:"
count = len(all_genes) + 1

# Below is of the form gene_name -> ((alpha_per_population), (ni_per_population),
#   ([P_s, P_n, D_s, D_n]_per_population))
mk_results = {}
for gene in gene_to_data:
    count -= 1
    print count

    # If this gene was not included for some reason
    if gene_to_data[gene] is None:
        print "Gene", gene, "not included."
        continue

    # Assigning initial values of 1 to provide buffer values against
    # division by 0
    P_s = [1] * num_populations
    P_n = [1] * num_populations
    D_s = [1] * num_populations
    D_n = [1] * num_populations
    for gene_location_data in gene_to_data[gene]:
        position = gene_location_data[POSITION_IDX]
        base_in_human = gene_location_data[BASE_IN_HUMAN_IDX]
        base_in_chimp = gene_location_data[BASE_IN_CHIMP_IDX]
        synonymity = gene_location_data[SYNONYMITY_IDX]

        # Check first whether the position is a polymorphism
        if position in gene_to_snp_locations[gene]:
            frequencies = population_frequencies[gene][position]
            # Iterate over each population's allele frequencies
            for population_i in xrange(num_populations):
                population_frequency = frequencies[population_i]
                # Only consider this population's location polymorphic if
                # it's above the set threshold
                if population_frequency > FREQUENCY_THRESHOLD:
                    minor_allele = minor_alleles[gene][position]
                    is_synonymous = is_polymorphic_idx_synonymous(
                        gene, position, minor_allele, codon_data)
                    if synonymity == gene_data.SYNONYMOUS:
                        P_s[population_i] += 1
                    elif synonymity == gene_data.NONSYNONYMOUS:
                        P_n[population_i] += 1

        # If the position is not a polymorphism, consider whether it is
        # a substitution
        else:
            if synonymity == gene_data.SYNONYMOUS:
                for idx in xrange(num_populations):
                    D_s[idx] += 1
            elif synonymity == gene_data.NONSYNONYMOUS:
                for idx in xrange(num_populations):
                    D_n[idx] += 1

    # Different alpha and NI for each population
    # Proportion of substitutions driven by positive selection
    alpha_per_population = [
        1 - (D_s[idx]*P_n[idx] / float(D_n[idx]*P_s[idx]))
        for idx in xrange(num_populations)]

    # Neutrality index
    neutrality_index_per_population = [
        (P_n[idx]/float(P_s[idx])) / (D_n[idx]/float(D_s[idx]))
        for idx in xrange(num_populations)]

    # Frequency parameters
    frequencies_per_population = [
        [P_s[idx], P_n[idx], D_s[idx], D_n[idx]]
        for idx in xrange(num_populations)]


    mk_results[gene] = (
        alpha_per_population,
        neutrality_index_per_population,
        frequencies_per_population)

print "Writing MK data to file..."
threshold_str = str(FREQUENCY_THRESHOLD)
threshold_str = threshold_str[threshold_str.index('.')+1:]
full_filename = MK_RESULTS_FILENAME_PREFIX + threshold_str + '.p'
pickle.dump(mk_results, open(full_filename, 'w+b'))
print "Wrote to file", full_filename + "."
