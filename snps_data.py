import pickle

SNPS_FILENAME = 'all.snps.in.nearby.gene.exon.distances.bed'
SNP_LOC_FILENAME = 'gene-to-snp-locations.p'
FREQUENCIES_FILENAME = 'population-frequencies.p'
MINOR_ALLELES_FILENAME = 'snp-minor-alleles.p'

GENE_IDX = 3
LOCATION_IDX = 4
MINOR_ALLELE_IDX = 6
POPULATION_IDXS = (7, 8, 9)


included_genes = set()
# Below is of the form gene_name -> {snp_locations}
gene_to_snp_locations = dict()
# Below is of the form gene_name -> snp_index -> (population frequencies)
population_frequencies = dict()
# Below is of the form gene_name -> snp_index -> minor_allele
minor_alleles = dict()

print "Opening SNPs file..."
snps_file = open(SNPS_FILENAME, 'rb')
num_snp_idxs = sum(1 for line in snps_file)
snps_file.seek(0)
print "Opened."
print "SNP locations read:", num_snp_idxs

print "Generating SNP data..."
print "Remaining SNP locations:"
count = num_snp_idxs + 1
for line in snps_file:
    count -= 1
    if count % 10000 == 0:
        print count

    # Generate snp attributes
    snp_attributes = line.split('\t')
    gene = snp_attributes[GENE_IDX]
    location = int(snp_attributes[LOCATION_IDX])
    frequencies = tuple([float(snp_attributes[idx]) for idx in POPULATION_IDXS])
    minor_allele = snp_attributes[MINOR_ALLELE_IDX]

    # Assign SNP values
    included_genes.add(gene)
    # If this is the first time this gene has been encountered
    if gene not in gene_to_snp_locations:
        gene_to_snp_locations[gene] = {location}
        population_frequencies[gene] = {location: frequencies}
        minor_alleles[gene] = {location: minor_allele}
    else:
        gene_to_snp_locations[gene].add(location)
        population_frequencies[gene][location] = frequencies
        minor_alleles[gene][location] = minor_allele

print ""
print "Generated."
num_genes = len(included_genes)
num_snps = sum([len(gene_to_snp_locations[gene]) for gene in included_genes])
avg_snps_per_gene = round(num_snps / float(num_genes), 0)
print "Number of genes considered:", len(included_genes)
print "Number of SNPs considered:", num_snps
print "Average number of SNPs per gene:", avg_snps_per_gene
print "Writing SNP information to file..."
pickle.dump(gene_to_snp_locations, open(SNP_LOC_FILENAME, 'w+b'))
pickle.dump(population_frequencies, open(FREQUENCIES_FILENAME, 'w+b'))
pickle.dump(minor_alleles, open(MINOR_ALLELES_FILENAME, 'w+b'))
print "Wrote to files:", ', '.join([
    SNP_LOC_FILENAME, FREQUENCIES_FILENAME, MINOR_ALLELES_FILENAME]) + '.'
