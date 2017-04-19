"""
Darius Douglas

Gene Set Enrichment Analysis
17 March 2017

Gene Set Enrichment Analysis
"""

from collections import defaultdict


def check_extreme_deviation(new_value, old_value, is_high_deviation):
    """
    Helper for calculate_enrichment_score(*) function
    Determine highest and lowest deviations from zero
    :param new_value: new deviation from zero
    :param old_value:
    :param is_high_deviation: boolean: check if value is a high deviation from zero or lower
    :return: if is_high deviation, return highest deviation from zero
             otherwise, return lowest deviation from zero
    """
    if is_high_deviation:
        return max(new_value, old_value)
    else:
        return min(new_value, old_value)


def calculate_enrichment_score(gene_expression_profile, gene_set, phenotype_label_data):
    """
    Walk through the gene_expression_profile
    Increase running_sum-statistic when encountering a gene in gene_set
    Decrease the running sum when encountering a gene not in gene_set
    Magnitude of the increment depends on the correlation of the gene with the phenotype
    Magnitude is obtained from phenotye_label_data
    :param gene_expression_profile: genes from a collection of samples
    :param gene_set: genes to check for in expression profile
    :param phenotype_label_data: file containing: (gene) (correlation with phenotype)*
    :return: Enrichment score as int: most extreme deviation from zero
    """
    running_sum_statistic = 0
    lowest_deviation_from_zero = 0
    highest_deviation_from_zero = 0

    for gene in gene_set:
        if gene in gene_expression_profile:
            # gene is in list, increment running sum
            print(gene)
            try:
                running_sum_statistic += phenotype_label_data[gene][0]

                # check if running_sum_statistic is highest deviation from zero
                highest_deviation_from_zero = check_extreme_deviation(
                    running_sum_statistic,
                    highest_deviation_from_zero,
                    True)

            except IndexError:
                pass  # correlation value was blank
        else:
            # gene isn't in list, decrement running sum
            try:
                running_sum_statistic -= phenotype_label_data[gene][0]

                # check if running_sum_statistic is lowest deviation from zero
                lowest_deviation_from_zero = check_extreme_deviation(
                    running_sum_statistic,
                    lowest_deviation_from_zero,
                    False)

            except IndexError:
                pass  # correlation value was blank

    # if lowest deviation is further from zero than highest deviation, return lowest
    # otherwise return highest
    return lowest_deviation_from_zero if abs(lowest_deviation_from_zero) > highest_deviation_from_zero \
        else highest_deviation_from_zero


def read_genes(genes_file):
    """
    :param words_file: file containing genes
    :return: list of genes containing in genes_file
    """
    return [gene for line in open(genes_file, 'r') for gene in line.split()]


def get_phenotype_label_data_from_file(filename):
    """
    Phenotype labels determines a gene's correlation with the phenotype
    File format: (gene) (tab) (correlationWithPhenotypeValue)
    :param filename: phennotype_label file
    :return: dictionary with gene(string) and correlationWithPhenotypeValue(int)
             i.e. {'STAT1', -13}
    """

    gene_to_phenotype_correlation_value = defaultdict(list)

    with open(filename) as f:
        for line in f:
            current_gene_info = line.split()
            # (gene, correlation_value) = line.split()
            gene_to_phenotype_correlation_value[current_gene_info[0]] = [int(item) for item in current_gene_info[1:]]

    return gene_to_phenotype_correlation_value


if __name__ == "__main__":
    """
    -Calculate enrichment score

    TODO:  Estimation of Significance Level of ES
           Adjustment for Multiple Hypothesis Testing
    """
    gene_to_enrichment_score = {}  #

    phenotype_label_data = get_phenotype_label_data_from_file('phenotype_label')

    gene_expression_profile = read_genes('gene_expression_profile.txt')

    gene_set = read_genes('gene_sets.txt')

    enrichment_score = calculate_enrichment_score(gene_expression_profile, gene_set, phenotype_label_data)
