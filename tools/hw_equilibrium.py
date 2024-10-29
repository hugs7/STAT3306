"""
Tool to calculate the test statistic and p-value for the Hardy-Weinberg equilibrium test.
Uses chi-squared with 1 degree of freedom.
"""

from typing import Dict, Tuple, Union
from utils.logger_helper import init_logger

logger = init_logger()

DEGREES_OF_FREEDOM = 1
GENOTYPES = ["AA", "AB", "BB"]
ALLELES = ["A", "B"]


def get_input() -> int:
    while True:
        try:
            value = int(input("> "))
        except ValueError:
            logger.error("Please enter a valid integer")
            return

        if value < 0:
            logger.error("Please enter a non-negative integer")
            return

        return value


def calculate_probability(observed: Dict[str, int], total: int, allele: str) -> float:
    """
    Calculate probability based on observed counts and total.

    If genotype is provided, return the probability of that genotype.
    Otherwise, return the probability of alleles A and B.

    Args:
        observed: observed counts for each genotype
        total: total number of individuals
        allele: allele to calculate probability for

    Returns:
        probability of the given genotype or alleles A and B
    """
    allele = allele.upper()
    if len(allele) == 1:
        homozygous = allele * 2
        mixed = "AB"
        return (2 * observed[homozygous] + observed[mixed]) / (2 * total)
    elif len(allele) == 2:
        return observed[allele] / total

    logger.error("Invalid allele: %s", allele)


def calcualate_expected_count(genotype: str, allele_probabilities: Dict[str, float], total: int) -> float:
    """
    Calculate the expected count for a given genotype.

    Args:
        genotype: genotype to calculate expected count for
        allele_probabilities: probabilities of alleles A and B
        total: total number of individuals

    Returns:
        expected count for the given genotype
    """
    allele1, allele2 = genotype
    if allele1 == allele2:
        # Homozygous
        return allele_probabilities[allele1] ** 2 * total
    else:
        # Heterozygous
        return 2 * allele_probabilities[allele1] * allele_probabilities[allele2] * total


def compute_chi_squared(observed: Dict[str, int], expected_counts: Dict[str, float], total: int) -> float:
    """
    Compute the chi-squared test statistic.

    Args:
        observed: observed counts for each genotype
        expected_counts: expected probabilities for each genotype
        total: total number of individuals

    Returns:
        chi-squared test statistic
    """
    chi_squared = 0
    for genotype in GENOTYPES:
        expected_count = expected_counts[genotype] * total
        chi_squared += (observed[genotype] -
                        expected_count) ** 2 / expected_count

    logger.info("Chi-squared: %.2f", chi_squared)
    return chi_squared


def main():
    logger.info("Running HW equilibrium test")

    logger.info("Observed values")
    observed_values: Dict[str, int] = {}
    for i in GENOTYPES:
        logger.info(f"Enter observed count with genotype {i}: ")
        observed = get_input()
        observed_values[i] = observed

    obs_total = sum(observed_values.values())
    logger.info("Total observed count: %d", obs_total)

    logger.info("Computing expected values")
    allele_probabilities = {allele: calculate_probability(
        observed_values, obs_total, allele) for allele in ALLELES}

    genotype_probabilities = {genotype: calculate_probability(
        observed_values, obs_total, genotype) for genotype in GENOTYPES}

    expected_counts = {genotype: calcualate_expected_count(
        genotype, allele_probabilities, obs_total) for genotype in GENOTYPES}
    chi_squared = compute_chi_squared(
        observed_values, expected_counts, obs_total)


if __name__ == "__main__":
    main()
