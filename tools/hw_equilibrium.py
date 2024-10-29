"""
Tool to calculate the test statistic and p-value for the Hardy-Weinberg equilibrium test.
Uses chi-squared with 1 degree of freedom.
"""

from utils.logger_helper import init_logger

logger = init_logger()

DEGREES_OF_FREEDOM = 1


def get_input():
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


def main():
    logger.info("Running HW equilibrium test")

    logger.info("Observed values")
    observed_values = {}
    for i in ["AA", "AB", "BB"]:
        logger.info(f"Enter observed count with genotype {i}: ")
        observed = get_input()
        observed_values[i] = observed

    obs_total = sum(observed_values.values())
    logger.info("Total observed count: %d", obs_total)

    logger.info("Computing expected values")
    pr_a = (2 * observed_values["AA"] +
            observed_values["AB"]) / (2 * obs_total)
    pr_b = (2 * observed_values["BB"] +
            observed_values["AB"]) / (2 * obs_total)

    pr_aa = observed_values["AA"] / obs_total
    pr_ab = observed_values["AB"] / obs_total
    pr_bb = observed_values["BB"] / obs_total


if __name__ == "__main__":
    main()
