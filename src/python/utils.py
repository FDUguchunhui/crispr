import numpy as np
import pandas as pd

def median_normalization(df):
    """
    Normalize the sgRNA read counts in a dataframe using the median ratio method.

    Parameters:
    df (pd.DataFrame): Dataframe containing sgRNA read counts with rows as sgRNAs and columns as experiments.

    Returns:
    pd.DataFrame: Normalized read counts with the same structure as the input dataframe.
    """

    # Calculate the geometric mean for each sgRNA across all experiments
    geometric_means = df.apply(lambda x: np.exp(np.log(x).mean()), axis=1)

    # Calculate the size factors for each experiment
    size_factors = df.divide(geometric_means, axis=0).apply(np.median, axis=0)

    # Normalize the read counts using the size factors
    normalized_counts = df.divide(size_factors, axis=1)
    normalized_counts_rounded = np.round(normalized_counts)

    return normalized_counts_rounded.astype(int)

# Example usage:
# Assuming 'data' is a pandas DataFrame with sgRNA read counts where rows are sgRNAs and columns are experiments.
# normalized_data = median_normalization(data)
