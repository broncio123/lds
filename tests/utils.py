# Utilites functions

import numpy as np

def find_max_min_nan_mat(input_mat):
    """Find maximum and minimum values in a matrix with NaN values"""
    idx_row, idx_col = np.where(~np.isnan(MMesh))

    MMesh_stripped = MMesh[idx_row, idx_col]

    max_val = np.max(MMesh_stripped)
    min_val = np.min(MMesh_stripped)
                
    return min_val, max_val

