# Test suite for obtaining unstable periodic orbits in the uncoupled 
# quartic Hamiltonian

import numpy as np

import h5py

import unittest

import ldds
from ldds.base import compute_lagrangian_descriptor
from ldds.vector_fields import quadratic_normalform_saddlecenter
from ldds.hamiltonians import quadratic_normal_form_saddlecenter_ham


import os
# path_to_data = os.path.join(os.path.dirname(os.path.dirname(__file__)), \
                            # 'data/')

from scipy.spatial.distance import directed_hausdorff


class TestContourMap(unittest.TestCase):
    """ 
    Test Lagrangian descriptor contour map for uncoupled quadratic normal form Hamiltonian with an index one saddle.

    1. Uses numpy assert to compare contour map values
    2. Uses Hausdorff metric to verify 
    """

    def test_uncoupled_ham2dof(self):
        """ 
        Compare with provided data for the system on (q1,p1) and (q1,q2) pair of coordinates 
        """
    
        Hamiltonian = quadratic_normal_form_saddlecenter_ham
        vector_field = quadratic_normalform_saddlecenter

        H0 = 0.1 # energy of the initial conditions

        # Integration parameters
        tau = 10

        # LDp, p-value
        p_value = 0.5

        # Mesh parameters
        x1_min, x1_max = [-1.0, 1.0]
        x2_min, x2_max = [-1.0, 1.0]
        x1_res, x2_res = [101, 101]
        slice_parameters = [[x1_min, x1_max, x1_res],[x2_min, x2_max, x2_res]]

        dims_fixed = [0,1,0,0] # Variable ordering (x1 x2 y1 y2)
        dims_fixed_values = [0] # This can also be an array of values
        dims_slice = [1,0,1,0] # Visualisation slice
        momentum_sign = 1 # Direction of momentum that defines the slice - (1) positive / (-1) negative
        grid_parameters = {
            'slice_parameters' : slice_parameters,
            'dims_slice' : dims_slice,
            'dims_fixed' : dims_fixed,
            'dims_fixed_values' : dims_fixed_values,
            'momentum_sign' : momentum_sign,
            'Hamiltonian': Hamiltonian,
            'energy_level': H0
        }
        

        # Obtain LD from the package
        forward_ld = compute_lagrangian_descriptor(grid_parameters, quadratic_normalform_saddlecenter, tau, rtol=1.0e-8)
        backward_ld = compute_lagrangian_descriptor(grid_parameters, quadratic_normalform_saddlecenter, -tau, rtol=1.0e-8)

        # Load benchmark data
        hf_data = h5py.File('./benchmark_data/quadratic_ham2dof/test_M100x100_finalT10_E1e-1_q1p1.h5', 'r')
        ld_benchmark = np.array(hf_data.get('LD_q1p1_q2zero'))
        hf_data.close()


        np.testing.assert_array_almost_equal(forward_ld+backward_ld, 
                                            ld_benchmark, decimal = 3)


        hausd_dist = directed_hausdorff(forward_ld+backward_ld, 
                                        ld_benchmark)[0]
        print(hausd_dist)

        # dims_fixed = [0,0,1,0] # Variable ordering (x1 x2 y1 y2)
        # dims_fixed_values = [0] # This can also be an array of values
        # dims_slice = [1,1,0,0] # Visualisation slice
        # momentum_sign = 1 # Direction of momentum that defines the slice - (1) positive / (-1) negative
        # grid_parameters = {
        #     'slice_parameters' : slice_parameters,
        #     'dims_slice' : dims_slice,
        #     'dims_fixed' : dims_fixed,
        #     'dims_fixed_values' : dims_fixed_values,
        #     'momentum_sign' : momentum_sign,
        #     'Hamiltonian': Hamiltonian,
        #     'energy_level': H0
        # }
        

        # # Obtain LD from the package
        # forward_ld = compute_lagrangian_descriptor(grid_parameters, quadratic_normalform_saddlecenter, tau, rtol=1.0e-8)
        # backward_ld = compute_lagrangian_descriptor(grid_parameters, quadratic_normalform_saddlecenter, -tau, rtol=1.0e-8)

        # # Load benchmark data
        # hf_data = h5py.File('./benchmark_data/quadratic_ham2dof/test_M100x100_finalT10_E1e-1_q1q2.h5', 'r')
        # ld_benchmark = np.array(hf_data.get('LD_q1q2_p1zero'))
        # hf_data.close()

        # np.testing.assert_array_almost_equal(forward_ld+backward_ld, 
                                            # ld_benchmark, decimal = 3)



    # def test_hausdorff_distance(self):
    #     """Obtain Hausdorff distance between the benchmark LD data and computed LD"""

    #     total_energy = deltaE_val + parameters[2]
    #     # y, py = uncoupled.upo_analytical(total_energy, t, parameters)

    #     numerical_orbit = np.array([orbit[:,1], orbit[:,3]])
    #     analytical_orbit = np.array([y, py])

    #     hausd_dist = directed_hausdorff(numerical_orbit, analytical_orbit)[0]

        


if __name__ == "__main__":
    unittest.main()








