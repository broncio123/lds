# -*- coding: utf-8 -*-
"""
Module description ...

Reference:
- Surename1, Forename1 Initials., Surename2, Forename2 Initials, YEAR. Publication/Book title
Publisher, Number(Volume No), pp.142-161.
"""
import os
import sys
import ldds
import numpy as np
import numpy.ma as ma
import h5py
import pathlib
from scipy.integrate import solve_ivp
from functools import reduce
from scipy.interpolate import RectBivariateSpline, CubicSpline
from scipy.interpolate import interp1d, interp2d
from scipy.optimize import brentq
from ldds.hamiltonians import Hamiltonian_from_potential

def remaining_coordinate_quadratic(phase_space_axes, H0, Hamiltonian, momentum_sign):
    """
    Returns a 1D array of values for the remaining momentum coordinate, assuming the kinetic energy
    is a sum-of-squares function of the momenta (0.5*px**2+0.5*py**2+...). The sign of the  returned
    values is given by `momentum_sing`.

    Parameters
    ----------
        phase_space_axes: ndarray,
            Array of coordinate values, where the remaining momentum is 0.

        H0: float,
            Value of energy.

        Hamiltonian: function,
            Function of (t,y) that returns the value of the Hamiltonian at time t and point y.

        momentum_sign: int,
            +1 returns positive momentum values, -1 negative.

    Returns
    -------
        points_dims_remaining: ndarray,
            Array of remaining momentum values.
    """

    energy = H0-Hamiltonian(0, phase_space_axes)
    masked_energy = np.copy(energy)
    masked_energy[np.argwhere(energy<0)]=np.nan
    points_dims_remaining = momentum_sign * np.sqrt(2*masked_energy)
    return points_dims_remaining

def remaining_coordinate_value(u, ind_remaining, remaining_coordinate_bounds, H0, Hamiltonian):
    """
    Returns a 1D array of values for the remaining coordinate (not necessarily momentum) using
    Brent’s method in the bracketing interval `remaining_coordinate_bounds`.

    Parameters
    ----------
        u: ndarray, shape(n,)
            Coordinate values of a single point, where the remaining momentum .

        remaining_coordinate_bounds: ndarray, shape(2,)
            Bracketing interval for root-finding method.

        H0: float,
            Value of energy.

        Hamiltonian: function,
            Function of (t,y) that returns the value of the Hamiltonian at time t and point y.

    Returns
    -------
        point_dim_remaining: float,
            Value of coordinate or, if root-finding is unsuccessfull, np.nan.
    """

    def remaining_energy(guess):
        u[ind_remaining] = guess
        return H0 - Hamiltonian(0, u)
    try:
        point_dim_remaining = brentq(remaining_energy, remaining_coordinate_bounds[0], remaining_coordinate_bounds[1])
        return point_dim_remaining
    except:
        return np.nan

def generate_points(grid_parameters):
    """
    Returns a 1D array of all points from a on a uniform grid with dimensions and size defined by list of input parameters.
    An additional dimension initiallised with zeros is added for the calculation of Lagrangian Descriptors.
    NOTE: For n-DoF systems, currently energy conservation is only used to determine momenta dimensions.

    Parameters
    ----------
    grid_parameters : list (1-DoF systems) or dict (n-DoF systems)
        if 1-DoF, list should have two 3-tuples of floats
        entries are input parameters of limits and size of mesh per axis

        if n-DoF, dict should have the following keys
        * 'slice_parameters' : ndarray or list, should have two 3-tuples of floats, for a 2D slice
        * 'dims_slice' : ndarray or list of 0 and 1, ones indicate slice axes
        * 'dims_fixed' : ndarray or list of 0 and 1, ones indicate fixed axes
        * 'dims_fixed_values' : ndarray or list of values on fixed axes
        * 'energy_level' : float, energy value for energy conservation condition
        and either one of
        * 'Hamiltonian' : function for the Hamiltonian
        * 'potential_energy' : potential energy function, sum-of-squares kinetic energy will be assumed
        and either one of
        * 'momentum_sign' : int, -1 / 1, for negative/positive momentum for remaining axis
        * 'remaining_coordinate_bounds' : ndarray or list of values that bracket the values of the remaining coordinate (not necessarily momentum coordinate)
    Returns
    -------
    mesh : ndarray,
        Flattened array of initial conditions.
    mask : ndarray,
        Masks Nan values in further calculations.
    """
    if type(grid_parameters) == dict:
        # Unpack extra grid parameters
        slice_parameters = grid_parameters['slice_parameters']
        dims_slice = np.array(grid_parameters['dims_slice'])
        dims_fixed = np.array(grid_parameters['dims_fixed'])
        dims_fixed_values = np.array(grid_parameters['dims_fixed_values'])
        H0 = grid_parameters['energy_level']

        # If Hamiltonian is provided, use Hamiltonian, otherwise define Hamiltonian
        # from provided potential and sum-of-squares kinetic energy.
        try:
            Hamiltonian = grid_parameters['Hamiltonian']
        except:
            potential_energy = grid_parameters['potential_energy']
            Hamiltonian = Hamiltonian_from_potential(potential_energy)


        N_dim = len(dims_slice)  # Phase space dimensions

        # Check N DoF is even.
        if N_dim % 2 != 0:
            error_mssg = ("Error: Number of phase space dimensions not even. ",
                          "Check your grid parameters")
            print(error_mssg)
            sys.exit()

        # Determine number of dimensions for energy conservation.
        # There must be only one.
        dims_remaining = N_dim - np.sum(dims_fixed + dims_slice)
        if  dims_remaining > 1:
            error_mssg = ("Error: More than one remaing dimension. ",
                          "Cannot uniquelly determine values of remaining coordinates.")
            print(error_mssg)
            sys.exit()

        # Check if remaining dimension falls in configuration space
        if  np.sum(dims_slice) > 2:
            error_mssg = ("Error: 3+ dimensional slices are not available.")
            print(error_mssg)
            sys.exit()

        # Axes of visualisation slice
        ax1 = np.linspace(*slice_parameters[0])
        ax2 = np.linspace(*slice_parameters[1])

        Ax1, Ax2 = np.meshgrid(ax1, ax2)
        dims_slice_axes = np.column_stack((Ax1.ravel(),Ax2.ravel()))

        # define array phase_space_axes containing phase space points
        N_points_slice = len(dims_slice_axes)
        phase_space_axes = np.zeros((N_points_slice,N_dim))

        # fill 'slice' values of phase_space_axes
        ind_slice = np.argwhere(dims_slice==1).squeeze()
        phase_space_axes[:,ind_slice] = dims_slice_axes

        # fill 'fixed' values of phase_space_axes
        ind_fixed = np.argwhere(dims_fixed==1).squeeze()
        phase_space_axes[:,ind_fixed] = dims_fixed_values

        # determine remaining axis index
        ind_remaining = np.argwhere(dims_fixed+dims_slice==0).squeeze()
        if not (ind_remaining.shape == ()):
            error_mssg = ("Error: More than one remaing dimension. ",
                          "Cannot uniquelly determine values of remaining coordinates.")
            print(error_mssg)
            sys.exit()

        # If momentum sign is provided, determine remaining momentum values. This assumes
        # sum-of-squares kinetic energy.
        # Otherwise determine values of remaining coordinate from Hamiltonian using
        # Brent’s method in the bracketing interval `remaining_coordinate_bounds`.
        try:
            momentum_sign = grid_parameters['momentum_sign']
            phase_space_axes[:,ind_remaining] = remaining_coordinate_quadratic(
                phase_space_axes, H0, Hamiltonian, momentum_sign)
        except:
            remaining_coordinate_bounds = np.array(grid_parameters['remaining_coordinate_bounds'])
            def f_remaining_coordinate_value(u):
                return remaining_coordinate_value(u, ind_remaining, remaining_coordinate_bounds, H0, Hamiltonian)
            phase_space_axes[:,ind_remaining] = np.array(list(    \
                            map(f_remaining_coordinate_value, phase_space_axes)))

        mask = np.isnan(phase_space_axes[:,ind_remaining]) # Mask grid points
        phase_space_axes[:,ind_remaining] = np.nan_to_num(phase_space_axes[:,ind_remaining])

        # Return array of mesh points for integrator
        lagrangian_descriptor_axis = np.zeros((N_points_slice,1))
        mesh = np.column_stack((phase_space_axes,lagrangian_descriptor_axis))

        return mesh.flatten(), mask

    else:
        if len(grid_parameters) > 2:
            error_mssg = ("ERROR: grid_parameters must be a list for 2D slices for 1DoF systems. ",
                          "Provide a dictionary for a higher-dimensional systems.")
        else:
            x_min, x_max, Nx = grid_parameters[0]
            y_min, y_max, Ny = grid_parameters[1]
            points_x = np.linspace(x_min, x_max, Nx)
            points_y = np.linspace(y_min, y_max, Ny)
            X, Y = np.meshgrid(points_x, points_y)  # Grid in phase space
            # 2D grid + a zero column for LDs
            mesh = np.transpose([X.ravel(), Y.ravel(), np.zeros(Nx*Ny)])
            mask = False

            return mesh.flatten(), mask

def perturb_field(vector_field, perturbation):
    """
    Returns the vector field function with a linearly added pertubation
    Both input function should input (t, u), with t: float, and u: ndarray
    Also, the output of these funcs must be ndarrays of the same shape

    Parameters
    ----------
        vector_field: function
            unperturbed vector field

        perturbation: function
            forcing added to the vector field

    Returns
    -------
        perturbed function
    """
    return lambda t, u: vector_field(t, u) + perturbation(t, u)

def check_if_points_escape_box(u, box_boundaries):
    """
    Determine if points in phase space u have scaped box with user-defined defined dimensions

    Parameters
    ----------
    u : ndarray, shape(n, )
        points in phase space to check if outside box boundaries

    box_boundaries : list of 2-tuples of floats
        box lower and upper limits along X and Y axes

    Returns
    -------
    u_indices : ndarray, shape(n, )
        array of True/False bool values if points inside/outside the box
    """
    N_dim = u.shape[-1]
    points_positions = u.T[:int(N_dim/2)]

    if len(points_positions) == len(box_boundaries):
        check = lambda x, box_axis_limits: (box_axis_limits[0]<=x)&(x<=box_axis_limits[1])
        positions_within_box = [check(points_positions[i], box_boundaries[i]) for i in range(int(N_dim/2))]
        u_indices = reduce(lambda x, y: x&y, positions_within_box)
        return u_indices
    else:
        error_mssg = ("ERROR: Number of box axes and configuration space axes do not match. ",
                      "Check the defintion of your box boundaries.")
        sys.exit()

def lagrangian_descriptor(u, v, p_value = 0.5):
    """
    Vector field equation for Lagrangian descriptor.

    Parameters
    ----------
    v : ndarray, shape(n,2)
        Vector field at given point.

    p_value : float, optional
        Exponent in Lagrangian descriptor definition.
        0 is the acton-based LD,
        0 < p_value < 1 is the Lp quasinorm,
        1 <= p_value < 2 is the Lp norm LD,
        2 is the arclength LD.
        The default is 0.5.

    Returns
    -------
    LD : ndarray, shape(n,1)
        Vector field for Lagrangian descriptor dimension.
    """
    if p_value == 0:
        LD = np.abs(u[:,1]*v[:,0])
    elif p_value>0:
        LD = np.sum(np.abs(v)**p_value, axis=1)
    else:
        LD = np.zeros(len(u[:,0]))
    return LD

def vector_field_flat(t, points, vector_field, p_value, box_boundaries):
    """
    Returns vector field values for integration of flattened input array.

    Parameters
    ----------
    t : float
        time

    points : ndarray, shape(n,3)

    vector_field: function
        User defined vector field.

    p_value : float, optional
        Exponent in Lagrangian descriptor definition.
        0 is the acton-based LD,
        0 < p_value < 1 is the Lp quasinorm,
        1 <= p_value < 2 is the Lp norm LD,
        2 is the arclength LD.
        The default is 0.5.

    box_boundaries : list of 2-tuples, optional
        box boundaries for escape condition of variable time integration
        boundaries are infinite by default.

    Returns
    -------
    1d array
        y0 values for integrator
    """
    N_mesh_axes = 2*len(box_boundaries)+1
    u = points.reshape((-1,N_mesh_axes))
    u = u[:,:-1] #remove LD-values axis

    # Apply Escape condition
    u_inbox = check_if_points_escape_box(u, box_boundaries)

    # Define output vector field in combination with escape condition
    v = np.zeros(u.shape)
    v[u_inbox] = vector_field(t, u[u_inbox])

    # Calculate LD vector field
    LD_vec = np.zeros(len(u))
    LD_vec [u_inbox] = lagrangian_descriptor(u[u_inbox], v[u_inbox], p_value)

    # Add LD
    v_out=np.column_stack((v, LD_vec))
    return v_out.flatten()


def fit_pes(filename, clip_max = False):
    """
    Returns a 1- or 2-dimensional spline function (potential energy surface) fitted to data located in pylds/pes_files/filename.hdf5.

    Parameters
    ----------
    filename : string
        Name of file containing data:
            coords : list of ndarrays
                [x] or [x,y] contain coordinates.
            pes_data : ndarray, shape(len(x)) or shape(len(y),len(x))
                Array of potential energy values.

    clip_max : float
        Limit for clipping potential values that are not of interest.
        The default is False.

    Returns
    -------
    fspline : function
        fspline returns the potential at (x0) or (x0,y0).
    """
    dirname = "pes_files"
    dirpath = os.path.join(pathlib.Path(__file__).parent.absolute(), dirname)
    filepath = os.path.join(dirpath, filename+'.hdf5')
    
    hf = h5py.File(filepath,'r')
    coords = np.array(hf.get('coords'))
    pes_data = np.array(hf.get('pes_data'))
    hf.close()

    if clip_max:
        pes_data=np.clip(pes_data, a_min=-np.inf, a_max=clip_max)

    if len(coords) == 1:
        spline = CubicSpline(coords, pes_data)

        def fspline(positions):
            potential = np.array(list(map(spline,positions)))
            return potential

    elif len(coords) == 2:
        x, y = coords
        spline = RectBivariateSpline(x,y,pes_data.T)

        def spline_wrap(v):
            return spline(v[0],v[1]).squeeze()

        def fspline(positions):
            potential = np.array(list(map(spline_wrap,positions)))
            return potential
    else:
        print('splines in +3D are not implemented yet')
        fspline = np.nan

    return fspline

def fit_vector_field(filename):
    """
    Returns a 2-dimensional function (vector field) fitted to 
    data located in pylds/vector_field_files/filename.hdf5.

    Parameters
    ----------
    filename : string
        Name of file containing data:
        sample_time_points: 1d array,
            time-points in a sample time-interval.

        sample_coords : list of ndarrays,
            [x,y] contain coordinates.
            
        vector_field_data : ndarray, len(t),
            Array of function/vector field values.
    
    Returns
    -------
    vector_field_interpolated : function
        returns a vector field function to be evaluated at (t, u).
        with t a float and u (n,2)-array.
    """
    dirname = "vector_field_files"
    dirpath = os.path.join(pathlib.Path(__file__).parent.absolute(), dirname)
    filepath = os.path.join(dirpath, filename+'.hdf5')
    
    hf = h5py.File(filepath,'r')
    
    #extract data
    time = np.array(hf.get('sample_time'))
    coords = np.array(hf.get('sample_coords'))
    data = np.array(hf.get('vector_field_data'))
    
    hf.close()
    
    #interpolate data in time 
    v_interp_t = lambda t: interp1d(time, data.T, kind='cubic')(t).T
    
    def vector_field_wrap(v, u):
        return v(u[0],u[1])
    
    def vector_field_interpolated(t, u):
        #evaluate components of interpolated field at t
        x, y = coords
        v_t_eval_x = v_interp_t(t).T[0].reshape(len(x), len(y))
        v_t_eval_y = v_interp_t(t).T[1].reshape(len(x), len(y))
        
        #interpolate above data in space
        v_x = interp2d(x, y, v_t_eval_x, kind='cubic')
        v_y = interp2d(x, y, v_t_eval_y, kind='cubic')
        
        #evaluate components of interpolated field at u
        v_x_eval = np.array(list(map(lambda a: vector_field_wrap(v_x, a), u)))
        v_y_eval = np.array(list(map(lambda a: vector_field_wrap(v_y, a), u)))
        
        return np.column_stack([v_x_eval, v_y_eval])
    
    return vector_field_interpolated

def EulerMaruyama_solver(t_initial, u_initial, vector_field, time_step, noise_amplitude=[0, 0], noise_type="additive"):
    """
    Returns next time and state in the evolution of a stochastic ODE system via the Euler-Maruyama method.

    Euler-Maruyama scheme:

    t_next = t_initial + dt
    u_next = u_initial + v(t_initial, u_initial)*dt + b*dW

    with u = (x, y)

    NOTE: Depending on the 'noise_type' b will be a random constant (additive) or a vector (multiplicative).

    Currently only implemented for 2D vector fields.

    Parameters
    ----------
    t_initial : float
        initial time-point of all initial points in phase space.

    u_initial : array_like, shape(n,)
        initial points in phase space to determine their evolution at time t_next.

    vector_field: function
        vector field over phase space.

    time_step : float

    noise_amplitude : list of floats
        amplitude values multiplying Weiner process.

    noise_type : string
        options 'additive' (default) or 'multiplicative'.

    Returns
    -------
    t_next : float
        next time-point in the evolution of stochastic system.

    u_next : array_like, shape(n,)
        next states in the evolution of stochastic system.
    """
    #solver parameters
    dt = time_step
    v  = vector_field
    b  = np.array(noise_amplitude)

    #define Weinner process
    if noise_type == "additive":
        N_dims = u_initial.shape[1] # phase space dim
        dW = np.sqrt(abs(dt))*np.random.randn(N_dims)*np.ones(u_initial.shape)

    elif noise_type == "multiplicative":
        dW = np.sqrt(abs(dt))*np.random.randn(*u_initial.shape)

    else:
        error_mssg = ("ERROR: noise_type uknown. "
                      "Set as 'additive' or 'multiplicative'")
        print(error_mssg)

    #Euler-Maruyama iterative solver
    t_next = t_initial + dt
    u_next = u_initial + v(t_initial, u_initial)*dt + b*dW

    return t_next, u_next

def compute_lagrangian_descriptor(grid_parameters, vector_field, tau, p_value=0.5, box_boundaries=False, rtol=1.0e-4):
    """
    Returns the values of the LD function from integrated trajectories from initial conditions in phase space.

    Parameters
    ----------
    grid_parameters : list of 3-tuples of floats
        input parameters of limits and size of mesh per axis

    vector_field: function
        vector field over phase space

    tau : float
        Upper limit of integration.

    p_value : float, optional
        Exponent in Lagrangian descriptor definition.
        0 is the acton-based LD,
        0 < p_value < 1 is the Lp quasinorm,
        1 <= p_value < 2 is the Lp norm LD,
        2 is the arclength LD.
        The default is 0.5.

    box_boundaries : list of 2-tuples, optional
        Box boundaries for escape condition of variable time integration.
        Boundaries are infinite by default.

    rtol : float,
        Relative tolerance of integration step.

    Returns
    -------
    LD : ndarray, shape (Nx, Ny)
        Array of computed Lagrangian descriptor values for all initial conditions.
    """
    #get visualisation slice parameters and Number of DoF
    if type(grid_parameters) == dict:
        #n-DoF systems
        slice_parameters = np.array(grid_parameters['slice_parameters']) # 2n-D grid
        N_dim = len(grid_parameters['dims_slice'])
    else:
        #1-DoF systems
        slice_parameters = np.array(grid_parameters) # 2-D grid
        N_dim = len(slice_parameters)

    #set boundaries for escape-box condition, if not defined
    if not box_boundaries:
        box_boundaries = int(N_dim/2)*[[-np.infty, np.infty]] #restricted to configuration space

    #solve initial value problem
    f = lambda t, y: vector_field_flat(t, y, vector_field, p_value, box_boundaries)
    y0, mask = generate_points(grid_parameters)

    #mask y0 values
    if type(mask) == np.ndarray:
        mask_y0 = np.transpose([mask for i in range(N_dim+1)]).flatten()
        y0 = ma.masked_array(y0, mask=mask_y0)

    solution = solve_ivp(f, [0,tau], y0, t_eval=[tau], rtol=rtol, atol=1.0e-12)

    LD_values = solution.y[N_dim::N_dim+1] #values corresponding to LD
    LD_values[mask] = np.nan #mask LD values for slice

    N_points_slice_axes = slice_parameters[:,-1].astype('int')
    LD = np.abs(LD_values).reshape(*N_points_slice_axes) #reshape to 2-D array

    if p_value<=1:
        return LD
    else:
        return LD**(1/p_value)


__author__ = 'Broncio Aguilar-Sanjuan, Victor-Jose Garcia-Garrido, Vladimir Krajnak, Shibabrat Naik'
__status__ = 'Development'
