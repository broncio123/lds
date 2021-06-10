import numpy as np
from operator import itemgetter
from ldds.base import generate_points, lagrangian_descriptor

def check_if_points_escape_box(u, box_boundaries):
    """
    Determine if points u in 2D plane have escaped from box_boundaries limits.
    
    Parameters
    ----------
    u : ndarray, shape(n, 2),
        Points in plane.
    
    box_boundaries : list of 2 tuples of floats,
        Values are interpreted as [[x_min,x_max], [y_min, y_max]].
        
    Returns
    -------
    u_indices : ndarray of bools, shape(n, 2),
        True/False for points inside/outside the box_boundaries respectively.
    """
    x, y = u.T
    # Escape condition
    box_x_min, box_x_max = box_boundaries[0]
    box_y_min, box_y_max = box_boundaries[1]
    u_indices = (x >= box_x_min) & (x <= box_x_max) & (y >= box_y_min) & (y <= box_y_max)
    return u_indices

def pbc_correction_coords_single_axis(x, box_origin, box_length):
    """
    Correct single coordinate on a periodic domain.

    Parameters
    ----------
    x : ndarray, shape(n,)
        Coordinate values.
    
    box_origin : float or False(bool),
        Values of perdiodic domain origin. If False, no correction is applied.
        
    box_length : float or False(bool)
        Length of periodic doamin. If False, no correction is applied
        
    Returns
    -------
    x_pbc : ndarray, shape(n,)
        Corrected coordinate values.
    """
    x_pbc = x
    x0 = box_origin
    L = box_length
    if not isinstance(x0, bool) or not isinstance(L, bool):
        #apply PBC correction
        x = x + L/2 - x0
        x = np.mod(x + 2*L, L) 
        x_pbc = x - L/2 + x0
        
    return x_pbc

def pbc_correction_coords(u, periodic_boundaries):
    """
    Correct coordinates on a periodic domain.
    The periodic box can be a rectangle (torus) or an infinite band (cylinder) in the plane.
    
    Parameters
    ----------
    u : ndarray, shape(n, 2)
        Points in 2D plane to be correct.
    
    periodic_boundaries : list of 2-tuples of floats or False (bools)
        The first tuple must be the box origin coordinates in the plane.
        The second tuple must be the lengths of the box edges.
        For a torus PBC, both tuples entries must be floats.
        For a cylinder PBC, tuple entries for the non-periodic axis must be False.
        If no PBCs, this must be a False (single bool).
        
    Returns
    -------
    u_pbc : ndarray, shape(n, 2)
        Array of corrected coordinate values inside periodic domain.
    """
    N_dims = len(periodic_boundaries)
    pbc = periodic_boundaries
    f = pbc_correction_coords_single_axis
    u_pbc = [f(u.T[i], *pbc[i]) for i in range(N_dims)]
    
    return np.column_stack(u_pbc)

def pbc_correction_distances_single_axis(dx, box_length):
    """
    Correct the distances between a point and its image along a single coordinate on a periodic domain.

    Parameters
    ----------
    dx : ndarray, shape(n,)
        Distances between a point and its image along a single coordinate.
        
    box_length : float or False(bool)
        Length of periodic doamin. If False, no correction is applied.
        
    Returns
    -------
    dx_pbc : ndarray, shape(n,)
        Array of corrected distances on a periodic domain.
    """
    dx_pbc = dx
    L = box_length
    if not isinstance(L, bool):
        L = abs(L)
        nint = lambda x: np.round(x).astype(int) #nearest integer
        dx_pbc = dx - nint(dx/L) #minimum image criterion
    
    return dx_pbc
    
def pbc_correction_distances(du, periodic_boundaries):
    """
    Correct the distances between a point and its image in a periodic domain in 2D.
    The periodic domain can be a rectangle (torus) or an infinite band (cylinder) in the plane.

    Parameters
    ----------
    du : ndarray, shape(n, 2)
        Distances between a point and its image to be corrected.
    
    periodic_boundaries : list of 2-tuples of floats or False (bools)
        The first tuple must be the box origin coordinates in the plane.
        The second tuple must be the lengths of the box edges.
        For a torus PBC, both tuples entries must be floats.
        For a cylinder PBC, tuple entries for the non-periodic axis must be False.
        If no PBCs, this must be a False (single bool).
        
    Returns
    -------
    du_pbc : ndarray, shape(n, 2)
        Array of corrected distances in a periodic domain.
    """
    N_dims = len(periodic_boundaries)
    pbc = periodic_boundaries
    f = pbc_correction_distances_single_axis
    du_pbc = [f(du.T[i], pbc[i][-1]) for i in range(N_dims)]
    
    return np.column_stack(du_pbc)

def compute_lagrangian_descriptor(grid_parameters, discrete_map, N_iterations, p_value=0.5, box_boundaries=False, periodic_boundaries=False):
    """
    Returns Lagrangian descriptor values calculated along trajectories of a discrete system in plane.
    
    Parameters
    ----------
    grid_parameters : list of 3-tuples of floats,
        Limits and density from which a uniform grid of initial conditions will be generated.
    
    discrete_map: function,
        2D map defining a discrete dynamical system.
        
    N_iterations : int,
        Number of iterations (discrete-time equivalent of integration time).
        
    p_value : float, optional,
        Exponent in Lagrangian descriptor definition.
        0 is the acton-based LD,
        0 < p_value < 1 is the Lp quasinorm,
        1 <= p_value < 2 is the Lp norm LD,
        2 is the arclength LD.
        The default is 0.5.
    
    box_boundaries : list of 2-tuples, optional,
        Limits to prevent system from reaching infinite values (finite time blow-up).
        Trajectories that escape box_boundaries will not be further iterated.
        Default is False.
        
    perodic_boundaries: list of floats, optional,
        Values for periodic domain.
        Default is False.
    
    Returns
    -------
    LD : ndarray, shape (Nx, Ny),
        Array of Lagrangian descriptor values for all initial conditions on a regular grid.
    """
    t0 = 0 #initial iteration
    
    N_mesh_axes = len(grid_parameters)+1
    y0, mask = generate_points(grid_parameters)
    y0 = y0.reshape(-1,N_mesh_axes)
    y0 = y0[:,:-1] # exclude LD-axis
    
    f = discrete_map

    LD_values = np.zeros(len(y0))
    for i in range(N_iterations):
        t, y = f(t0, y0)
        # Escape box condition
        if box_boundaries:
            y_inbox = check_if_points_escape_box(y, box_boundaries)
            y[y_inbox == False] = y0[y_inbox == False]
        
        dt = t-t0
        dy = y-y0
        
        # Periodic Boundary conditions
        if periodic_boundaries:
            dy = pbc_correction_distances(dy, periodic_boundaries)
            y0 = pbc_correction_coords(y0, periodic_boundaries)
            y  = pbc_correction_coords(y , periodic_boundaries)
                
        LD_values = LD_values + lagrangian_descriptor(y0, dy, p_value)
        
        t0 = t0+dt
        y0 = y

    N_points_slice_axes = [x[-1] for x in grid_parameters] #take number of points
    LD = LD_values.reshape(*N_points_slice_axes) #reshape to 2-D array  

    if p_value<=1:
        return LD
    else:
        return LD**(1/p_value)
