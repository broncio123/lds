import numpy as np

def HamCenter1D_Hamiltonian(t, u, PARAMETERS = [1]):
    """
    Returns Hamiltonian for a 1DoF centre, for an array of points in phase space.
    Number of model parameters: 1 . PARAMETERS = [omega].
    Functional form: H = omega/2*(y**2 + x**2), with u = (x, y).

    Parameters
    ----------
    t : float,
        Time. (This Hamiltonian is independent of time.)

    u : ndarray, shape(n,),
        Points in phase space.

    PARAMETERS : list of floats,
        vector field parameters.

    Returns
    -------
    H : ndarray, shape(n,),
        Hamiltonian at points u, in phase space at time t.
    """
    x, y = u.T
    # Hamiltonian Model Parameter
    omega, = PARAMETERS
    return 0.5*omega*(y*y + x*x)

def HamSaddle1D_Hamiltonian(t, u, PARAMETERS = [1]):
    """
    Returns Hamiltonian for 1DoF saddle, for an array of points in phase space.
    Number of model parameters: 1 . PARAMETERS = [lamda].
    Functional form: H = omega/2*(y**2 - x**2), with u = (x, y).

    Parameters
    ----------
    t : float,
        Time. (This Hamiltonian is independent of time.)

    u : ndarray, shape(n,),
        Points in phase space.

    PARAMETERS : list of floats,
        vector field parameters.

    Returns
    -------
    H : ndarray, shape(n,),
        Hamiltonian at points u, in phase space at time t.
    """
    x, y = u.T
    # Hamiltonian Model Parameter
    lamda, = PARAMETERS
    return 0.5*lamda*(y*y - x*x)

def Duffing1D_Hamiltonian(t, u, PARAMETERS = [1, 1]):
    """
    Returns Hamiltonian for the Duffing oscillator. PARAMETERS = [alpha, beta]
    Functional form: H = 1/2*p_x**2 - alpha*x**2 + beta*x**4), with u = (x, p_x).

    Parameters
    ----------
    t : float
        Time. (This Hamiltonian is independent of time.)

    u : ndarray, shape(n,)
        Points in phase space.

    PARAMETERS : list of floats, optional
        Hamiltonian parameters [alpha, beta]. Default is [1, 1].

    Returns
    -------
    H : ndarray, shape(n,)
        Hamiltonian at points u, in phase space at time t.
    """
    x, p_x = u.T
    alpha, beta = PARAMETERS
    return 0.5*(p_x**2 - alpha*x**2 + 0.5*beta*x**4)

def HamSN1D_Hamiltonian(t, u):
    """
    Returns Hamiltonian for the 1DoF saddle-node model.
    Functional form: H = y**2/2 + x**3/3 + x**2/2, with u = (x, y).

    Parameters
    ----------
    t : float
        Time. (This Hamiltonian is independent of time.)

    u : ndarray, shape(n,)
        Points in phase space.

    Returns
    -------
    H : ndarray, shape(n,)
        Hamiltonian at points u, in phase space at time t.
    """
    x, y = u.T
    return 0.5*y*y + x**3/3 + 0.5*x*x

def HenonHeiles_Hamiltonian(t, u):
    """
    Hamiltonian for 2DoF Henon-Heiles system.
    Functional form: H = 1/2(p_x**2 + p_y**2) + 1/2(x**2 + y**2) + yx**2 - y**3/3, with u = (x, y, p_x, p_y).

    Parameters
    ----------
    t : float
        Time. (This Hamiltonian is independent of time.)

    u : ndarray, shape(n,)
        Points in phase space.

    Returns
    -------
    H : ndarray, shape(n,)
        Hamiltonian at points u, in phase space at time t.

    """
    points_positions = u.T[:2]
    points_momenta = u.T[2:4]
    x, y = points_positions
    p_x, p_y = points_momenta
    # Potential energy function
    H = 0.5*(x**2 + y**2) + (y * x**2) - (1/3)*y**3 + 0.5*(p_x**2 + p_y**2)
    return H

def HenonHeiles_potential(positions, PARAMETERS = None):
    """
    Potential energy function for 2DoF Henon-Heiles system.
    Functional form: V = 1/2(x**2 + y**2) + yx**2 - y**3/3, with positions = (x, y).

    Parameters
    ----------
    positions : ndarray, shape(n,)
        Points in phase space.

    Returns
    -------
    V : ndarray, shape(n,)
        Potential energy at points u.

    """
    x, y = positions.T
    V = (1/2)*(x**2 + y**2) + (y * x**2) - (1/3)*y**3
    return V

def quadratic_normal_form_saddlecenter_ham(t, u):
    """
    Hamiltonian for a 2DoF Saddle.
    Functional form: H = 1/2(p_x**2 + p_y**2) + 1/2(y**2 - x**2), with u = (x, y, p_x, p_y).

    Parameters
    ----------
    t : float
        Time. (This Hamiltonian is independent of time.)

    u : ndarray, shape(n,)
        Points in phase space.

    Returns
    -------
    H : ndarray, shape(n,)
        Hamiltonian at points u, in phase space at time t.

    """
    points_positions = u.T[:2]
    points_momenta = u.T[2:4]
    x, y = points_positions
    p_x, p_y = points_momenta
    H = 0.5*(p_x**2 + p_y**2 + y**2 - x**2)
    return H

def NFSaddle_potential(positions, PARAMETERS = None):
    """
    Potential energy function for a 2DoF Saddle.
    Functional form: V = 1/2(y**2 - x**2), with positions = (x, y).
    ----------
    positions : ndarray, shape(n,)
        Points in phase space.

    Returns
    -------
    V : ndarray, shape(n,)
        Potential energy at points u.
    """
    x, y = positions.T
    V = 0.5*(y*y - x*x)
    return V

def quadratic_normal_form_saddlecentercenter_ham(t, u):
    """
    Hamiltonian for a 3DoF Saddle.
    Functional form: H = 1/2(p_x**2 + p_y**2 + p_z**2) + 1/2(y**2 + z**2 - x**2), with u = (x, y, z, p_x, p_y, p_z).

    Parameters
    ----------
    t : float
        Time. (This Hamiltonian is independent of time.)

    u : ndarray, shape(n,)
        Points in phase space.

    Returns
    -------
    H : ndarray, shape(n,)
        Hamiltonian at points u, in phase space at time t.

    """
    points_positions = u.T[:3]
    points_momenta = u.T[3:6]
    x, y, z = points_positions
    p_x, p_y, p_z = points_momenta
    H = 0.5*(p_x**2 + p_y**2 + p_z**2 + y**2 +z**2 - x**2)
    return H

def kinetic_squares(t,points_momenta):
    """
    Kinetic energy of the form sum-of-squares.

    Parameters
    ----------
    t : float
        Time. (This Hamiltonian is independent of time.)

    points_momenta : ndarray, shape(n,)
        Points in phase space.

    Returns
    -------
    T : ndarray, shape(n,)
        Kinetic energy at points u, in phase space at time t.

    """
    return 0.5* np.sum(points_momenta**2, axis=1)

def Hamiltonian_from_potential(potential):
    """
    Forms a Hamiltonian by combining input potential function with kinetic energy 'kinetic_squares'.

    Parameters
    ----------
    potential : function

    Returns
    -------
    Hamiltonian : function
        Sum of 'potential' and 'kinetic_squares'.

    """
    def Hamiltonian(t,u):
        N_dim = int(u.shape[1]/2)
        points_positions = u[:,:N_dim]
        points_momenta = u[:,N_dim:2*N_dim]
        return potential(points_positions) + kinetic_squares(t,points_momenta)
    return Hamiltonian


__author__ = 'Broncio Aguilar-Sanjuan, Victor-Jose Garcia-Garrido, Vladimir Krajnak, Shibabrat Naik'
__status__ = 'Development'
