import numpy as np

def HamCenter1D(t, u, PARAMETERS = [1]):
    """
    Returns vector field for a 1DoF centre at time t, for an array of points in phase space.
    Number of model parameters: 1 . PARAMETERS = [omega]
    Functional form: v = (omega*y, - omega*x), with u = (x, y)

    Parameters
    ----------
    t : float
        Time. (This vector field is independent of time.)

    u : ndarray, shape(n,)
        Points in phase space.

    PARAMETERS : list of floats
        Vector field parameters.

    Returns
    -------
    v : ndarray, shape(n,)
        Vector field at points u and time t.
    """
    x, y = u.T
    # Hamiltonian Model Parameter
    omega, = PARAMETERS
    v = np.column_stack([ omega * y, - omega * x])
    return v

def HamSaddle1D(t, u, PARAMETERS = [1]):
    """
    Returns vector field for a 1DoF saddle at time t, for an array of points in phase space.
    Number of model parameters: 1 . PARAMETERS = [lamda]
    Functional form: v = (lamda*y, lamda*x), with u = (x, y)

    Parameters
    ----------
    t : float
        Time. (This vector field is independent of time.)

    u : ndarray, shape(n,)
        Points in phase space.

    PARAMETERS : list of floats
        Vector field parameters.

    Returns
    -------
    v : ndarray, shape(n,)
        Vector field at points u and time t.
    """
    x, y = u.T
    # Hamiltonian Model Parameter
    lamda, = PARAMETERS
    v = np.column_stack([ lamda * y, lamda * x])
    return v

def Duffing1D(t, u, PARAMETERS = [1, 1]):
    """
    Returns vector field for the Duffing oscillator.
    Number of model parameters: 2 . PARAMETERS = [alpha, beta]
    Functional form: v = (y, alpha*x - beta*x**3), with u = (x, y)

    Parameters
    ----------
    t : float
        Time. (This vector field is independent of time.)

    u : ndarray, shape(n,)
        Points in phase space.

    PARAMETERS : list of floats, optional
        Vector field parameters [alpha, beta]. Default is [1, 1].

    Returns
    -------
    v : ndarray, shape(n,)
        Vector field at points u and time t..
    """
    x, y = u.T
    # Hamiltonian Model Parameter
    alpha, beta = PARAMETERS
    v = np.column_stack([ y, alpha*x - beta*x**3])
    return v

def HamSN1D(t, u, PARAMETERS = [None]):
    """
    Returns vector field for the 1DoF saddle-node model.
    Number of model parameters: 0 . PARAMETERS = [None]
    Functional form: v = (y, -x -x**2), with u = (x, y)

    Parameters
    ----------
    t : float
        Time. (This vector field is independent of time.)

    u : ndarray, shape(n,)
        Points in phase space.

    PARAMETERS : list of floats
        Vector field parameters.

    Returns
    -------
    v : ndarray, shape(n,)
        Vector field at points u and time t.
    """
    x, y = u.T
    # Hamiltonian Model Parameter
    v = np.column_stack([ y, -x -x**2])
    return v

def forcing(t, u, perturbation_params = [0, 1, 0.15, 0.5]):
    """
    Returns vector field for a perturbation.
    Number of model parameters: 3. perturbation_params = [phase_shift, perturbation_type, amplitude, frequency]
    Functional form: v = (, ), with u = (x, y)

    Parameters
    ----------
    t : float
        Time. (This vector field is independent of time.)

    u : ndarray, shape(n,)
        Points in phase space.

    perturbation_params : list of floats, [phase_shift, perturbation_type, amplitude, frequency]
        Perturbation parameters.

    Returns
    -------
    v : ndarray, shape(n,)
        Vector field at points u and time t.
    """
    x, y = u.T
    perturbation = np.zeros(u.shape)

    # Perturbation parameters
    phase_shift, perturbation_type, amplitude, freq = perturbation_params
    time = t + phase_shift

    if perturbation_type == 1:
        perturbation = perturbation + np.array([0, amplitude * np.sin(freq*time)])
    elif perturbation_type == 2:
        perturbation = perturbation + np.array([0, amplitude * np.sin(freq*time)/np.cosh(time)])

    return perturbation

def HenonHeiles_vector_field(t, u):
    """
    Returns Henon-Heiles vector field (2DoF).
    Functional form: v = (p_x, p_y, -x - 2*x*y, -x**2 -y + y**2), with u = (x, y, p_x, p_y)

    Parameters
    ----------
    t : float
        Time. (This vector field is independent of time.)

    u : ndarray, shape(n,)
        Points in phase space.

    Returns
    -------
    v : ndarray, shape(n,)
        Vector field at points u and time t.
    """
    points_positions = u.T[:2]
    points_momenta = u.T[2:4]
    x, y = points_positions
    p_x, p_y = points_momenta

    # Vector field defintion
    v_x   =  p_x
    v_y   =  p_y
    v_p_x = -x - 2*x*y
    v_p_y = -x**2 -y + y**2
    v = np.column_stack([v_x, v_y, v_p_x, v_p_y])
    return v

def quadratic_normalform_saddlecenter(t, u, PARAMETERS = [1,1]):
    """
    Returns vector field for a 2D index-1 saddle.
    Functional form: v = (p_x, p_y, x, -y), with u = (x, y, p_x, p_y)

    Parameters
    ----------
    t : float
        Time. (This vector field is independent of time.)

    u : ndarray, shape(n,)
        Points in phase space.

    PARAMETERS : list of floats
        Vector field parameters.

    Returns
    -------
    v : ndarray, shape(n,)
        Vector field at points u and time t.
    """
    N_dim = u.shape[-1]
    points_positions = u.T[:int(N_dim/2)]
    points_momenta = u.T[int(N_dim/2):]
    x, y = points_positions
    p_x, p_y = points_momenta

    # Hamiltonian Model Parameter
    # None

    # Vector field defintion
    v_x   =  p_x
    v_y   =  p_y
    v_p_x = PARAMETERS[0]*x
    v_p_y = -PARAMETERS[1]*y
    v = np.column_stack([v_x, v_y, v_p_x, v_p_y])
    return v


def DoubleGyre(t, u, PARAMETERS = [0, 0.25, 2*np.pi, 0, 0, 1, 0.25]):
    """
    Returns 2D Double Gyre vector field at time t, for an array of points in phase space.
    Number of model parameters: 6 . PARAMETERS = [phase_shift, A, phi, psi, mu, s, epsilon]
    Functional form: 
    
    vx = -pi*A*sin(pi*f(t + phase_shift, x)/s)*cos(pi*y/s) - mu*x
    vy =  pi*A*cos(pi*f(t + phase_shift, x)/s)*sin(pi*y/s)*df(t + phase_shift,x)/dx - mu*y
    
    with
    
    f(t, x)    = epsilon*sin(phi*t + psi)*x**2 + (1 - 2*epsilon*sin(phi*t + psi))*x
    df/dx(t,x) = 2*epsilon*sin(phi*t + psi)*x + (1 - 2*epsilon*sin(phi*t + psi))
    u = (x, y)
    Parameters
    ----------
    t : float
        fixed time-point of vector field, for all points in phase space.
    u : array_like, shape(n,)
        points in phase space to determine vector field at time t.
    PARAMETERS : list of floats
        vector field parameters
    Returns
    -------
    v : array_like, shape(n,)
        vector field corresponding to points u, in phase space at time t
    """
    x, y = u.T
    # model parameter
    phase_shift, A, phi, psi, mu, s, epsilon = PARAMETERS
    
    time = t + phase_shift
    # vector field components
    def f(t, x): return epsilon*np.sin(phi*t + psi)*x**2 + (1-2*epsilon*np.sin(phi*t + psi))*x
    def df_dx(t,x): return 2*epsilon*np.sin(phi*t + psi)*x + (1-2*epsilon*np.sin(phi*t + psi))
    v_x = -np.pi*A*np.sin(np.pi*f(time, x)/s)*np.cos(np.pi*y/s) - mu*x
    v_y =  np.pi*A*np.cos(np.pi*f(time, x)/s)*np.sin(np.pi*y/s)*df_dx(time,x) - mu*y
    v   = np.column_stack([v_x, v_y])
    return v


def quadratic_normalform_saddlecentercenter(t, u, PARAMETERS = [1,1,1]):
    """
    Returns vector field for a 3D index-1 saddle.
    Functional form: v = (p_x, p_y, p_z, x, -y, -z), with u = (x, y, z, p_x, p_y, p_z)

    Parameters
    ----------
    t : float
        Time. (This vector field is independent of time.)

    u : ndarray, shape(n,)
        Points in phase space.

    PARAMETERS : list of floats
        Vector field parameters.

    Returns
    -------
    v : ndarray, shape(n,)
        Vector field at points u and time t.
    """
    N_dim = u.shape[-1]
    points_positions = u.T[:int(N_dim/2)]
    points_momenta = u.T[int(N_dim/2):]
    x, y, z = points_positions
    p_x, p_y, p_z = points_momenta

    # Hamiltonian Model Parameter
    # None

    # Vector field defintion
    v_x   =  p_x
    v_y   =  p_y
    v_z   =  p_z
    v_p_x = PARAMETERS[0]*x
    v_p_y = -PARAMETERS[1]*y
    v_p_z = -PARAMETERS[2]*z
    v = np.column_stack([v_x, v_y, v_z, v_p_x, v_p_y, v_p_z])
    return v

__author__ = 'Broncio Aguilar-Sanjuan, Victor-Jose Garcia-Garrido, Vladimir Krajnak, Shibabrat Naik'
__status__ = 'Development'
