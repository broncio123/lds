Examples
========

.. examples::


Discrete systems
----------------

1. Standard map

The standard map (kicked rotator) is a two-dimensional map used in dynamical systems to study a periodically 
kicked pendulum. Its equations of motion are given by the expressions:

.. math::
   \begin{align}
    x_{n+1} &= x_{n} + y_{n} - \dfrac{K}{2\pi} \sin(2\pi x_{n)),
    y_{n+1} &= y_{n} - \dfrac{K}{2\pi} \sin(2\pi x_{n)),
   \end{align}
where :math:`K` is the parameter that controls the forcing strength of the perturbation.
   
The inverse map is described by:

.. math::
   \begin{align}
    x_{n} = x_{n+1} - y_{n+1},
    y_{n} = y_{n+1} + \dfrac{K}{2\pi} \sin(2\pi (x_{n+1} - y_{n+1})).
   \end{align}

   
2. Hénon map

The Hénon map was introduced by Michel Hénon as a simplified model of the Poincaré section 
of the Lorenz model. The map equations are as follows:

.. math::
   \begin{align}
    x_{n+1} = a - x_{n}^2 + b y_{n},
    y_{n+1} = x_{n},
   \end{align}
   
where :math:`a,b` are the model parameters.

The inverse Hénon map is:

.. math::
   \begin{align}
    x_{n} = y_{n+1},
    y_{n} = \dfrac{x_{n+1} - a + y_{n+1}^2}{b}.
   \end{align}


Continuous systems with one degree of freedom
---------------------------------------------

1. Hamiltonian center

The Hamiltonian function:

.. math::
   H(x,p_x) = \dfrac{\omega}{2} \left( p_x^2 + x^2 \right), \label{eqn:ham_center1dof}

defines the normal form of a 1 DoF system with a center equilibrium at the origin. The associated equations of motion are:

.. math::
   \begin{align}
   \dot{x} &= \dfrac{\partial H}{\partial p_x} =  \omega p_x, \\
   \dot{p}_x &= -\dfrac{\partial H}{\partial x} = -\omega x.
   \end{align}
   

2. `Hamiltonian saddle <https://champsproject.github.io/lagrangian_descriptors/content/chapter2_1.html#one-degree-of-freedom-hyperbolic-equilibrium-point>`_

The Hamiltonian function:

.. math::
   H(x,p_x) = \dfrac{\lambda}{2} \left( p_x^2 - x^2 \right), \label{eqn:ham_saddle1dof}

defines the normal form of a 1 DoF system with a saddle equilibrium at the origin. The associated equations of motion are:

.. math::
   \begin{align}
   \dot{x} &= \dfrac{\partial H}{\partial p_x} = \lambda p_x, \\
   \dot{p}_x &= -\dfrac{\partial H}{\partial x} = \lambda x.
   \end{align}

3. Duffing oscillator 

a. Unforced

The Hamiltonian function:

.. math::
   H(x,p_x,t) = \dfrac{1}{2}p_x^2 - \dfrac{\alpha}{2}x^2 + \dfrac{\beta}{4}x^4, \label{eqn:ham_duff}

with :math:`\alpha,\beta>0`describes the Duffing oscillator with the associated equations of motion

.. math::
   \begin{align}
   \dot{x} &= \dfrac{\partial H}{\partial p_x} =  p_x, \\
   \dot{p}_x &= -\dfrac{\partial H}{\partial x} =  \alpha x - \beta x^3.
   \end{align}

b. Forced 
Time dependent Hamiltonian function

.. math::
   H(x,p_x,t) = \dfrac{1}{2}p_x^2 - \dfrac{\alpha}{2}x^2 + \dfrac{\beta}{4}x^4 - f(t) x, \label{eqn:ham_duff_forced}

defines the Duffing oscillator with time dependent forcing :math:`f(t)`. This package offers two predefined options for the external forcing, namely :math:`f(t) = A\mathrm{sech}(t)\sin(\omega t)` and :math:`f(t) = A\sin(\omega t)`. Other versions can be added manually by the user in the forcing function of the vector_fields.py file.

The corresponding equations of motion are:

.. math::
   \begin{align}
   \dot{x} &= \dfrac{\partial H}{\partial p_x} =  p_x, \\
   \dot{p}_x &= -\dfrac{\partial H}{\partial x} =  \alpha x - \beta x^3 + f(t).
   \end{align}



c. Inverted 

The inverted Duffing oscillator can be obtained from Hamiltonian ~\eqref{eqn:eqn:ham_duff}, by setting the parameters :math:`\alpha = \beta = - 1`.

4. Saddle-node Hamiltonian 

This system is defined by the Hamiltonian:

.. math::
    H(x,p_x) = \dfrac{1}{2}p_x^2 + \dfrac{1}{2}x^2 + \dfrac{1}{3}x^3, \label{eqn:ham_saddnode}

and its associated equations of motion are:

.. math::
    \begin{align}
    \dot{x} = \dfrac{\partial H}{\partial p_x} =  p_x, \\
    \dot{p}_x = -\dfrac{\partial H}{\partial x} =  -x - x^2.
    \end{align} 

5. Non-autonomous double-gyre flow

The double-gyre flow is a classical system popular in geophysical fluid dynamics. This non-autonomous two-dimensional dynamical system is defined by the equations:

.. math::
   \begin{align}
   \dot{x} &= -\pi A \sin\left(\dfrac{\pi f(x,t)}{s}\right) \cos\left(\dfrac{\pi y}{s}\right) - \mu x, \\[.2cm]
   \dot{y} &= \pi A \cos\left(\dfrac{\pi f(x,t)}{s}\right) \sin\left(\dfrac{\pi y}{s}\right) \dfrac{\partial f}{\partial x}\left(x,t\right) - \mu y,
   \end{align} 

where we have that :math:`f(x,t) = \varepsilon \sin(\omega t + \phi) x^2 + \left(1-2\varepsilon \sin(\omega t + \phi)\right) x`.

Continuous systems with two degrees of freedom
----------------------------------------------

1. `Saddle-center <https://champsproject.github.io/lagrangian_descriptors/content/chapter2_1.html#two-degrees-of-freedom-and-the-hyperbolic-periodic-orbit>`_ 

The Hamiltonian function:

.. math::
   H(x,y,p_x,p_y) = \dfrac{1}{2} \left( p_x^2 + p_y^2 + y^2 - x^2) \right),  \label{eqn:ham_saddle2dof}

is the normal form of a 2 DoF system with a saddle-center equilibrium point at the origin. The dynamics of any 2 DoF dynamical system near a potential index-1 saddle point is conjugate to this system.
The associated equations of motion are:

.. math::
   \begin{align}
   \dot{x} &= \dfrac{\partial H}{\partial p_x}  = p_x, \\
   \dot{y} &= \dfrac{\partial H}{\partial p_y}  = p_y, \\
   \dot{p}_x &= -\dfrac{\partial H}{\partial x}  = x, \\
   \dot{p}_y &= -\dfrac{\partial H}{\partial y}  = - y.
   \end{align}

2. Hénon-Heiles

The Hamiltonian for the Hénon-Heiles system is given:

.. math::
   H(x,y,p_x,p_y) = \dfrac{1}{2} \left( p_x^2 + p_y^2 \right) + \dfrac{1}{2} \left( x^2 + y^2 \right) + yx^2 - \dfrac{1}{3} y^3, \label{eqn:ham_hh}

and Hamilton's equations of motion are:

.. math::
   \begin{align}
   \dot{x} &= \dfrac{\partial H}{\partial p_x}  = p_x, \\
   \dot{y} &= \dfrac{\partial H}{\partial p_y} = p_y, \\
   \dot{p}_x &= -\dfrac{\partial H}{\partial x} = - x - 2xy, \\
   \dot{p}_y &= -\dfrac{\partial H}{\partial y} = - x^2 - y + y^2.
   \end{align}

This system is a fundamental system for studying complex dynamics. Depending on the value of total energy, it can exhibit different dynamical behaviour ranging from near-integrable to completely chaotic.

Continuous systems with three degrees of freedom
------------------------------------------------

1. `Saddle-center-center <https://champsproject.github.io/lagrangian_descriptors/content/chapter2_1.html#three-and-more-degrees-of-freedom-and-nhims>`_

The Hamiltonian function:

.. math::
   H(x,y,z,p_x,p_y,p_z) = \dfrac{1}{2} \left( p_x^2 + p_y^2+ p_z^2 - x^2 + y^2 + z^2) \right),  \label{eqn:ham_saddle3dof}

is the normal form of a 3 DoF system with a saddle-center-center equilibrium point at the origin (also referred to as an index-1 saddle).
The associated equations of motion are:

.. math::
   \begin{align}
   \dot{x} &= \dfrac{\partial H}{\partial p_x} = p_x, \\
   \dot{y} &= \dfrac{\partial H}{\partial p_y} = p_y, \\
   \dot{z} &= \dfrac{\partial H}{\partial p_z} = p_z, \\
   \dot{p}_x &= -\dfrac{\partial H}{\partial x} = x, \\
   \dot{p}_y &= -\dfrac{\partial H}{\partial y}= - y, \\
   \dot{p}_z &= -\dfrac{\partial H}{\partial z}= - z.
   \end{align}


.. automodule:: vector_fields
   :members:


