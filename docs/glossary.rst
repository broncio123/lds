========================
Glossary of Common Terms
========================

This glossary gives the mathematical definition of dynamical systems concepts used in the software for computing Lagrangian descriptors. More detailed description can be found in books listed in the references; Ref. [wiggins2003]_. 


.. glossary::

Dynamical system
----------------

   A system that changes in time; usually described by differential equations (continuous time) or difference equations (sometimes called *maps*) (discrete time), or, possibly, some combination of the two.


Vector field, map, phase space, trajectory
------------------------------------------

   The equation of the form

   .. math::
      \begin{equation}
      \dot{x} = f(x,t; \mu)
      \end{equation}

   with :math:`x \in U \subset \mathbb{R}^n, t \in \mathbb{R}^1`, and :math:`\mu \in V \subset \mathbb{R}^p` where :math:`U, V` are open sets in :math:`\mathbb{R}^n, \mathbb{R}^p`, respectively, is an **ordinary differential equation** and :math:`f` is a **vector field**. The above form is referred to as a non-autonomous or time-dependent ODE.

   The equation of the form

   .. math::
      \begin{equation}
      x \mapsto g(x; \mu)
      \end{equation}

   as a **map or difference equation**. 

   We refer to the space of dependent variables as the **phase space**. Both of these form of equations are referred to as dynamical system. Solution of the ordinary differential equations, :math:`x(t,t_0,x_0)` will be referred to as the **trajectory or phase curve** through the point :math:`x_0` at :math:`t = t_0`.

   The goal of dynamical systems analysis is to understand the geometry of solution curves in phase space. 



Equilibrium solutions, linear stability
---------------------------------------

   We give the definition of equilibrium solutions and their stability for an autonomous vector field of the form

   .. math::
      \begin{equation}
      \dot{x} = f(x; \mu), \qquad x \in \mathbb{R}^n
      \end{equation}

   The **equilibrium solution** is a point :math:`\bar{x} \in \mathbb{R}^n` such that 
   
   .. math::
      \begin{equation}
      f(\bar{x}; \mu) = 0.
      \end{equation}

   This constant solution of a vector field, that is a point in phase space where the vector field is zero, is also referred to as a "fixed point", "steady state", "stationary point", "rest point".

   Cautionary note: In case of non-autonomous vector fields, this definition of equilibrium solutions do not apply to the "frozen" (fixing :math:`t = \bar{t}`) vector field. Further explanation based on an example is given in the Chapter 1 of Ref. [wiggins2003]_ 

   We give the definition for the type of linear stability which is most relevant for this software.
   
   For a linear, autonomous vector field on :math:`\mathbb{R}^n`:
   
   .. math::
      \begin{equation}
      \dot{x} = A x, \qquad x(0) = x_0, \qquad x \in \mathbb{R}^n
      \end{equation}

   If :math:`A` has no eigenvalues with zero real part, the linear stability of the origin is determined by the real part of the eigenvalues :math:`A`. 

   If all of the real parts of the eigenvalues are strictly less than zero, then the origin is asymptotically **stable**. If at least one of the eigenvalues of :math:`A` has real part strictly larger than zero, then the origin is **unstable**.
   
   The origin of is said to be **hyperbolic** if none of the real parts of the eigenvalues of :math:`A` have zero real parts. Hyperbolic equilibria of linear, autonomous vector fields on :math:`\mathbb{R}^N` can be either sinks, sources, or saddles.
   
   More details on classification of stability of equilibria can be found in Ref. [wiggins2017]_.


Periodic orbits and invariant manifolds
---------------------------------------

   A solution :math:`x` of a dynamical system passing through the point :math:`x(0)=x_0` is said to be a **periodic solution** of period :math:`T` if there exists :math:`T > 0` such that :math:`x(t) = x(t + T)` for all :math:`t \in \mathbb{R}`.

   In a discrete system, an orbit of a point :math:`x_0 \in \mathbb{R}^n` is said to be a **periodic orbit** of period :math:`k` if :math:`g^k(x_0) = x_0`.

   In the context of this software, it is sufficient to know that a **manifold** is a set which *locally* has the structure of Euclidean space. For the "typical" dynamical systems in applications, a manifold is a :math:`m`-dimensional surface embedded in :math:`\mathbb{R}^n`.

   **Invariant manifolds** are a set of points, which is also a **manifold**, in phase space such that trajectories with initial conditions in the set remain in the set forever. In the context of this software, we are interested in invariant manifolds of hyperbolic equilibrium points or hyperbolic periodic orbits. Thus, when trajectories on an invariant manifold asymptotically approach the equilibrium point or the periodic orbit for :math:`t \rightarrow \infty` (or :math:`t \rightarrow -\infty`), the invariant manifold is called a **stable** (or **unstable**) invariant manifold. Mathematical definitions of these manifolds can be found in Chapter 3 of [wiggins2003]_ and Chapter 6 of [wiggins2017]_.

   Lagrangian descriptors identify stable (or unstable) invariant manifolds when the initial conditions are integrated for the time interval :math:`[t_0, t_0 + \tau]` (or :math:`[t_0, t_0 - \tau]`), respectively. 


==========
References
==========
   
.. [wiggins2003] Wiggins, S., (2003) Introduction to applied nonlinear dynamical systems and chaos (Vol. 2). Springer Science & Business Media.

.. [wiggins2017] Wiggins, S. (2017) Ordinary Differential Equations (Ver. 2), https://figshare.com/articles/Ordinary_Differential_Equations/5311612,  Figshare.



