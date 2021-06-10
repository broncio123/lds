# Two degrees of freedom quadratic normal form Hamiltonian 

The two degrees of freedom quadratic normal form Hamiltonian gives a linear system with an index one saddle. 

## Uncoupled system

.. math::
   \begin{equation}
   H_2(q_1, q_2, p_1, p_2) = \underbrace{\frac{\lambda}{2} (p_1^2 - q_1^2)}_\text{$H_r$} + \underbrace{\frac{\omega_2}{2}(q_2^2 + p^2_2)}_\text{$H_b$}\quad \lambda,\omega_2 > 0 \label{eqn:ham_nf_2dof}
   \end{equation}

with the corresponding vector field given by

.. math::
   \begin{equation}
   \begin{aligned}
   \dot{q}_1 = & \frac{\partial H_2}{\partial p_1} &= \lambda p_1, \\
   \dot{q}_2 = & \frac{\partial H_2}{\partial p_2} &= \omega_2 p_2,\\
   \dot{p}_1 = & -\frac{\partial H_2}{\partial q_1} &= \lambda q_1,\\
   \dot{p}_2 = & -\frac{\partial H_2}{\partial q_2} &= -\omega_2 q_2\\
   \end{aligned}
   \label{eqn:eom_nf_2dof}
    \end{equation}

## Coupled system


