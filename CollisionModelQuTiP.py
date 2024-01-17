from qutip import *

import numpy as np

import matplotlib.pyplot as plt

"""Switched to Jupyter notebook"""

def Hxyz(N, Jx, Jy, Jz, b, periodic=False):
    """Generates an XYZ spin chain Hamiltonian with N spins"""
    J_list = [Jx, Jy, Jz]
    sigma_list = [sigmax(), sigmay(), sigmaz()]
    if N <= 1:
        return qeye(2**N)
    Ham = 0
    for i in range(N-1):
        for j in range(3):
            spin_list = [sigma_list[j] if k == i or k == i+1 else qeye(2) for k in range(N)]
            Ham += J_list[j]*tensor(spin_list)
    for i in range(N):
        spin_list = [sigmaz() if k == i else qeye(2) for k in range(N)]
        Ham += b*tensor(spin_list)
    if periodic:
        for j in range(3):
            spin_list = [sigma_list[j] if k == 0 or k == N-1 else qeye(2) for k in range(N)]
            Ham += J_list[j]*tensor(spin_list)
    return Ham

def collision_model_dynamics(num_col, rho_S_0, rho_A, expec_list_S, expec_list_SA, H_S, H_int, t, gamma, dt):
    """Simulates collision model dynamics and returns expectation values at each step"""
    rho_S = rho_S_0
    state = tensor(rho_S_0, rho_A)
    H_col = tensor(H_S, qeye(2)) + H_int
    free_steps = np.ceil(np.random.exponential(gamma/dt, size = num_col))
    col_step = int(np.ceil(t/dt))

    times_free = np.linspace(dt, free_steps[0]*dt, int(free_steps[0]))
    times_col = np.linspace(dt, col_step*dt, col_step)

    res_free = mesolve(H_S, rho_S, times_free, [], expec_list_S)
    res_col = mesolve(H_col, state, times_col, [], expec_list_SA)
    return res_free, res_col

    





bfield = 0.5
J = 1.
gamma = 1.
Delta = 0.

num_collision = 400
collision_duration = 0.5
timestep = 0.01

num_spins = 4

psi_0 = tensor([rand_ket_haar(2) for _ in range(num_spins)])
rho_S_0 = psi_0 * psi_0.dag()
rho_A = basis(2,0) * basis(2,0).dag()

expec_list_S = [tensor([sigmax() if k == i else qeye(2) for k in range(num_spins)]) 
              for i in range(num_spins)]
expec_list_SA = [tensor([sigmax() if k == i else qeye(2) for k in range(num_spins )] + [qeye(2)]) 
              for i in range(num_spins)]

H_S = -Hxyz(num_spins, J, J*gamma, J*Delta, bfield, periodic=True)
H_SA = 0.5*(tensor(sigmax(), *[qeye(2) for _ in range(num_spins - 1)], sigmax()) 
          + tensor(sigmay(), *[qeye(2) for _ in range(num_spins - 1)], sigmay()))

out_free, out_col = collision_model_dynamics(num_collision, rho_S_0, rho_A, expec_list_S, expec_list_SA, 
                               H_S, H_SA, collision_duration, gamma, timestep)

print(out_free.states[-1])
#print(np.average(out), gamma/timestep)