from qutip import *

import numpy as np

import matplotlib.pyplot as plt

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

def collision_model_dynamics(N, rho_S_0, rho_A, expec_list, H_A, H_S, H_int, t, gamma, dt):
    """Simulates collision model dynamics and returns expectation values at each step"""
    state = tensor(rho_S_0, rho_A)
    





bfield = 0.5
J = 1.
gamma = 1.
Delta = 0.

num_collision = 400
collision_duration = 0.5
timestep = 0.01

num_spins = 4

psi_0 = rand_ket_haar(2 ** num_spins)
rho_S_0 = psi_0 * psi_0.dag()

rho_A = basis(2,0) * basis(2,0).dag()

expec_list = [tensor([sigmax() if k == i else qeye(2) for k in range(num_spins)]) for i in range(num_spins)]

H_S = -Hxyz(num_spins, J, J*gamma, J*Delta, bfield, periodic=True)

print(expec_list[0])