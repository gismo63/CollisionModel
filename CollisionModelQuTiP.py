from qutip import *

import numpy as np

import matplotlib.pyplot as plt



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

expec_list = [tensor(identity(2**i), sigmax(), identity(2**(num_spins - i - 1))) for i in range(num_spins)]

print(len(expec_list))