import numpy as np
import sys
import math


class Box:
    deltaX = 0

    def __init__(self, L=0, N=0):
        self.L = L
        self.N = N

    def read_parameters(self, file_name):
        with open(file_name, 'r') as inp:
            parameters = inp.read().splitlines()
        self.L = int(parameters[0])
        self.N = int(parameters[1])

    def start_position(self):
        deltaX = 1. / self.N
        psi_R = 0
        psi_I = 0
        psi = []
        for k in range(self.N):
            x_k = (k * deltaX)
            psi_R = math.sqrt(2) * math.sin(math.pi * x_k)
            psi.append([psi_R, psi_I])
        return Particle(psi)

    def hamiltonian_R(self, psi_R, psi_R_b, psi_R_a):  # dorobrobic k
        return -1. / 2 * (psi_R_a + psi_R_b - 2 * psi_R) / self.deltaX ** 2

    def hamiltonian_I(self, psi_I, psi_I_b, psi_I_a):  # dorobrobic k
        return -1. / 2 * (psi_I_a + psi_I_b - 2 * psi_I) / self.deltaX ** 2

    def hamiltonian(self, particle):
        hamiltonian = [[0, 0]]
        for k in range(self.N):
            if k != 0 or k != self.N - 1:
                hamiltonian.append([self.hamiltonian_R(particle.psi[k][0], particle.psi[k-1][0], particle.psi[k+1][0]),
                                    self.hamiltonian_I[particle.psi[k][1], particle.psi[k-1][1], particle.psi[k+1][1]]])
        hamiltonian.append([0, 0])
        print(hamiltonian)


class Particle:
    def __init__(self, psi=[], hamiltonian=0):
        self.psi = psi
        self.hamiltonian = hamiltonian


def main():
    try:
        sys.argv[1]
    except IndexError:
        print("insert file into argument")
        return 0
    box = Box()
    box.read_parameters(sys.argv[1])
    particle = box.start_position()
    box.hamiltonian(particle)


if __name__ == "__main__":
    main()
