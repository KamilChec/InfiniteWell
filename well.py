import numpy as np
import sys
import math


class Box:
    deltaX = 0

    def __init__(self, L=0, N=0, tau=0):
        self.L = L
        self.N = N
        self.tau = tau

    def read_parameters(self, file_name):
        with open(file_name, 'r') as inp:
            parameters = inp.read().splitlines()
        self.L = int(parameters[0])
        self.N = int(parameters[1])
        self.tau = float(parameters[2])

    def print_data(self, n, x, e, i):
        if i == 0:
            out_file = open("out.txt", "w")
        else:
            out_file = open("out.txt", "a")
        out_file.write(str(n) + "\t" + str(x) + "\t" + str(e) + "\n")

    def start_position(self):
        self.deltaX = 1. / self.N
        psi_I = 0
        psi = []
        for k in range(self.N):
            x_k = (k * self.deltaX)
            psi_R = math.sqrt(2) * math.sin(math.pi * x_k)
            psi.append([psi_R, psi_I])
        return Particle(psi)

    def hamiltonian_R(self, psi_R, psi_R_b, psi_R_a):  # dorobrobic k
        return -0.5 * (psi_R_a + psi_R_b - 2 * psi_R) / self.deltaX ** 2

    def hamiltonian_I(self, psi_I, psi_I_b, psi_I_a):  # dorobrobic k
        return -0.5 * (psi_I_a + psi_I_b - 2 * psi_I) / self.deltaX ** 2

    def start_hamiltonian(self, particle):
        hamiltonian = [[0, 0]]
        for k in range(self.N):
            if k != 0 and k != self.N - 1:
                hamiltonian.append([self.hamiltonian_R(particle.psi_R(k), particle.psi_R(k-1), particle.psi_R(k+1)),
                                    self.hamiltonian_I(particle.psi_I(k), particle.psi_I(k-1), particle.psi_I(k+1))])
        hamiltonian.append([0, 0])
        particle.hamiltonian = hamiltonian

    def count_hamiltonian_R(self, particle, k):
        if k == 0 or k == self.N - 1:
            particle.setHamiltonian_R(0, k)
        else:
            particle.setHamiltonian_R(self.hamiltonian_R(particle.psi_R(k), particle.psi_R(k-1), particle.psi_R(k+1)),
                                      k)

    def count_hamiltonian_I(self, particle, k):
        if k == 0 or k == self.N - 1:
            particle.setHamiltonian_I(0, k)
        else:
            particle.setHamiltonian_I(self.hamiltonian_I(particle.psi_I(k), particle.psi_I(k-1), particle.psi_I(k+1)),
                                      k)

    def simulation(self, particle):
        for i in range(4):
            for k in range(self.N):
                psi_R_halfTau = particle.psi_R(k) + particle.hamiltonian_I(k) * self.tau / 2
                particle.setPsi_R(psi_R_halfTau, k)
            for k in range(self.N):
                self.count_hamiltonian_R(particle, k)
            for k in range(self.N):
                psi_I = particle.psi_I(k) - particle.hamiltonian_R(k) * self.tau
                particle.setPsi_I(psi_I, k)
            for k in range(self.N):
                self.count_hamiltonian_I(particle, k)
            for k in range(self.N):
                psi_R = psi_R_halfTau + particle.hamiltonian_I(k) * self.tau / 2
                particle.setPsi_R(psi_R, k)
            # if i % 100 == 0:
            #     squareSum_psi = 0
            #     squareSum_psi_dotX = 0
            #     hamiltonianSum_psi = 0
            #     for k in range(self.N):
            #         squareSum_psi += particle.psi_R(k)**2 + particle.psi_I(k)**2
            #         x_k = (k * self.deltaX)
            #         squareSum_psi_dotX += x_k * (particle.psi_R(k)**2 + particle.psi_I(k)**2)
            #         hamiltonianSum_psi += particle.psi_R(k) * particle.hamiltonian_R(k) + \
            #                               particle.psi_I(k) * particle.hamiltonian_I(k)
            #     n = self.deltaX * squareSum_psi
            #     x = self.deltaX * squareSum_psi_dotX
            #     e = self.deltaX * hamiltonianSum_psi
            #     self.print_data(n, x, e, i)


class Particle:
    def __init__(self, psi=[], hamiltonian=[]):
        self.psi = psi
        self.hamiltonian = hamiltonian

    def setHamiltonian(self, hamiltonian, k):
        self.hamiltonian[k] = hamiltonian

    def setHamiltonian_R(self, hamiltonian_R, k):
        self.hamiltonian[k][0] = hamiltonian_R

    def setHamiltonian_I(self, hamiltonian_I, k):
        self.hamiltonian[k][1] = hamiltonian_I

    def setPsi_R(self, psi_R, k):
        self.psi[k][0] = psi_R

    def setPsi_I(self, psi_I, k):
        self.psi[k][1] = psi_I

    def psi_R(self, k):
        return self.psi[k][0]

    def psi_I(self, k):
        return self.psi[k][1]

    def hamiltonian_R(self, k):
        return self.hamiltonian[k][0]

    def hamiltonian_I(self, k):
        return self.hamiltonian[k][1]


def main():
    try:
        sys.argv[1]
    except IndexError:
        print("insert file into argument")
        return 0
    box = Box()
    box.read_parameters(sys.argv[1])
    particle = box.start_position()
    box.start_hamiltonian(particle)
    box.simulation(particle)
    print(particle.psi)


if __name__ == "__main__":
    main()
