import numpy as np
from numpy import cos, sin, exp, pi

exi = lambda p: exp(1j*p)
class ArbitraryRetarder():
    def __init__(self, connection = None, name = "arb_1", theta = 0.0, eta = 0.0, phi = 0.0):
        self.connection = connection
        self.name = name
        self.theta = theta
        self.eta = eta
        self.phi = phi

    def polarize_state(self, state):
        return -1

    def pol_matrix(self, phi: float, theta: float, eta: float):
        c = np.cos(theta)
        s = np.sin(theta)
        gp = exi(-eta / 2)

        off_diag_factor = (1 - exi(eta)) * c * s
        J00 = c ** 2 + exi(eta) * (s ** 2)
        J01 = off_diag_factor * exi(-1 * phi)
        J10 = off_diag_factor * exi(phi)
        J11 = exi(eta) * c ** 2 + (s ** 2)

        J00 *= gp
        J01 *= gp
        J10 *= gp
        J11 *= gp
        matrix = [[J00,J01],[J10,J11]]

        return matrix
    def pol_piecewise(self,):
        matrix = self.pol_matrix(phi=self.phi, theta = self.theta, eta = self.eta)

        return matrix[0][0], matrix[0][1], matrix[1][0], matrix[1][1]

class WavePlate(ArbitraryRetarder):
    def __init__(self, connection = None, name = "plate_1", plate_val = 1/4, eta = 0, theta = 0):
        super().__init__(connection = connection, name = name, eta = plate_val*2*pi , phi = 0, theta = theta)

        if connection is not None:
            self.theta = self._position
        else:
            self.theta = self.position_dry

    def position_dry(self):
        return self.theta
    def _position(self):
        return self.connection.get_position()



waveplate = WavePlate(eta=1/4, theta=2)

print(waveplate.pol_piecewise())

