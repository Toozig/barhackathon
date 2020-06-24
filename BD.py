import argparse
# import seaborn as sbn
import matplotlib.pyplot as plt

import numpy as np

DEFAULT_SHAPE = (5, 5)
DT = 0.1
BD_FACTOR = 6
FORCE_VEC = np.asarray([-0.1, -0.5])
UPPER_BOUND = 4.5
LOWER_BOUND = -0.5


class BD:
    def __init__(self, n, kt, dt, d, x0, force_f):
        self.__n = n
        self.__kt = kt
        self.__dt = dt
        self.__d = d
        self.__x0 = x0
        self.__force_func = force_f
        self.__r = 0


    def BD_algorithm():
        """
        The Brownian dynamic algorithm implementation,  as described in our paper.
        :param n: Number of iterations
        :param kT: kT
        :return: The configuration after n iterations.
        """
        result = np.zeros((self.__n, len(self.__x0)))
        result[0] = self.__x0
        # x = np.random.uniform(LOWER_BOUND, UPPER_BOUND)
        # y = np.random.uniform(LOWER_BOUND, UPPER_BOUND)
        cur = self.__x0
        print(cur)
        while i < range(self.__n):
            r = np.random.normal(0, 1, 2)
            # result[round(cur[0]), round(cur[1])] += 1
            addition = self.__dt * self.__force_func(cur) + r * np.sqrt(BD_FACTOR * self.__kt * self.__dt)
            temp = cur + addition

            # todo : try to write the function again
            if LOWER_BOUND < temp[0] < UPPER_BOUND and LOWER_BOUND < temp[1] < UPPER_BOUND:
                cur = temp
                result[i] = cur
                i += 1
            else:
                temp -= addition

        return result

    def Brownian_dynamics(n, R, t, kT):
        grid = np.zeros((5, 5))
        x = np.random.uniform(-0.5, 4.5, 2)
        force = np.array([-1, -0.5])
        # cur_x = np.random.uniform(-0.5, 4.5)
        # cur_y = np.random.uniform(-0.5, 4.5)
        for i in range(n):
            grid[int(round(x[0])), int(round(x[1]))] += 1
            f_x = force @ x
            r_vec = np.random.normal(size=2)
            x = x + t * f_x + np.sqrt(6 * kT * t) * r_vec
            # if x[0] < -0.5: x[0] = 2 + (-0.5 + x[0])
            # if x[1] < -0.5: x[1] = 2 + (-0.5 + x[1])
            # if x[0] > 4.5: x[0] = 2 + (4.5 - x[0])
            # if x[1] > 4.5: x[1] = 2 + (4.5 - x[1])
            if x[0] < -0.5: x[0] = -0.5
            if x[1] < -0.5: x[1] = -0.5
            if x[0] > 4.5: x[0] = 4.5
            if x[1] > 4.5: x[1] = 4.5


def get_cmdline_parser():
    parser = argparse.ArgumentParser(
        description='Run BD on a n configuration space.')
    parser.add_argument('n', type=int, default=1000, nargs='?',
                        help='number of iterations of MCMC optimization')
    parser.add_argument('kT', type=float, default=1.0, nargs='?',
                        help='kT - the denominator for the metropolis criterion'
                             ' (Boltzmann constant times temperature)')
    return parser


if __name__ == '__main__':
    parser = get_cmdline_parser().parse_args()
    plt.figure()
    for i in [0.1, 1]:
        result = BD_algorithm(1000000, i)
        print(result)
        sbn.heatmap(result / np.sum(result), annot=True, fmt=".2f", linewidths=.5)
        plt.title("Configuration frequency heatmap for BD algorithm (kT = %.2f)" % i)
        plt.legend("frequency of visits")
        plt.show()
        plt.savefig("BD_alg_kt_%.2f.png" % i)
