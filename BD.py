import argparse
import seaborn as sbn
import matplotlib.pyplot as plt

import numpy as np

DEFAULT_SHAPE = (5, 5)
DT = 0.1
BD_FACTOR = 6
FORCE_VEC = np.asarray([-0.1, -0.5])
UPPER_BOUND = 4.5
LOWER_BOUND = -0.5


def get_cmdline_parser():
    parser = argparse.ArgumentParser(
        description='Run BD on a n configuration space.')
    parser.add_argument('n', type=int, default=1000, nargs='?',
                        help='number of iterations of MCMC optimization')
    parser.add_argument('kT', type=float, default=1.0, nargs='?',
                        help='kT - the denominator for the metropolis criterion'
                             ' (Boltzmann constant times temperature)')
    return parser


def BD_algorithm(n, kT):
    """
    The Brownian dynamic algorithm implementation,  as described in our paper.
    :param n: Number of iterations
    :param kT: kT
    :return: The configuration after n iterations.
    """
    result = np.zeros(DEFAULT_SHAPE)
    x = np.random.uniform(LOWER_BOUND, UPPER_BOUND)
    y = np.random.uniform(LOWER_BOUND, UPPER_BOUND)
    cur = np.asarray([x, y])
    print(cur)
    for i in range(n):
        r = np.random.normal(0, 1, 2)
        result[round(cur[0]), round(cur[1])] += 1
        addition = DT * FORCE_VEC + r * np.sqrt(BD_FACTOR * kT * DT)
        temp = cur + addition
        if LOWER_BOUND < temp[0] < UPPER_BOUND and LOWER_BOUND < temp[1] < UPPER_BOUND:
            cur = temp
        else:
            temp -= addition
        if LOWER_BOUND < temp[0] < UPPER_BOUND and LOWER_BOUND < temp[1] < UPPER_BOUND:
            cur = temp

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
