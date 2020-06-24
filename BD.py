import argparse
import seaborn as sbn
import matplotlib.pyplot as plt

import numpy as np

# DEFAULT_SHAPE = (5, 5)
DT = 0.1
BD_FACTOR = 6
# FORCE_VEC = np.asarray([-0.1, -0.5])
UPPER_BOUND = 4.5
LOWER_BOUND = -0.5
MAX_DEPTH = 5


class BD:
    def __init__(self, n, kt, dt, x0,dim):
        """
        :param n: iteration number
        :param kt: kT
        :param dt: time interval length
        :param d: D
        :param x0: initial state
        :param force_f: function that calculate the forces that work on the proteins
        :param dim: the dimension of the rectangle that we looking at (2 or 3)
        """
        self.__n = n
        self.__kt = kt
        self.__dt = dt
        self.__x0 = x0
        self.__force_func = None
        self.__dim = dim
        self.__r = 0

    def __vec_in_rec(self, vec, x_0, y_0, z_0, x_1, y_1, z_1):
        """
        The function checks if the vactor values are in the rectangle
        :param vec: the vector to check
        :param x_0: lower x bound
        :param y_0: lower y bound
        :param z_0: lower z bound
        :param x_1: upper x bound
        :param y_1: upper y bound
        :param z_1: upper z bound
        :return:
        """
        for i in vec:
            if self.__dim == 3:
                if not ((x_0 <= i[0] <= x_1) or (y_0 <= i[1] <= y_1) or (z_0 <= i[2] <= z_1)):
                    return False
            elif self.__dim == 2:
                if not ((x_0 <= i[0] <= x_1) or (y_0 <= i[1] <= y_1)):
                    return False
        return True

    def __build_distance_matrix(self, protein_vec):
        real_res = []
        for i in range(len(protein_vec)):
            res = []
            for j in range(len(protein_vec)):
                res.append(np.linalg.norm(protein_vec[i]-protein_vec[j]))
            real_res.append(res)
        return real_res

    def BD_algorithm(self):
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
        i , j = 0, 0
        while i < self.__n:
            if j > MAX_DEPTH:
                print("error, stack outside of the square")
                return
            r = np.random.normal(0, 1, 2)
            # result[round(cur[0]), round(cur[1])] += 1
            m = self.__build_distance_matrix(cur)
            addition = self.__dt * self.__force_func(m) + r * np.sqrt(BD_FACTOR * self.__kt * self.__dt)
            temp = cur + addition

            # todo : try to write the function again
            if self.__vec_in_rec(cur, LOWER_BOUND, LOWER_BOUND, 0, UPPER_BOUND, UPPER_BOUND, 0): #todo: decide the rectangle boundaries
                cur = temp
                result[i] = cur
                i += 1
                j = 0
            else:
                temp -= addition
                j += 1

        return result


# def get_cmdline_parser():
#     parser = argparse.ArgumentParser(description='Run BD on a n configuration space.')
#     parser.add_argument('n', type=int, default=1000, nargs='?',
#                         help='number of iterations of MCMC optimization')
#     parser.add_argument('kT', type=float, default=1.0, nargs='?',
#                         help='kT - the denominator for the metropolis criterion'
#                              ' (Boltzmann constant times temperature)')
#     return parser


if __name__ == '__main__':
    print("hi")
    # parser = get_cmdline_parser().parse_args()
    # plt.figure()
    # for i in [0.1, 1]:
    #     result = BD_algorithm(1000000, i)
    #     print(result)
    #     sbn.heatmap(result / np.sum(result), annot=True, fmt=".2f", linewidths=.5)
    #     plt.title("Configuration frequency heatmap for BD algorithm (kT = %.2f)" % i)
    #     plt.legend("frequency of visits")
    #     plt.show()
    #     plt.savefig("BD_alg_kt_%.2f.png" % i)
