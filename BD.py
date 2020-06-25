import argparse
import seaborn as sbn
import matplotlib.pyplot as plt
import pandas as pd

import numpy as np

# DEFAULT_SHAPE = (5, 5)
DT = 0.1
BD_FACTOR = 6
# FORCE_VEC = np.asarray([-0.1, -0.5])
UPPER_BOUND = 4.5
LOWER_BOUND = -0.5
MAX_DEPTH = 5


class BD:
    def __init__(self, n, kt, dt, x0, dim, indexs, protein_vec, matrix):
        """
        :param n: iteration number
        :param kt: kT
        :param dt: time interval length
        :param d: D
        :param x0: initial state
        :param dim: the dimension of the rectangle that we looking at (2 or 3)
        """
        self.__n = n
        self.__kt = kt
        self.__dt = dt
        self.__dim = dim
        self.__r = 0
        self.__inx = indexs
        self.__protein_vec = protein_vec
        self.__matrix = matrix
        self.__x0 = self.__create_x0(x0)


    def __create_x0(self, c):
        print(c)
        l = (len(self.__protein_vec))
        new_x0 = [0] * l
        if c is not tuple:
            pairs = []
        elif type(c[0]) is int:
            pairs = list(c)
        else:
            pairs = [x for y in c for x in y]
        for i in range(len(new_x0)):
            if new_x0[i] != 0:
                continue
            x = np.random.uniform(LOWER_BOUND, UPPER_BOUND)
            y = np.random.uniform(LOWER_BOUND, UPPER_BOUND)
            if i in pairs:
                ind = pairs.index(i)
                if ind % 2 == 0:
                    new_x0[i] = [x,y]
                    new_x0[pairs[ind + 1]] = [x,y]
                else:
                    new_x0[i] = (x,y)
                    new_x0[pairs[ind - 1]] = [x,y]
            else:
                new_x0[i] = [x,y]
        return new_x0


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
                res.append(np.linalg.norm(np.asarray(protein_vec[i]) - np.asarray(protein_vec[j])))
            real_res.append(res)
        return real_res

    def BD_algorithm(self):
        """
        The Brownian dynamic algorithm implementation,  as described in our paper.
        :param n: Number of iterations
        :param kT: kT
        :return: The configuration after n iterations.
        """
        result = [0] * self.__n
        result[0] = self.__x0
        # x = np.random.uniform(LOWER_BOUND, UPPER_BOUND)
        # y = np.random.uniform(LOWER_BOUND, UPPER_BOUND)
        cur = self.__x0
        # print(cur)
        i , j = 0, 0
        while i < self.__n:
            if j > MAX_DEPTH:
                print("error, stack outside of the square")
                result[i] = cur
                j = 0
                i += 1
                continue
            r = np.random.normal(0, 1, len(cur))
            # result[round(cur[0]), round(cur[1])] += 1
            m = np.matrix(self.__build_distance_matrix(cur))
            force_vec = self.__force_func(m)
            dt = self.__dt
            kt = self.__kt
            addition = dt * force_vec + r * np.sqrt(BD_FACTOR * kt * dt)
            temp = cur
            for l in range(len(temp)):
                temp[l][0] += addition[l]
                temp[l][1] += addition[l]
            # temp = cur + addition

            # todo : try to write the function again
            if self.__vec_in_rec(temp, LOWER_BOUND, LOWER_BOUND, 0, UPPER_BOUND, UPPER_BOUND, 0): #todo: decide the rectangle boundaries
                cur = temp
                result[i] = cur
                i += 1
                j = 0
            else:
                for k in range(len(temp)):
                    temp[k][0] -= addition[k]
                    temp[k][1] -= addition[k]
                j += 1
        return result

    def __force_func(self, dist_met):
        distance = self.__matrix
        new_mat = []
        for i in range(dist_met.shape[0]):
            new_row = []
            for j in range(dist_met.shape[1]):
                x_pro = self.__inx[i]
                y_pro = self.__inx[j]
                if dist_met[i,j] != 0:
                    x = 1/ dist_met[i,j]
                else:
                    x = 0
                new_row.append(x * distance[x_pro][y_pro])
            new_mat.append(new_row)
        return np.array(new_mat).sum(axis=0)

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
