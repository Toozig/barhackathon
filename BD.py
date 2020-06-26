import argparse
import matplotlib.pyplot as plt
import pandas as pd

import numpy as np

DT = 0.1
BD_FACTOR = 6
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
        """
        A function that initial coordination for the initial state
        :param c: the binding state
        :return: the coordination as np.array
        """
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
        """
        A function that gets protein vector and calculate the distances between every two proteins
        :param protein_vec: a list of all the proteins
        :return: distance matrix as numpy matrix
        """
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
            if self.__vec_in_rec(temp, LOWER_BOUND, LOWER_BOUND, 0, UPPER_BOUND, UPPER_BOUND, 0):
                cur = temp
                result[i] = cur
                i += 1
                j = 0
            else:
                for k in range(len(temp)):
                    temp[k][0] -= addition[k]
                    temp[k][1] -= addition[k]
                j += 1
        return result, self.__calculate_bars(result)

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
                    # x = dist_met[i,j]
                else:
                    x = 0
                calc = 1 if distance[x_pro][y_pro] else -1
                new_row.append(x * calc)
            new_mat.append(new_row)
        forces = np.array(new_mat).sum(axis=0)
        return forces


    def __calculate_bars(self, frames):
        """
        A function that gets the BD results and calculate the high of the fret bars
        :param frames: matrix of proteins coordinate
        :return: list which every index represent the light levels
        """
        errors = np.asarray([0 for i in range(len(frames))])
        dd = []
        for frame in range(len(frames)):
            d = []
            for i in range(len(frames[frame])):
                for j in range(i, len(frames[frame])):
                    if i == j:
                        continue
                    dist = np.linalg.norm(np.asarray(frames[frame][i]) - np.asarray(frames[frame][j]))
                    d.append(dist)
                    if dist < 1:
                        errors[frame] += 1
            dd.append(np.mean(d))
        print(np.mean(dd))
        return errors
