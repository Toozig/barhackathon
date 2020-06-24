import argparse
import numpy as np
import pandas as pd
from BD import *


class MCHC:
    def __init__(self, n, kt, m, protein_vec, dt, compare):
        """

        :param n: proteins number
        :param kt: kT
        :param m: iteration number
        :param protein_vec: list of proteins
        :param dt: time interval
        :param compare: real result for comparision
        """
        self.__n = n
        self.__kt = kt
        self.__m = m
        self.__protein_vec = protein_vec
        self.__dt = dt
        self.__compare = compare

    def is_valid(self, c):
        ''' Return True if c is a valid 2-D coordinate on an n x n grid    with 0-based indices '''
        if len(c) != 2:
            return False
        return c[0] >= 0 and c[1] >= 0 and c[0] < self.__n and c[1] < self.__n

    def __get_p_accept_metropolis(self, dE, p_forward, p_backward):
        '''    return the probability to accept the metropolis criteria
            for the specified conditions
            dE - change in energy from current to proposed configuration
            kT - the factor of Boltzmann constant (kB) and the temperature (T)
            p_forward - probability to propose a move from current to proposed configuration
            p_backward - probability to propose a move from proposed to current configuration    '''
        p = np.exp(-dE / self.__kt) * p_backward / p_forward
        return min(p, 1.0)

    def __loss(self, frames):
        """
        calculate the loss score of the results from BD simulator
        :param frames: the result of BD
        :return: loss score
        """
        errors = np.asarray([0 for i in range(len(frames))])
        for frame in frames:
            for i in range(len(frame)):
                for j in range(len(i , frame)):
                    if np.linalg.norm(frame[i] - frame[j]) < 0.1:   #todo
                        errors[frame] += 1
        err = (np.asarray(self.__compare) - errors) / np.asarray(self.__compare)
        return np.mean(err)



    def __E(self, c):
        """
        E function that used BD algorithm and calculate the loss of the result
        :param c: the initial state
        :return: the loss score
        """
        bd = BD(10000, self.__kt,  self.__dt, c, 2)
        res = bd.BD_algorithm()
        return self.__loss(res)

    def __get_neighbours(self, c): #todo
        ''' get up/down/left/right neighbours on an n x n grid with 0-based indices'''
        orig = pd.read_excel("matrix.xlsx")
        ret_value = []
        if c[0] > 0:
            ret_value.append((c[0] - 1, c[1]))
        if c[0] < self.__n - 1:
            ret_value.append((c[0] + 1, c[1]))
        if c[1] > 0:
            ret_value.append((c[0], c[1] - 1))
        if c[1] < self.__n - 1:
            ret_value.append((c[0], c[1] + 1))
        return ret_value


    def __create_configuration_space(self): #todo
        filtered = 0
        return filtered


    def mcmc(self):
        configuration = self.__create_configuration_space()
        result = [0 for i in range(len(configuration))]
        # matrix = [(i, j) for j in range(self.__n) for i in range(self.__n)]
        random_state = configuration[np.random.randint(len(configuration))]
        for iter in range(self.__m):
            result[random_state] += 1
            neighbours = self.__get_neighbours(random_state)
            chosen = neighbours[np.random.randint(len(neighbours))]
            dE = self.__E(chosen) - self.__E(random_state)
            p = self.__get_p_accept_metropolis(dE, self.__kt, 1 / len(neighbours), 1 / len(self.__get_neighbours(chosen)))
            random_state = chosen if np.random.binomial(1, p) else random_state
        return result


# def get_cmdline_parser():
#     parser = argparse.ArgumentParser(
#         description='Run MCMC on a discrete n x n configuration space.')
#     parser.add_argument('n', type=int, default=5, nargs='?',
#                         help='number of rows/columns of grid')
#     parser.add_argument('m', type=int, default=1000, nargs='?',
#                         help='number of iterations of MCMC optimization')
#     parser.add_argument('kT', type=float, default=1.0, nargs='?',
#                         help='kT - the denominator for the metropolis criterion'
#                              ' (Boltzmann constant times temperature)')
#     return parser
if __name__ == '__main__':
    print("hi")
    # parser = get_cmdline_parser().parse_args()
    # print(mcmc(parser))


