import argparse
import numpy as np
import pandas as pd
from BD import *
from itertools import chain, combinations
from generate_data import *

PROTEIN_DICT = {0: "CD28", 1: "CD80", 2: "CD86" , 3: "CTLA-4", 4: "PD-1", 5: "PD-L1", 6: "PD-L2", 7: "IL-2RA",
                8: "IL-12R", 9: "IL-2", 10: "IL-12", 11: "empty state"}
CONVERT_DICT = {"CD28":0, "CD80":1, "CD86":2, "CTLA-4":3,  "PD-1":4, "PD-L1":5, "PD-L2":6, "IL-2RA":7, "IL-12R":8,
                "IL-2":9, "IL-12":10, "empty state":11}
class MCMC:

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
        self.__df, self.__protein_dict = self.__read_excel()
        self.__space = self.__create_configuration_space()


    def __read_excel(self):
        df = pd.read_csv("newMatrix.csv", index_col=0)
        # colnames = df.columns
        # col_ind = [np.argwhere([i is in self.__protein_vec for i in colnames])]
        filtered = df.loc[self.__protein_vec, self.__protein_vec]
        # filtered = filtered.iloc[col_ind]
        filtered['empty state'] = 1
        protein_dict = {}
        for i in range(filtered.shape[0] + 1):
            if i == filtered.shape[0]:
                protein_dict[i] = 'empty state'
            else:
                protein_dict[i] = filtered.columns[i]
        return filtered, protein_dict


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
        for frame in range(len(frames)):
            for i in range(len(frames[frame])):
                for j in range(i, len(frames[frame])):
                    if i == j:
                        continue
                    dist = np.linalg.norm(np.asarray(frames[frame][i]) - np.asarray(frames[frame][j]))
                    if  dist < 5:   #todo
                        errors[frame] += 1
        err = np.mean((np.asarray(self.__compare) - errors))
        s = sum(self.__compare)
        return 1 if s == 0 else err / s



    def __E(self, c):
        """
        E function that used BD algorithm and calculate the loss of the result
        :param c: the initial state
        :return: the loss score
        """
        bd = BD(10, self.__kt,  self.__dt, c, 2, self.__protein_dict, self.__protein_vec, self.__df)
        res = bd.BD_algorithm()
        return self.__loss(res)

    def __get_neighbours(self, c):
        """

        :param c:
        :return:
        """
        return self.__space


    def __create_configuration_space(self):
        """
        the function build the configuration space from the given input
        :return: the configuration space as list
        """
        configurations = []
        i = 0
        while i < self.__n:
            if i == 0:
                configurations.append([(self.__n)])
            elif i == 2:
                temp = []
                for protein in range(self.__n):
                    for sec_protein in range(protein, self.__n):
                        if self.__df[self.__protein_dict[protein]][self.__protein_dict[sec_protein]] == 1:
                            temp.append((protein, sec_protein))
                configurations.append(temp)
            elif i % 2 == 0:
                temp2 = []
                for tup in configurations[-1]:
                    seen = self.__plat_the_tup(tup)
                    for pair in configurations[1]:
                        if pair[0] not in seen and pair[1] not in seen:
                            temp2.append((tup, pair))
                configurations.append(temp2)
            i +=1
        result = (lambda l: [item for sublist in l for item in sublist])(configurations)
        return result

    def __plat_the_tup(self, tup):
        if type(tup[0]) is int:
            seen = list(tup)
        else:
            seen = [x for y in tup for x in y]
        return seen

    def mcmc(self):
        result = [0 for i in range(len(self.__space))]
        # matrix = [(i, j) for j in range(self.__n) for i in range(self.__n)]
        random_ind = np.random.randint(len(self.__space))
        random_state = self.__space[random_ind]
        for iter in range(self.__m):
            result[random_ind] += 1
            neighbours = self.__get_neighbours(random_ind)
            chosen = neighbours[np.random.randint(len(neighbours))]
            c_E =  self.__E(chosen)
            r_E = self.__E(random_state)
            dE =  c_E - r_E
            p = self.__get_p_accept_metropolis(dE, 1 / len(neighbours), 1 / len(self.__get_neighbours(chosen)))
            random_state = chosen if np.random.binomial(1, p) else random_state
        return result



def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return list(chain.from_iterable(list(combinations(s, r)) for r in range(len(s)+1)))






# def findsubsets(s, n):
#     base = list(itertools.combinations(s, n))
#     level = []
#     for pair in base:
#         temp = [i for i in base if i[0] != pair[0] and i[0] != pair[1] and i[1] != pair[0] and i[1] != pair[1]]
#         for j in temp:
#             level.append((pair, j))
#     return level

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
    # print(findsubsets([1, 2, 3, 4], 2))
    # parser = get_cmdline_parser().parse_args()
    # print(mcmc(parser))
    # l = [1,2,3,4]
    # pairs = list(combinations(l, 2))
    # psets = powerset(pairs)
    # valid = []
    # for pset in psets:
    #     already_seen = []
    #     add = True
    #     for pair in pset:
    #         if [x for x in pair if x in already_seen]:
    #             add = False
    #         already_seen.extend(pair)
    #
    #     if add:
    #         valid.append(pset)
    # print(valid)
    data, n = generate_histogram(10)
    proteins = [PROTEIN_DICT[i] for i in range(n)]
    sim = MCMC(n, 1.0, 1000, proteins, 0.1, data)
    res = sim.mcmc()
    print(res)



