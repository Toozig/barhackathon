from BD import *
from generate_data import *

EMPTY_STATE = 'empty state'

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
        filtered = df.loc[self.__protein_vec, self.__protein_vec]
        filtered[EMPTY_STATE] = 1
        protein_dict = {}
        for i in range(filtered.shape[0] + 1):
            if i == filtered.shape[0]:
                protein_dict[i] = EMPTY_STATE
            else:
                protein_dict[i] = filtered.columns[i]
        return filtered, protein_dict


    def __get_p_accept_metropolis(self, dE, p_forward, p_backward):
        '''    return the probability to accept the metropolis criteria
            for the specified conditions
            dE - change in energy from current to proposed configuration
            kT - the factor of Boltzmann constant (kB) and the temperature (T)
            p_forward - probability to propose a move from current to proposed configuration
            p_backward - probability to propose a move from proposed to current configuration    '''
        p = np.exp(-dE / self.__kt) * p_backward / p_forward
        return min(p, 1.0)

    def __loss(self, bars):
        """
        calculate the loss score of the results from BD simulator
        :param frames: the result of BD
        :return: loss score
        """
        errors = bars
        err = np.mean((np.asarray(self.__compare) - errors))
        s = sum(self.__compare)
        return 1 if s == 0 else abs(err / s)


    def __E(self, c):
        """
        E function that used BD algorithm and calculate the loss of the result
        :param c: the initial state
        :return: the loss score
        """
        bd = BD(10, self.__kt,  self.__dt, c, 2, self.__protein_dict, self.__protein_vec, self.__df)
        res, bars = bd.BD_algorithm()
        return self.__loss(bars)

    def __get_neighbours(self, c):
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
        """
        A function that get tuple and returns it as plat list
        :param tup: the tuple to plat
        :return: list
        """
        if type(tup[0]) is int:
            seen = list(tup)
        else:
            seen = [x for y in tup for x in y]
        return seen

    def get_iter(self):
        return self.__m

    def mcmc(self):
        result = [0 for i in range(len(self.__space))]
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
            coin = np.random.binomial(1, p)
            random_state = chosen if coin else random_state
            conf = self.__space
            if type(chosen) is not tuple:
                chosen_ind = 0
            else:
                chosen_ind = conf[1:].index(chosen)
            random_ind = chosen_ind if coin else random_ind
        best_configuration = self.__space[np.argmax(np.asarray(result))]
        bd = BD(10, self.__kt,  self.__dt, best_configuration, 2, self.__protein_dict, self.__protein_vec, self.__df)
        return result, best_configuration, bd.BD_algorithm()




