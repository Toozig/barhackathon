import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
HIGHEST_NUMBER_OF_MOL = 500

LOWEST_NUM_OF_MOL = 50


def generate_histogram(amount_of_frames:int) -> list:
    """
    This function generate random lightening histogram
    :param amount_of_frames: the length of the wanted histogram
    :return: histogram
    """
    # choose random number of molecules
    num_mol = np.random.randint(LOWEST_NUM_OF_MOL, HIGHEST_NUMBER_OF_MOL, 1)[0]
    # calculates the max number of connections
    max_con = np.floor(num_mol / 2)
    # 1st value of the result
    cur_val = np.random.randint(0, max_con , 1)[0]
    result = [cur_val]
    i = 1
    # each iteration a factor is generated from normal dist with a
    # "cooling effect"
    while len(result) < amount_of_frames:
        sigma = max_con * ((amount_of_frames - len(result) + 1) / (amount_of_frames + i))
        factor =  round(np.random.normal(0, sigma, 1)[0])
        val = cur_val + factor
        i +=1
        if  0 >  val or val > max_con :
            continue
        cur_val = val
        result.append(cur_val)

#     return result
#
# res = []
# fig = go.Figure()
# for i in range(5):
#     res = generate_histogram(50)
#     fig.add_trace(go.Scatter(y=generate_histogram(50), x=list(range(50)), mode='lines', name='lines', orientation='h'))
#
# fig.show()