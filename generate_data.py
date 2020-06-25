import numpy as np
# import plotly.express as px
# import plotly.graph_objects as go
import pandas as pd
# from ipywidgets import widgets

HIGHEST_NUMBER_OF_MOL = 10
LOWEST_NUM_OF_MOL = 3
f_amount = 10
s_amount = 3


# frame_amount = widgets.IntSlider(
#     value=10,
#     min=1,
#     max=1000,
#     description='Number of frames in each sample:',
#     disabled=False)
#
# samples_amount = widgets.IntSlider(
#     value=3,
#     min=1,
#     max=110,
#     description='Number of samples:',
#     disabled=False)
#
# g = go.FigureWidget(data=[],
#                     layout=go.Layout(
#                         title=dict(
#                             text='NYC FlightDatabase'
#                         ),
#                         barmode='overlay'
#                     ))
#
#
# container = widgets.HBox(children=[samples_amount, frame_amount])


def generate_histogram(amount_of_frames:int) -> (list, int):
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
    # "cooling effect" such that each iteration the sigma is lower.
    # The factor is added to the current
    while len(result) < amount_of_frames:
        sigma = max_con * ((amount_of_frames - len(result) + 1) / (amount_of_frames + i))
        factor =  round(np.random.normal(0, sigma, 1)[0])
        val = cur_val + factor
        i +=1
        if 0 > val or val > max_con :
            continue
        cur_val = val
        result.append(cur_val)

    return result, num_mol


# def generate_scatter(change):
#     res = []
#     fig = go.Figure()
#     if frame_amount.value != f_amount:
#         for i in range(frame_amount.value):
#             res, num_mol = generate_histogram(samples_amount.value)
#             res.append(fig.add_trace(go.Scatter(y=res, x=list(range(samples_amount.value)), mode='lines', name='%i molecules' % num_mol, )))
#         f_amount = frame_amount.value
#         with g.batch_update():
#             g.data = res
#
#     elif samples_amount.value < len(g.data):
#         for i in range(samples_amount.value - len(g.data)):
#             res, num_mol = generate_histogram(frame_amount.value)
#             res.append(fig.add_trace(go.Scatter(y=res, x=list(range(frame_amount.value)), mode='lines', name='%i molecules' %num_mol,)))
#
#         with g.batch_update():
#             g.data = list(g.data) + res
#
#     else:
#         with g.batch_update():
#             g.data = g.data[:samples_amount.value]




# frame_amount.observe(generate_scatter, names="value")
# samples_amount.observe(generate_histogram, names='value')
# t = widgets.VBox([container])
# y = 45
