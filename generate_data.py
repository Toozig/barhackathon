#https://docs.google.com/document/d/19dhv9HE3doeHZAL-hx5jIrQxsNzB1egbFNLy7T8QWZk/edit#

import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
from ipywidgets import widgets

HIGHEST_NUMBER_OF_MOL = 500
LOWEST_NUM_OF_MOL = 50


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

    return go.Scatter(y=result, x=list(range(amount_of_frames)), name="%i molecules" %num_mol, visible=False)

fig = go.Figure()

for i in range(10):
    fig.add_trace(generate_histogram(1000))

for i in range(3):
    fig.data[i].visible = True

steps = []
for i in range(len(fig.data)):
    sample = dict(
        method="update",
        args=[{"visible": [False] * len(fig.data)},
              {"title": "Synthetic Data cooling effect (" + str(i) +" Samples)",
               }],  # layout attribute

    )
    for j in range(i):
        sample["args"][0]["visible"][j] = True  # Toggle i'th trace to "visible"
        steps.append(sample)


sliders = [dict(
    active=1,
    currentvalue={"prefix": "Number of samples: "},
    steps=steps

)]

fig.update_layout(
    sliders=sliders
)
fig.update_layout(title="Synthetic Data cooling effect",
                  xaxis_title="Time",
                  yaxis_title="# connections")

fig.write_html("S_data.html")
# fig.show()
