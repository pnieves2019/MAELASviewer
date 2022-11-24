# -*- coding: utf-8 -*-

# Run this app with `python3 maelasviewer.py` and
# visit http://127.0.0.1:8050/ in your web browser.

import os
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objects as go
from plotly.subplots import make_subplots
#import numpy as np
from dash.dependencies import Input, Output


app = dash.Dash(__name__)

server = app.server

app.config.suppress_callback_exceptions = True

app.layout = html.Div(children=[
    html.Div([html.Img(src=app.get_asset_url('logo_maelasviewer.png'))]),
    html.H1(children='MAELASviewer: Online visualization of magnetostriction'),
    html.Hr(),
    html.H5('Authors: P. Nieves, S. Arapan, A.P. Kądzielawa and D. Legut'),
    html.H6('Flagship Material Design via Exascale Computing at the IT4Innovations National Supercomputing Center at VŠB - Technical University of Ostrava, Czech Republic'),
    html.Hr(),

    html.H3("Introduction"),
    dcc.Markdown('''A magnetostrictive material is one which changes in size due to a change of state of magnetization. The main magnetostriction effects are the Joule effect (length change induced by a linear magnetic field),  Villari  effect  (magnetization  change  due  to  a  mechanical  stress  applied), Wiedemann effect (twisting of a magnetostrictive cylinder when helical magnetic field is applied to the material) and Matteucci effect (induced helical magnetization by a torsion). Presently, MAELASviewer contains interactive applets to simulate the field-induced magnetostriction effects (Joule effect and Wiedemann effect).'''),


    html.Hr(),
    html.H3('Select a magnetostriction effect:'),
    dcc.Dropdown(
        id='effect',
        options=[
            {'label': 'Joule Effect', 'value': 'Joule'},
            {'label': 'Wiedemann Effect', 'value': 'Wied'},
        ],
        placeholder="Select a magnetostriction effect"
    ),
    html.Div(id='effect_output'),








    html.Hr(),
    html.H3("Bibliography"),
    html.H6(" [1] Nieves, P.; Arapan, S.; Kądzielawa, A.P.; Legut, D. MAELASviewer: An Online Tool to Visualize Magnetostriction. Sensors 2020, 20, 6436. "),
    dcc.Markdown('''[https://www.mdpi.com/1424-8220/20/22/6436/pdf](https://www.mdpi.com/1424-8220/20/22/6436/pdf)'''),
    html.H6(" [2] P. Nieves, S. Arapan, S.H. Zhang, A.P. Kądzielawa, R.F. Zhang and D. Legut, MAELAS: MAgneto-ELAStic properties calculation via computational high-throughput approach, Comput. Phys. Commun. 264, 107964 (2021). "),
    dcc.Markdown('''[https://doi.org/10.1016/j.cpc.2021.107964](https://doi.org/10.1016/j.cpc.2021.107964)'''),
    html.H6(" [3] P. Nieves, S. Arapan, S.H. Zhang, A.P. Kądzielawa, R.F. Zhang and D. Legut, MAELAS 2.0: A new version of a computer program for the calculation of magneto-elastic properties, Comput. Phys. Commun. 271, 108197 (2022). "),
    dcc.Markdown('''[https://doi.org/10.1016/j.cpc.2021.108197](https://doi.org/10.1016/j.cpc.2021.108197)'''),
    html.H3("Source files"),
    dcc.Markdown('''[https://github.com/pnieves2019/MAELASviewer](https://github.com/pnieves2019/MAELASviewer)'''),
    html.Hr(),


])



if __name__ == '__main__':
        app.run_server(debug=True)
