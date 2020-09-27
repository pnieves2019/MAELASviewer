# -*- coding: utf-8 -*-

# Run this app with `python3 maelasviewer.py` and
# visit http://127.0.0.1:8050/ in your web browser.

import os
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from dash.dependencies import Input, Output

app = dash.Dash(__name__)

server = app.server

app.layout = html.Div(children=[
    html.Div([html.Img(src=app.get_asset_url('logo_maelasviewer.png'))]),
    html.H1(children='MAELASviewer: Online visualization of magnetostriction'),
    html.Hr(),
    html.H5('Authors: P. Nieves, S. Arapan, A.P. Kądzielawa and D. Legut'),
    html.H6('Flagship Material Design via Exascale Computing at the IT4Innovations National Supercomputing Center at VŠB - Technical University of Ostrava'),
    html.Hr(),

    html.H3("Introduction"),
    html.H6("A magnetostrictive material is one which changes in size due to a change of state of magnetization. The main magnetostriction effects are the Joule effect (length change induced by a linear magnetic field),  Villari  effect  (magnetization  change  due  to  a  mechanical  stress  applied), Wiedemann effect (twisting of a magnetostrictive cylinder when helical magnetic field is applied to the material) and Matteucci effect (induced helical magnetization by a torsion). Presently, MAELASviewer contains interactive applets to simulate the field-induced magnetostriction effects (Joule effect and Wiedemann effect)."),
    
    
    
    html.Hr(),
    html.H3("The Joule effect"),
    html.H6("This interactive applet shows the magnetostriction due to Joule effect for some crystal systems. You can visualize the relative length change (\u0394l/lo=[l-lo]/lo) of the material along an arbitrary direction (β) as a function of the external magnetic field (H) and magnetostrictive coefficients (λ). The magnitude of the external magnetic field is assumed to be strong enough to saturate the magnetization (α) along the magnetic field (α||H). The length lo corresponds to the size of the magnetic material in a demagnetized state along direction \u03B2=(sinθ*cosφ, sinθ*sinφ, cosθ), where θ and φ are the polar and azimuthal angles, respectively. In the simulation the length lo of the material in any direction β is represented by a sphere with radius 1 (lo=1), so that it does not depend on the geometry of the material. However, note that the simulations can also be understood as the real shape deformation of a magnetic spherical nanoparticle with diameter larger than the domain wall width (in order to allow the formation of magnetic domains)."),
    html.Div([html.Img(src=app.get_asset_url('diagram_online.png'))]),
    html.Hr(),
    html.H3("Available systems"),
    html.H6("1. Single crystal: Cubic I (space group numbers  207-230)"),
    html.H6("2. Polycrystal: Cubic I (space group numbers  207-230)"),
    html.H6("3. Single crystal: Hexagonal I (space group numbers  177-194)"),
    html.H6("4. Single crystal: Trigonal I (space group numbers  149-167)"),
    html.H6("5. Single crystal: Tetragonal I (space group numbers  89-142)"),
    html.H6("6. Single crystal: Orthorhombic (space group numbers  16-74)"),
    html.Hr(),
    html.H3("1. Single crystal: Cubic I (space group numbers  207-230)"),
    html.H4("1.1 Theory"),
    html.H6("The relative length change for cubic (I) systems is given by:"),
    html.Div([html.Img(src=app.get_asset_url('eq_cub.png'))]),
    html.H6("where αi and βi (i=x,y,z) are the direction of magnetization (parallel to the external magentic field H) and the measured length direction, respectively."),

    html.H5("Scaling Factor:"),
    html.H6("Magentostriction is a small effect that is hard to visualize. To facilitate its visualization in the simulation, we multiply the right hand side of this equation by a scaling factor parameter (s) that can be modified by the user"),
    html.Div([html.Img(src=app.get_asset_url('eq_cub_s1.png'))]),
    html.H6("Note that this scaling preserve the ratio between the magnetostrictive coefficients. Obviously, the case with s=1 corresponds to the real situation. The total length l is"),
    html.Div([html.Img(src=app.get_asset_url('eq_cub_s2.png'))]),
    html.H6("where we took into account that lo=1. Similar procedure is applied to the other supported crystal systems."),

    html.H4("1.2 Parameters of the simulation"),
    html.H6("(Press Enter after changing any input to update the figures)"),
    html.Hr(),
    html.H6("Direction of the external magnetic field:"),
    html.Div(['H',html.Sub('x')," = ",
              dcc.Input(id='fieldx', value=1.0, type='number', debounce=True, step=0.1)]),
    html.Div(['H',html.Sub('y')," = ",
              dcc.Input(id='fieldy', value=0.0, type='number', debounce=True, step=0.1)]),
    html.Div(['H',html.Sub('z')," = ",
              dcc.Input(id='fieldz', value=0.0, type='number', debounce=True, step=0.1)]),
    html.H6("Magnetostrictive coefficients:"),
    html.Div(['\u03BB',html.Sup("\u03B1")," = ",
              dcc.Input(id='L0', value=0.0, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sub("001")," = ",
              dcc.Input(id='L1', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sub("111")," = ",
              dcc.Input(id='L2', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(["Scaling factor = ",
              dcc.Input(id='scale', value=1.0, type='number', debounce=True)]),

    html.Div(id='my-output'),
    html.Hr(),
    html.H4("1.3 Simulation"),
    html.H6("The distance between a point on the surface and the origin (0,0,0) describes the length l along direction \u03B2=(sinθ*cosφ, sinθ*sinφ, cosθ), where θ and φ are the polar and azimuthal angles, respectively. The color of the surface corresponds to the relative length change multiplied by the scaling factor along direction \u03B2."),
    dcc.Graph(id='cub_3D'),
    dcc.Graph(id='cub_2D'),
    html.H4("1.4 Magnetostriction of some cubic crystals"),

    html.Table([
        html.Tr([html.Td('Material'), html.Td('Space group'),html.Td('Temperature (K)'),html.Td(html.Div(['    \u03BB',html.Sup("\u03B1"),'    '])), html.Td(html.Div(['\u03BB',html.Sub("001")])),html.Td(html.Div(['\u03BB',html.Sub("111")]))]),
        html.Tr([html.Td('BCC Fe'), html.Td('229'),html.Td('0'),html.Td('    -   '),html.Td('0.000026'), html.Td('-0.00003')]),
        html.Tr([html.Td('BCC Fe'), html.Td('229'),html.Td('300'),html.Td('   -   '),html.Td('0.000021'), html.Td('-0.000021')]),
        html.Tr([html.Td('FCC Ni'), html.Td('225'),html.Td('0'),html.Td('   -   '),html.Td('-0.00006'), html.Td('-0.000035')]),
        html.Tr([html.Td('FCC Ni'), html.Td('225'),html.Td('300'),html.Td('   -   '),html.Td('-0.000046'), html.Td('-0.000024')]),
        html.Tr([html.Td(html.Div(['C15 SmFe',html.Sub("2")])), html.Td('227'),html.Td('0'),html.Td('   -   '),html.Td('0.00003'), html.Td('-0.0041')]),
        html.Tr([html.Td(html.Div(['C15 DyFe',html.Sub("2")])), html.Td('227'),html.Td('0'),html.Td('   -   '),html.Td('-0.00007'), html.Td('0.003')]),
        html.Tr([html.Td(html.Div(['C15 TbCo',html.Sub("2")])), html.Td('227'),html.Td('0'),html.Td('   -   '),html.Td('-0.0012'), html.Td('0.0045')]),
        html.Tr([html.Td(html.Div(['C15 ErCo',html.Sub("2")])), html.Td('227'),html.Td('0'),html.Td('   -   '),html.Td('-0.001'), html.Td('-0.0025')]),
    ]),
    html.Hr(),
####################### Cubic I polycrystal

    html.H3("2. Polycrystal: Cubic I (space group numbers  207-230)"),
    html.H4("2.1 Theory"),
    html.H6("The theory of magnetostriction for polycrystalline materials is more complex. A widely used approximation is to assume that the stress distribution is uniform through the material. In this case the relative change in length may be put into the form:"),
    html.Div([html.Img(src=app.get_asset_url('eq_cub_poly.png'))]),
    html.H6("where"),
    html.Div([html.Img(src=app.get_asset_url('eq_cub_poly_lmb_s.png'))]),

    html.H4("2.2 Parameters of the simulation"),
    html.H6("(Press Enter after changing any input to update the figures)"),
    html.Hr(),
    html.H6("Direction of the external magnetic field:"),
    html.Div(['H',html.Sub('x')," = ",
              dcc.Input(id='pfieldx', value=1.0, type='number', debounce=True, step=0.1)]),
    html.Div(['H',html.Sub('y')," = ",
              dcc.Input(id='pfieldy', value=0.0, type='number', debounce=True, step=0.1)]),
    html.Div(['H',html.Sub('z')," = ",
              dcc.Input(id='pfieldz', value=0.0, type='number', debounce=True, step=0.1)]),
    html.H6("Magnetostrictive coefficients:"),
    html.Div(['\u03BB',html.Sub("S")," = ",
              dcc.Input(id='LS', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(["Scaling factor = ",
              dcc.Input(id='pscale', value=1.0, type='number', debounce=True)]),

    html.Div(id='my-output-p'),
    html.Hr(),
    html.H4("2.3 Simulation"),
    html.H6("The distance between a point on the surface and the origin (0,0,0) describes the length l along direction \u03B2=(sinθ*cosφ, sinθ*sinφ, cosθ), where θ and φ are the polar and azimuthal angles, respectively. The color of the surface corresponds to the relative length change multiplied by the scaling factor along direction \u03B2."),

    dcc.Graph(id='graph-poly'),
    dcc.Graph(id='graph-poly2D'),
    html.H4("2.4 Magnetostriction of some polycrystal materials (cubic I)"),
    html.Table([
        html.Tr([html.Td('Material'), html.Td('Space group'),html.Td('Temperature (K)'),html.Td(html.Div(['\u03BB',html.Sub("S")]))]),
        html.Tr([html.Td('BCC Fe'), html.Td('229'),html.Td('300'),html.Td('-0.000007')]),
        html.Tr([html.Td('FCC Ni'), html.Td('225'),html.Td('300'),html.Td('-0.000034')]),
        html.Tr([html.Td(html.Div(['C15 TbFe',html.Sub("2")])), html.Td('227'),html.Td('300'),html.Td('0.001753')]),
    ]),
    html.Hr(),

######## Hex I

    html.H3("3. Single crystal: Hexagonal I (space group numbers  177-194)"),
    html.H4("3.1 Theory"),
    html.H6("The relative length change for hexagonal (I) systems is given by:"),
    html.Div([html.Img(src=app.get_asset_url('eq_hex.png'))]),
    html.H6("where αi and βi (i=x,y,z) are the direction of magnetization (parallel to the external magentic field H) and the measured length direction, respectively. Magentostriction is a small effect that is hard to visualize. To facilitate its visualization in the simulation, we multiply the right hand side of this equation by a scaling factor parameter which can be modified by the user."),

    html.H4("3.2 Parameters of the simulation"),
    html.H6("(Press Enter after changing any input to update the figures)"),
    html.Hr(),
    html.H6("Direction of the external magnetic field:"),
    html.Div(['H',html.Sub('x')," = ",
              dcc.Input(id='hfieldx', value=1.0, type='number', debounce=True, step=0.1)]),
    html.Div(['H',html.Sub('y')," = ",
              dcc.Input(id='hfieldy', value=0.0, type='number', debounce=True, step=0.1)]),
    html.Div(['H',html.Sub('z')," = ",
              dcc.Input(id='hfieldz', value=0.0, type='number', debounce=True, step=0.1)]),
    html.H6("Magnetostrictive coefficients:"),
    html.Div(['\u03BB',html.Sup("\u03B11,0")," = ",
              dcc.Input(id='hL01', value=0.0, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sup("\u03B12,0")," = ",
              dcc.Input(id='hL02', value=0.0, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sup("\u03B11,2")," = ",
              dcc.Input(id='hLa1', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sup("\u03B12,2")," = ",
              dcc.Input(id='hLa2', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sup("\u03B3,2")," = ",
              dcc.Input(id='hLg', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sup("\u03B5,2")," = ",
              dcc.Input(id='hLe', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(["Scaling factor = ",
              dcc.Input(id='hscale', value=1.0, type='number', debounce=True)]),

    html.Div(id='my-output-h'),
    html.Hr(),
    html.H4("3.3 Simulation"),
    html.H6("The distance between a point on the surface and the origin (0,0,0) describes the length l along direction \u03B2=(sinθ*cosφ, sinθ*sinφ, cosθ), where θ and φ are the polar and azimuthal angles, respectively. The color of the surface corresponds to the relative length change multiplied by the scaling factor along direction \u03B2."),
    dcc.Graph(id='hex_3D'),
    dcc.Graph(id='hex_2D'),
    html.H4("3.4 Magnetostriction of some hexagonal crystals"),

    html.Table([
        html.Tr([html.Td('Material'), html.Td('Space group'),html.Td('Temperature (K)'),html.Td(html.Div(['    \u03BB',html.Sup('\u03B11,0'),'    '])), html.Td(html.Div(['    \u03BB',html.Sup('\u03B12,0'),'    '])), html.Td(html.Div(['\u03BB',html.Sup('\u03B11,2')])), html.Td(html.Div(['\u03BB',html.Sup('\u03B12,2')])), html.Td(html.Div(['\u03BB',html.Sup('\u03B3,2')])),html.Td(html.Div(['\u03BB',html.Sup('\u03B5,2')]))]),
        html.Tr([html.Td('HCP Co'), html.Td('194'),html.Td('0'),html.Td('    -   '),html.Td('    -   '), html.Td('0.000095'), html.Td('-0.000126'),html.Td('0.000057'),html.Td('-0.000286')]),
        html.Tr([html.Td('HCP Gd'), html.Td('194'),html.Td('0'),html.Td('    -   '),html.Td('    -   '), html.Td('0.00014'), html.Td('-0.00013'),html.Td('0.00011'),html.Td('0.00002')]),
        html.Tr([html.Td('HCP Tb'), html.Td('194'),html.Td('0'),html.Td('    -   '),html.Td('    -   '), html.Td('-0.0026'), html.Td('0.009'),html.Td('0.0087'),html.Td('0.015')]),
    ]),
    html.Hr(),


######## Trig I

    html.H3("4. Single crystal: Trigonal I (space group numbers  149-167)"),
    html.H4("4.1 Theory"),
    html.H6("The relative length change for trigonal (I) systems is given by:"),
    html.Div([html.Img(src=app.get_asset_url('eq_trig.png'))]),
    html.H6("where αi and βi (i=x,y,z) are the direction of magnetization (parallel to the external magentic field H) and the measured length direction, respectively. Magentostriction is a small effect that is hard to visualize. To facilitate its visualization in the simulation, we multiply the right hand side of this equation by a scaling factor parameter which can be modified by the user."),

    html.H4("4.2 Parameters of the simulation"),
    html.H6("(Press Enter after changing any input to update the figures)"),
    html.Hr(),
    html.H6("Direction of the external magnetic field:"),
    html.Div(['H',html.Sub('x')," = ",
              dcc.Input(id='trfieldx', value=1.0, type='number', debounce=True, step=0.1)]),
    html.Div(['H',html.Sub('y')," = ",
              dcc.Input(id='trfieldy', value=0.0, type='number', debounce=True, step=0.1)]),
    html.Div(['H',html.Sub('z')," = ",
              dcc.Input(id='trfieldz', value=0.0, type='number', debounce=True, step=0.1)]),
    html.H6("Magnetostrictive coefficients:"),
    html.Div(['\u03BB',html.Sup("\u03B11,0")," = ",
              dcc.Input(id='trL01', value=0.0, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sup("\u03B12,0")," = ",
              dcc.Input(id='trL02', value=0.0, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sup("\u03B11,2")," = ",
              dcc.Input(id='trLa1', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sup("\u03B12,2")," = ",
              dcc.Input(id='trLa2', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sup("\u03B3,1")," = ",
              dcc.Input(id='trLg1', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sup("\u03B3,2")," = ",
              dcc.Input(id='trLg2', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sub("12")," = ",
              dcc.Input(id='trL12', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sub("21")," = ",
              dcc.Input(id='trL21', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(["Scaling factor = ",
              dcc.Input(id='trscale', value=1.0, type='number', debounce=True)]),


    html.Div(id='my-output-tr'),
    html.Hr(),
    html.H4("4.3 Simulation"),
    html.H6("The distance between a point on the surface and the origin (0,0,0) describes the length l along direction \u03B2=(sinθ*cosφ, sinθ*sinφ, cosθ), where θ and φ are the polar and azimuthal angles, respectively. The color of the surface corresponds to the relative length change multiplied by the scaling factor along direction \u03B2."),
    dcc.Graph(id='tri_3D'),
    dcc.Graph(id='tri_2D'),

    html.Hr(),

######## Tet I

    html.H3("5. Single crystal: Tetragonal I (space group numbers  89-142)"),
    html.H4("5.1 Theory"),
    html.H6("The relative length change for tetragoanl (I) systems is given by:"),
    html.Div([html.Img(src=app.get_asset_url('eq_tet.png'))]),
    html.H6("where αi and βi (i=x,y,z) are the direction of magnetization (parallel to the external magentic field H) and the measured length direction, respectively. Magentostriction is a small effect that is hard to visualize. To facilitate its visualization in the simulation, we multiply the right hand side of this equation by a scaling factor parameter which can be modified by the user."),

    html.H4("5.2 Parameters of the simulation"),
    html.H6("(Press Enter after changing any input to update the figures)"),
    html.Hr(),
    html.H6("Direction of the external magnetic field:"),
    html.Div(['H',html.Sub('x')," = ",
              dcc.Input(id='tefieldx', value=1.0, type='number', debounce=True, step=0.1)]),
    html.Div(['H',html.Sub('y')," = ",
              dcc.Input(id='tefieldy', value=0.0, type='number', debounce=True, step=0.1)]),
    html.Div(['H',html.Sub('z')," = ",
              dcc.Input(id='tefieldz', value=0.0, type='number', debounce=True, step=0.1)]),
    html.H6("Magnetostrictive coefficients:"),
    html.Div(['\u03BB',html.Sup("\u03B11,0")," = ",
              dcc.Input(id='teL01', value=0.0, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sup("\u03B12,0")," = ",
              dcc.Input(id='teL02', value=0.0, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sup("\u03B11,2")," = ",
              dcc.Input(id='teLa1', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sup("\u03B12,2")," = ",
              dcc.Input(id='teLa2', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sup("\u03B3,2")," = ",
              dcc.Input(id='teLg', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sup("\u03B4,2")," = ",
              dcc.Input(id='teLd', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sup("\u03B5,2")," = ",
              dcc.Input(id='teLe', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(["Scaling factor = ",
              dcc.Input(id='tescale', value=1.0, type='number', debounce=True)]),

    html.Div(id='my-output-te'),
    html.Hr(),
    html.H4("5.3 Simulation"),
    html.H6("The distance between a point on the surface and the origin (0,0,0) describes the length l along direction \u03B2=(sinθ*cosφ, sinθ*sinφ, cosθ), where θ and φ are the polar and azimuthal angles, respectively. The color of the surface corresponds to the relative length change multiplied by the scaling factor along direction \u03B2."),
    dcc.Graph(id='tet_3D'),
    dcc.Graph(id='tet_2D'),

    html.Hr(),


######## Orth

    html.H3("6. Single crystal: Orthorhombic (space group numbers  16-74)"),
    html.H4("6.1 Theory"),
    html.H6("This interactive applet shows the magnetostriction due to Joule effect for some crystal systems."),
    html.Div([html.Img(src=app.get_asset_url('eq_ort.png'))]),
    html.H6("where αi and βi (i=x,y,z) are the direction of magnetization (parallel to the external magentic field H) and the measured length direction, respectively. Magentostriction is a small effect that is hard to visualize. To facilitate its visualization in the simulation, we multiply the right hand side of this equation by a scaling factor parameter which can be modified by the user."),

    html.H4("6.2 Parameters of the simulation"),
    html.H6("(Press Enter after changing any input to update the figures)"),
    html.Hr(),
    html.H6("Direction of the external magnetic field:"),
    html.Div(['H',html.Sub('x')," = ",
              dcc.Input(id='ofieldx', value=1.0, type='number', debounce=True, step=0.1)]),
    html.Div(['H',html.Sub('y')," = ",
              dcc.Input(id='ofieldy', value=0.0, type='number', debounce=True, step=0.1)]),
    html.Div(['H',html.Sub('z')," = ",
              dcc.Input(id='ofieldz', value=0.0, type='number', debounce=True, step=0.1)]),
    html.H6("Magnetostrictive coefficients:"),
    html.Div(['\u03BB',html.Sup("\u03B11,0")," = ",
              dcc.Input(id='oL01', value=0.0, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sup("\u03B12,0")," = ",
              dcc.Input(id='oL02', value=0.0, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sup("\u03B13,0")," = ",
              dcc.Input(id='oL03', value=0.0, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sub("1")," = ",
              dcc.Input(id='oL1', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sub("2")," = ",
              dcc.Input(id='oL2', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sub("3")," = ",
              dcc.Input(id='oL3', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sub("4")," = ",
              dcc.Input(id='oL4', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sub("5")," = ",
              dcc.Input(id='oL5', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sub("6")," = ",
              dcc.Input(id='oL6', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sub("7")," = ",
              dcc.Input(id='oL7', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sub("8")," = ",
              dcc.Input(id='oL8', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(['\u03BB',html.Sub("9")," = ",
              dcc.Input(id='oL9', value=0.000001, type='number', debounce=True, step=0.000001)]),
    html.Div(["Scaling factor = ",
              dcc.Input(id='oscale', value=1.0, type='number', debounce=True)]),

    html.Div(id='my-output-o'),
    html.Hr(),
    html.H4("6.3 Simulation"),
    html.H6("The distance between a point on the surface and the origin (0,0,0) describes the length l along direction \u03B2=(sinθ*cosφ, sinθ*sinφ, cosθ), where θ and φ are the polar and azimuthal angles, respectively. The color of the surface corresponds to the relative length change multiplied by the scaling factor along direction \u03B2."),
    dcc.Graph(id='ort_3D'),
    dcc.Graph(id='ort_2D'),

    html.Hr(),
    
    html.H3("The Wiedemann effect"),
    html.H4("Theory"),
    
    
    html.H6("This interactive applet simulates the Wiedemann effect for a isotropic cylindrical rod (the twist of a rod induced by helical mangetic field). The helical magnetic field to twist a magnetic rod is achieved by applying a magnetic field along the rod height axis (longitudinal field H‖), and an electric current (I) through the rod which induces a circular magnetic field due to Ampère’s circuital law (perpendicular field H⊥)."),
    html.Div([html.Img(src=app.get_asset_url('torsion_rod.png'))]),
    html.H6("For an isotropic magnetic cylindrical rod aligned to the z-axis with height L and radius R, the twisted angle ϕ induced by a helical magnetic field"),
    html.Div([html.Img(src=app.get_asset_url('eq_twist.png'))]),
    html.H6("where \u03BBs is the isotropic magnetostrictive coefficient and A = 4πR^2 is the area of the cross section of the rod. The helical field-induced torque can be calculated as"),
    html.Div([html.Img(src=app.get_asset_url('eq_torque.png'))]),
    html.H6("where Y is the Young’s modulus and σ is the Poisson’s ratio. In the 3D visualization of the rod, it is also plotted the perpendicular and longitudinal magnetic fields in the exterior of the rod at z = L/2. The perpendicular field is calculated applying the Biot-Savart law for a finite wire. Note that the magnetic field generated by the magnetization of the rod is not plotted here."),

    html.H4("Parameters of the simulation"),
    html.H6("(Press Enter after changing any input to update the figures)"),
    html.Hr(),
    html.H6("Geometry of the magnetic rod:"),
    html.Div(['L(m) = ',
              dcc.Input(id='LL', value=0.0005, type='number', debounce=True, step=0.000000001)]),
    html.Div(['R(m) = ',
              dcc.Input(id='RR', value=0.00001, type='number', debounce=True, step=0.000000001)]),
    html.H6("Longitudinal external magnetic field:"),
    html.Div(['H',html.Sub('||'),"(A/m) = ",
              dcc.Input(id='hlong', value=0.000001, type='number', debounce=True, step=0.000000001)]),
    html.H6("Isotropic magnetostrictive coefficient:"),
    html.Div(['\u03BB',html.Sub('S')," = ",
              dcc.Input(id='lmbs', value=0.000001, type='number', debounce=True, step=0.0000001)]),
    html.H6("Range of applied electric current:"),
    html.Div(['I',html.Sub('min'),"(A) = ",
              dcc.Input(id='imin', value=-0.00000001, type='number', debounce=True, step=0.000000001)]),
    html.Div(['I',html.Sub('max'),"(A) = ",
              dcc.Input(id='imax', value=0.00000001, type='number', debounce=True, step=0.000000001)]),
    html.H6("Isotropic elastic properties:"),
    html.Div(['Young modulus Y(GPa) = ',
              dcc.Input(id='y', value=100.0, type='number', debounce=True, step=1.0)]),
    html.Div(['Poisson ratio σ = ',
              dcc.Input(id='sigma', value=0.33, type='number', debounce=True, step=0.001)]),
    html.Div(id='my-output-rod'),
    html.Hr(),
    
    html.H4("Simulation"),
    
    dcc.Graph(id='rod'),
    dcc.Graph(id='twist_angle'),
    dcc.Graph(id='torque'),
    
    
    
    
    
    
    
    
    html.Hr(),
    html.H3("Bibliography"),
    html.H6(" [1] P. Nieves, S. Arapan, A.P. Kądzielawa and D. Legut, MAELASviewer: an online tool to visualize magnetostriction, 2020, Submitted to Sensors"),
    html.H6(" [2] P. Nieves, S. Arapan, S.H. Zhang, A.P. Kądzielawa, R.F. Zhang and D. Legut, MAELAS: MAgneto-ELAStic properties calculation via computational high-throughput approach, 2020, arXiv:2009.01638"),
    html.H3("Source files"),
    html.H6("https://github.com/pnieves2019/MAELASwiewer"),
    html.Hr(),


])

##############Cubic I -3D

@app.callback(
    Output('cub_3D', 'figure'),
    [Input(component_id='fieldx', component_property='value'),
     Input(component_id='fieldy', component_property='value'),
     Input(component_id='fieldz', component_property='value'),
     Input(component_id='L0', component_property='value'),
     Input(component_id='L1', component_property='value'),
     Input(component_id='L2', component_property='value'),
     Input(component_id='scale', component_property='value'),
    ]
)


def update_fig(hx,hy,hz,lmb0,lmb1,lmb2,s):

    d=1.0
    h=np.sqrt(hx*hx+hy*hy+hz*hz)
    ax=hx/h
    ay=hy/h
    az=hz/h
    u, v = np.mgrid[0:np.pi:100j, 0:2*np.pi:100j]
    bx = np.sin(u)*np.cos(v)
    by = np.sin(u)*np.sin(v)
    bz = np.cos(u)
    f = lmb0+lmb1*1.5*(ax*ax*bx*bx+ay*ay*by*by+az*az*bz*bz-(1/3))+3*lmb2*(ax*ay*bx*by+ay*az*by*bz+ax*az*bx*bz)
    x = d*(1.0+s*f)*bx
    y = d*(1.0+s*f)*by
    z = d*(1.0+s*f)*bz
    dl_l=(np.sqrt(x**2+y**2+z**2)-d)/d

    fig = make_subplots(rows=1, cols=2,
                    specs=[[{'is_3d': True}, {'is_3d': True}]],
                    subplot_titles=['   Magnetic field (H), magnetization (\u03B1) and unit cell lattice vectors (a,b,c)           ', '        Color corresponds to (\u0394l/lo)*scaling_factor along direction \u03B2'],)

    fig.add_trace(go.Cone(x=[1], y=[1], z=[1], u=[ax], v=[ay], w=[az],name="H",colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(255,0,0)']]),1, 1)
    fig.add_trace(go.Cone(x=[1], y=[1], z=[0], u=[ax], v=[ay], w=[az],name="\u03B1",colorscale=[[0, 'rgb(255,127,14)'], [1, 'rgb(255,127,16)']]),1, 1)
    fig.add_trace(go.Cone(x=[0.4], y=[0], z=[0], u=[2], v=[0], w=[0],name="a",colorscale=[[0, 'rgb(255,0,0)'], [1, 'rgb(240,0,0)']]),1, 1)
    fig.add_trace(go.Cone(x=[0], y=[0.4], z=[0], u=[0], v=[2], w=[0],name="b",colorscale=[[0, 'rgb(0,255,0)'], [1, 'rgb(0,240,0)']]),1, 1)
    fig.add_trace(go.Cone(x=[0], y=[0], z=[0.4], u=[0], v=[0], w=[2],name="c",colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(0,0,240)']]),1, 1)

    fig.update_traces(hoverinfo="name", showscale=False)
    fig.add_trace(go.Surface(x=x, y=y, z=z, surfacecolor=dl_l, name="l"), 1, 2)

    fig.update_layout(transition_duration=500)

    fig.update_yaxes(automargin=True)
    fig.update_xaxes(automargin=True)

    return fig

######################################## Cub 2D

@app.callback(
    Output('cub_2D', 'figure'),
    [Input(component_id='fieldx', component_property='value'),
     Input(component_id='fieldy', component_property='value'),
     Input(component_id='fieldz', component_property='value'),
     Input(component_id='L0', component_property='value'),
     Input(component_id='L1', component_property='value'),
     Input(component_id='L2', component_property='value'),
     Input(component_id='scale', component_property='value'),
    ]
)


def update_figc2d(hx,hy,hz,lmb0,lmb1,lmb2,s):

    d=1.0
    h=np.sqrt(hx*hx+hy*hy+hz*hz)
    ax=hx/h
    ay=hy/h
    az=hz/h
    vxy = np.mgrid[0:2*np.pi:200j]
    bxxy = np.cos(vxy)
    byxy = np.sin(vxy)
    bzxy = 0.0
    fxy = lmb0+lmb1*1.5*(ax*ax*bxxy*bxxy+ay*ay*byxy*byxy+az*az*bzxy*bzxy-(1/3))+3*lmb2*(ax*ay*bxxy*byxy+ay*az*byxy*bzxy+ax*az*bxxy*bzxy)
    xxy = d*(1.0+s*fxy)*bxxy
    yxy = d*(1.0+s*fxy)*byxy
    zxy = d*(1.0+s*fxy)*bzxy

    vxz = np.mgrid[0:2*np.pi:200j]
    bxxz = np.sin(vxz)
    byxz = 0.0
    bzxz = np.cos(vxz)
    fxz = lmb0+lmb1*1.5*(ax*ax*bxxz*bxxz+ay*ay*byxz*byxz+az*az*bzxz*bzxz-(1/3))+3*lmb2*(ax*ay*bxxz*byxz+ay*az*byxz*bzxz+ax*az*bxxz*bzxz)
    xxz = d*(1.0+s*fxz)*bxxz
    yxz = d*(1.0+s*fxz)*byxz
    zxz = d*(1.0+s*fxz)*bzxz

    vyz = np.mgrid[0:2*np.pi:200j]
    bxyz = 0.0
    byyz = np.sin(vyz)
    bzyz = np.cos(vyz)
    fyz = lmb0+lmb1*1.5*(ax*ax*bxyz*bxyz+ay*ay*byyz*byyz+az*az*bzyz*bzyz-(1/3))+3*lmb2*(ax*ay*bxyz*byyz+ay*az*byyz*bzyz+ax*az*bxyz*bzyz)
    xyz = d*(1.0+s*fyz)*bxyz
    yyz = d*(1.0+s*fyz)*byyz
    zyz = d*(1.0+s*fyz)*bzyz

    figc2d = make_subplots(rows=1, cols=3,
                    specs=[[{'type': 'xy'}, {'type': 'xy'}, {'type': 'xy'}]],
                    subplot_titles=['Plot XY plane (z=0)', 'Plot XZ plane (y=0)', 'Plot YZ plane (x=0)'],)

    figc2d.add_trace(go.Scatter(x=xxy, y=yxy, mode='lines',name='(lx,ly)'), 1, 1)
    figc2d.add_trace(go.Scatter(x=bxxy, y=byxy, mode='lines',name='(lox,loy)'), 1, 1)

    figc2d.add_trace(go.Scatter(x=xxz, y=zxz, mode='lines',name='(lx,lz)'), 1, 2)
    figc2d.add_trace(go.Scatter(x=bxxz, y=bzxz, mode='lines',name='(lox,loz)'), 1, 2)

    figc2d.add_trace(go.Scatter(x=yyz, y=zyz, mode='lines',name='(ly,lz)'), 1, 3)
    figc2d.add_trace(go.Scatter(x=byyz, y=bzyz, mode='lines',name='(loy,loz)'), 1, 3)
    figc2d.update_layout(xaxis_title='X-axis')
    figc2d.update_layout(yaxis_title='Y-axis')
    figc2d.update_layout(transition_duration=500)

    figc2d.update_yaxes(automargin=True)
    figc2d.update_xaxes(automargin=True)

    return figc2d



##################### Cubic polycrystal 3D

@app.callback(
    Output('graph-poly', 'figure'),
    [Input(component_id='pfieldx', component_property='value'),
     Input(component_id='pfieldy', component_property='value'),
     Input(component_id='pfieldz', component_property='value'),
     Input(component_id='LS', component_property='value'),
     Input(component_id='pscale', component_property='value'),
    ]
)



def update_figure(hx,hy,hz,lmbs,s):

    d=1.0
    h=np.sqrt(hx*hx+hy*hy+hz*hz)
    ax=hx/h
    ay=hy/h
    az=hz/h
    u, v = np.mgrid[0:np.pi:100j, 0:2*np.pi:100j]
    bx = np.sin(u)*np.cos(v)
    by = np.sin(u)*np.sin(v)
    bz = np.cos(u)
    f = lmbs*1.5*((ax*bx+ay*by+az*bz)**2-(1/3))
    x = d*(1.0+s*f)*bx
    y = d*(1.0+s*f)*by
    z = d*(1.0+s*f)*bz

    figp = make_subplots(rows=1, cols=2,
                    specs=[[{'is_3d': True}, {'is_3d': True}]],
                    subplot_titles=['   Magnetic field (H), magnetization (\u03B1) and unit cell lattice vectors (a,b,c)           ', '        Color corresponds to (\u0394l/lo)*scaling_factor along direction \u03B2'],)

    figp.add_trace(go.Cone(x=[1], y=[1], z=[1], u=[ax], v=[ay], w=[az],name="H",colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(255,0,0)']]),1, 1)
    figp.add_trace(go.Cone(x=[1], y=[1], z=[0], u=[ax], v=[ay], w=[az],name="\u03B1",colorscale=[[0, 'rgb(255,127,14)'], [1, 'rgb(255,127,16)']]),1, 1)
    figp.add_trace(go.Cone(x=[0.4], y=[0], z=[0], u=[2], v=[0], w=[0],name="a",colorscale=[[0, 'rgb(255,0,0)'], [1, 'rgb(240,0,0)']]),1, 1)
    figp.add_trace(go.Cone(x=[0], y=[0.4], z=[0], u=[0], v=[2], w=[0],name="b",colorscale=[[0, 'rgb(0,255,0)'], [1, 'rgb(0,240,0)']]),1, 1)
    figp.add_trace(go.Cone(x=[0], y=[0], z=[0.4], u=[0], v=[0], w=[2],name="c",colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(0,0,240)']]),1, 1)

    figp.update_traces(hoverinfo="name", showscale=False)
    figp.add_trace(go.Surface(x=x, y=y, z=z, surfacecolor=(np.sqrt(x**2+y**2+z**2)-d)/d, name="\u03B2"), 1, 2)
    figp.update_layout(transition_duration=500)
    figp.update_yaxes(automargin=True)

    return figp

######################################## Cub poly 2D

@app.callback(
    Output('graph-poly2D', 'figure'),
    [Input(component_id='pfieldx', component_property='value'),
     Input(component_id='pfieldy', component_property='value'),
     Input(component_id='pfieldz', component_property='value'),
     Input(component_id='LS', component_property='value'),
     Input(component_id='pscale', component_property='value'),
    ]
)

def update_figcp2d(hx,hy,hz,lmbs,s):

    d=1.0
    h=np.sqrt(hx*hx+hy*hy+hz*hz)
    ax=hx/h
    ay=hy/h
    az=hz/h
    vxy = np.mgrid[0:2*np.pi:200j]
    bxxy = np.cos(vxy)
    byxy = np.sin(vxy)
    bzxy = 0.0
    fxy = lmbs*1.5*((ax*bxxy+ay*byxy+az*bzxy)**2-(1/3))
    xxy = d*(1.0+s*fxy)*bxxy
    yxy = d*(1.0+s*fxy)*byxy
    zxy = d*(1.0+s*fxy)*bzxy

    vxz = np.mgrid[0:2*np.pi:200j]
    bxxz = np.sin(vxz)
    byxz = 0.0
    bzxz = np.cos(vxz)
    fxz = lmbs*1.5*((ax*bxxz+ay*byxz+az*bzxz)**2-(1/3))
    xxz = d*(1.0+s*fxz)*bxxz
    yxz = d*(1.0+s*fxz)*byxz
    zxz = d*(1.0+s*fxz)*bzxz

    vyz = np.mgrid[0:2*np.pi:200j]
    bxyz = 0.0
    byyz = np.sin(vyz)
    bzyz = np.cos(vyz)
    fyz = lmbs*1.5*((ax*bxyz+ay*byyz+az*bzyz)**2-(1/3))
    xyz = d*(1.0+s*fyz)*bxyz
    yyz = d*(1.0+s*fyz)*byyz
    zyz = d*(1.0+s*fyz)*bzyz

    figcp2d = make_subplots(rows=1, cols=3,
                    specs=[[{'type': 'xy'}, {'type': 'xy'}, {'type': 'xy'}]],
                    subplot_titles=['Plot XY plane (z=0)', 'Plot XZ plane (y=0)', 'Plot YZ plane (x=0)'],)

    figcp2d.add_trace(go.Scatter(x=xxy, y=yxy, mode='lines',name='(lx,ly)'), 1, 1)
    figcp2d.add_trace(go.Scatter(x=bxxy, y=byxy, mode='lines',name='(lox,loy)'), 1, 1)

    figcp2d.add_trace(go.Scatter(x=xxz, y=zxz, mode='lines',name='(lx,lz)'), 1, 2)
    figcp2d.add_trace(go.Scatter(x=bxxz, y=bzxz, mode='lines',name='(lox,loz)'), 1, 2)

    figcp2d.add_trace(go.Scatter(x=yyz, y=zyz, mode='lines',name='(ly,lz)'), 1, 3)
    figcp2d.add_trace(go.Scatter(x=byyz, y=bzyz, mode='lines',name='(loy,loz)'), 1, 3)
    figcp2d.update_layout(xaxis_title='X-axis')
    figcp2d.update_layout(yaxis_title='Y-axis')
    figcp2d.update_layout(transition_duration=500)

    figcp2d.update_yaxes(automargin=True)
    figcp2d.update_xaxes(automargin=True)

    return figcp2d



############## hex-3D

@app.callback(
    Output('hex_3D', 'figure'),
    [Input(component_id='hfieldx', component_property='value'),
     Input(component_id='hfieldy', component_property='value'),
     Input(component_id='hfieldz', component_property='value'),
     Input(component_id='hL01', component_property='value'),
     Input(component_id='hL02', component_property='value'),
     Input(component_id='hLa1', component_property='value'),
     Input(component_id='hLa2', component_property='value'),
     Input(component_id='hLg', component_property='value'),
     Input(component_id='hLe', component_property='value'),
     Input(component_id='hscale', component_property='value'),
    ]
)


def update_hfig(hx,hy,hz,lmb01,lmb02,lmba1,lmba2,lmbg,lmbe,s):

    d=1.0
    h=np.sqrt(hx*hx+hy*hy+hz*hz)
    ax=hx/h
    ay=hy/h
    az=hz/h
    u, v = np.mgrid[0:np.pi:100j, 0:2*np.pi:100j]
    bx = np.sin(u)*np.cos(v)
    by = np.sin(u)*np.sin(v)
    bz = np.cos(u)
    f = lmb01*(bx**2+by**2)+lmb02*(bz**2)+lmba1*((az**2)-(1/3))*(bx**2+by**2)+lmba2*((az**2)-(1/3))*(bz**2)+lmbg*(0.5*(ax**2-ay**2)*(bx**2-by**2)+2*ax*ay*bx*by)+2*lmbe*(ax*az*bx*bz+az*ay*bz*by)
    x = d*(1.0+s*f)*bx
    y = d*(1.0+s*f)*by
    z = d*(1.0+s*f)*bz
    dl_l=(np.sqrt(x**2+y**2+z**2)-d)/d

    hfig = make_subplots(rows=1, cols=2,
                    specs=[[{'is_3d': True}, {'is_3d': True}]],
                    subplot_titles=['   Magnetic field (H), magnetization (\u03B1) and unit cell lattice vectors (a,b,c)           ', '        Color corresponds to (\u0394l/lo)*scaling_factor along direction \u03B2'],)

    hfig.add_trace(go.Cone(x=[1], y=[1], z=[1], u=[ax], v=[ay], w=[az],name="H",colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(255,0,0)']]),1, 1)
    hfig.add_trace(go.Cone(x=[1], y=[1], z=[0], u=[ax], v=[ay], w=[az],name="\u03B1",colorscale=[[0, 'rgb(255,127,14)'], [1, 'rgb(255,127,16)']]),1, 1)
    hfig.add_trace(go.Cone(x=[0.4], y=[0], z=[0], u=[2], v=[0], w=[0],name="a",colorscale=[[0, 'rgb(255,0,0)'], [1, 'rgb(240,0,0)']]),1, 1)
    hfig.add_trace(go.Cone(x=[-0.2], y=[0.34641016], z=[0], u=[-1], v=[1.7320508], w=[0],name="b",colorscale=[[0, 'rgb(0,255,0)'], [1, 'rgb(0,240,0)']]),1, 1)
    hfig.add_trace(go.Cone(x=[0], y=[0], z=[0.4], u=[0], v=[0], w=[2],name="c",colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(0,0,240)']]),1, 1)

    hfig.update_traces(hoverinfo="name", showscale=False)
    hfig.add_trace(go.Surface(x=x, y=y, z=z, surfacecolor=dl_l, name="l"), 1, 2)

    hfig.update_layout(transition_duration=500)

    hfig.update_yaxes(automargin=True)
    hfig.update_xaxes(automargin=True)

    return hfig

######################################## Hex 2D

@app.callback(
    Output('hex_2D', 'figure'),
    [Input(component_id='hfieldx', component_property='value'),
     Input(component_id='hfieldy', component_property='value'),
     Input(component_id='hfieldz', component_property='value'),
     Input(component_id='hL01', component_property='value'),
     Input(component_id='hL02', component_property='value'),
     Input(component_id='hLa1', component_property='value'),
     Input(component_id='hLa2', component_property='value'),
     Input(component_id='hLg', component_property='value'),
     Input(component_id='hLe', component_property='value'),
     Input(component_id='hscale', component_property='value'),
    ]
)


def update_figh2d(hx,hy,hz,lmb01,lmb02,lmba1,lmba2,lmbg,lmbe,s):

    d=1.0
    h=np.sqrt(hx*hx+hy*hy+hz*hz)
    ax=hx/h
    ay=hy/h
    az=hz/h
    vxy = np.mgrid[0:2*np.pi:200j]
    bxxy = np.cos(vxy)
    byxy = np.sin(vxy)
    bzxy = 0.0
    fxy = lmb01*(bxxy**2+byxy**2)+lmb02*(bzxy**2)+lmba1*((az**2)-(1/3))*(bxxy**2+byxy**2)+lmba2*((az**2)-(1/3))*(bzxy**2)+lmbg*(0.5*(ax**2-ay**2)*(bxxy**2-byxy**2)+2*ax*ay*bxxy*byxy)+2*lmbe*(ax*az*bxxy*bzxy+az*ay*bzxy*byxy)
    xxy = d*(1.0+s*fxy)*bxxy
    yxy = d*(1.0+s*fxy)*byxy
    zxy = d*(1.0+s*fxy)*bzxy

    vxz = np.mgrid[0:2*np.pi:200j]
    bxxz = np.sin(vxz)
    byxz = 0.0
    bzxz = np.cos(vxz)
    fxz = lmb01*(bxxz**2+byxz**2)+lmb02*(bzxz**2)+lmba1*((az**2)-(1/3))*(bxxz**2+byxz**2)+lmba2*((az**2)-(1/3))*(bzxz**2)+lmbg*(0.5*(ax**2-ay**2)*(bxxz**2-byxz**2)+2*ax*ay*bxxz*byxz)+2*lmbe*(ax*az*bxxz*bzxz+az*ay*bzxz*byxz)
    xxz = d*(1.0+s*fxz)*bxxz
    yxz = d*(1.0+s*fxz)*byxz
    zxz = d*(1.0+s*fxz)*bzxz

    vyz = np.mgrid[0:2*np.pi:200j]
    bxyz = 0.0
    byyz = np.sin(vyz)
    bzyz = np.cos(vyz)
    fyz = lmb01*(bxyz**2+byyz**2)+lmb02*(bzyz**2)+lmba1*((az**2)-(1/3))*(bxyz**2+byyz**2)+lmba2*((az**2)-(1/3))*(bzyz**2)+lmbg*(0.5*(ax**2-ay**2)*(bxyz**2-byyz**2)+2*ax*ay*bxyz*byyz)+2*lmbe*(ax*az*bxyz*bzyz+az*ay*bzyz*byyz)
    xyz = d*(1.0+s*fyz)*bxyz
    yyz = d*(1.0+s*fyz)*byyz
    zyz = d*(1.0+s*fyz)*bzyz

    figh2d = make_subplots(rows=1, cols=3,
                    specs=[[{'type': 'xy'}, {'type': 'xy'}, {'type': 'xy'}]],
                    subplot_titles=['Plot XY plane (z=0)', 'Plot XZ plane (y=0)', 'Plot YZ plane (x=0)'],)

    figh2d.add_trace(go.Scatter(x=xxy, y=yxy, mode='lines',name='(lx,ly)'), 1, 1)
    figh2d.add_trace(go.Scatter(x=bxxy, y=byxy, mode='lines',name='(lox,loy)'), 1, 1)

    figh2d.add_trace(go.Scatter(x=xxz, y=zxz, mode='lines',name='(lx,lz)'), 1, 2)
    figh2d.add_trace(go.Scatter(x=bxxz, y=bzxz, mode='lines',name='(lox,loz)'), 1, 2)

    figh2d.add_trace(go.Scatter(x=yyz, y=zyz, mode='lines',name='(ly,lz)'), 1, 3)
    figh2d.add_trace(go.Scatter(x=byyz, y=bzyz, mode='lines',name='(loy,loz)'), 1, 3)
    figh2d.update_layout(xaxis_title='X-axis')
    figh2d.update_layout(yaxis_title='Y-axis')
    figh2d.update_layout(transition_duration=500)

    figh2d.update_yaxes(automargin=True)
    figh2d.update_xaxes(automargin=True)

    return figh2d





############## tri-3D

@app.callback(
    Output('tri_3D', 'figure'),
    [Input(component_id='trfieldx', component_property='value'),
     Input(component_id='trfieldy', component_property='value'),
     Input(component_id='trfieldz', component_property='value'),
     Input(component_id='trL01', component_property='value'),
     Input(component_id='trL02', component_property='value'),
     Input(component_id='trLa1', component_property='value'),
     Input(component_id='trLa2', component_property='value'),
     Input(component_id='trLg1', component_property='value'),
     Input(component_id='trLg2', component_property='value'),
     Input(component_id='trL12', component_property='value'),
     Input(component_id='trL21', component_property='value'),
     Input(component_id='trscale', component_property='value'),
    ]
)


def update_trfig(hx,hy,hz,lmb01,lmb02,lmba1,lmba2,lmbg1,lmbg2,lmb12,lmb21,s):

    d=1.0
    h=np.sqrt(hx*hx+hy*hy+hz*hz)
    ax=hx/h
    ay=hy/h
    az=hz/h
    u, v = np.mgrid[0:np.pi:100j, 0:2*np.pi:100j]
    bx = np.sin(u)*np.cos(v)
    by = np.sin(u)*np.sin(v)
    bz = np.cos(u)
    f = lmb01*(bx**2+by**2)+lmb02*(bz**2)+lmba1*((az**2)-(1/3))*(bx**2+by**2)+lmba2*((az**2)-(1/3))*(bz**2)+lmbg1*(0.5*(ax**2-ay**2)*(bx**2-by**2)+ax*ay*bx*by)+lmbg2*(ax*az*bx*bz+ay*az*by*bz)+lmb12*(0.5*ay*az*(bx**2-by**2)+ax*az*bx*by)+lmb21*(0.5*(ax**2-ay**2)*by*bz+ax*ay*bx*bz)
    x = d*(1.0+s*f)*bx
    y = d*(1.0+s*f)*by
    z = d*(1.0+s*f)*bz
    dl_l=(np.sqrt(x**2+y**2+z**2)-d)/d

    trfig = make_subplots(rows=1, cols=2,
                    specs=[[{'is_3d': True}, {'is_3d': True}]],
                    subplot_titles=['   Magnetic field (H), magnetization (\u03B1) and unit cell lattice vectors (a,b,c)           ', '        Color corresponds to (\u0394l/lo)*scaling_factor along direction \u03B2'],)

    trfig.add_trace(go.Cone(x=[1], y=[1], z=[1], u=[ax], v=[ay], w=[az],name="H",colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(255,0,0)']]),1, 1)
    trfig.add_trace(go.Cone(x=[1], y=[1], z=[0], u=[ax], v=[ay], w=[az],name="\u03B1",colorscale=[[0, 'rgb(255,127,14)'], [1, 'rgb(255,127,16)']]),1, 1)
    trfig.add_trace(go.Cone(x=[0.4], y=[0], z=[0], u=[2], v=[0], w=[0],name="a",colorscale=[[0, 'rgb(255,0,0)'], [1, 'rgb(240,0,0)']]),1, 1)
    trfig.add_trace(go.Cone(x=[-0.2], y=[0.34641016], z=[0], u=[-1], v=[1.7320508], w=[0],name="b",colorscale=[[0, 'rgb(0,255,0)'], [1, 'rgb(0,240,0)']]),1, 1)
    trfig.add_trace(go.Cone(x=[0], y=[0], z=[0.4], u=[0], v=[0], w=[2],name="c",colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(0,0,240)']]),1, 1)

    trfig.update_traces(hoverinfo="name", showscale=False)
    trfig.add_trace(go.Surface(x=x, y=y, z=z, surfacecolor=dl_l, name="l"), 1, 2)

    trfig.update_layout(transition_duration=500)

    trfig.update_yaxes(automargin=True)
    trfig.update_xaxes(automargin=True)

    return trfig

######################################## Tri 2D

@app.callback(
    Output('tri_2D', 'figure'),
    [Input(component_id='trfieldx', component_property='value'),
     Input(component_id='trfieldy', component_property='value'),
     Input(component_id='trfieldz', component_property='value'),
     Input(component_id='trL01', component_property='value'),
     Input(component_id='trL02', component_property='value'),
     Input(component_id='trLa1', component_property='value'),
     Input(component_id='trLa2', component_property='value'),
     Input(component_id='trLg1', component_property='value'),
     Input(component_id='trLg2', component_property='value'),
     Input(component_id='trL12', component_property='value'),
     Input(component_id='trL21', component_property='value'),
     Input(component_id='trscale', component_property='value'),
    ]
)


def update_figtr2d(hx,hy,hz,lmb01,lmb02,lmba1,lmba2,lmbg1,lmbg2,lmb12,lmb21,s):

    d=1.0
    h=np.sqrt(hx*hx+hy*hy+hz*hz)
    ax=hx/h
    ay=hy/h
    az=hz/h
    vxy = np.mgrid[0:2*np.pi:200j]
    bxxy = np.cos(vxy)
    byxy = np.sin(vxy)
    bzxy = 0.0
    fxy = lmb01*(bxxy**2+byxy**2)+lmb02*(bzxy**2)+lmba1*((az**2)-(1/3))*(bxxy**2+byxy**2)+lmba2*((az**2)-(1/3))*(bzxy**2)+lmbg1*(0.5*(ax**2-ay**2)*(bxxy**2-byxy**2)+ax*ay*bxxy*byxy)+lmbg2*(ax*az*bxxy*bzxy+ay*az*byxy*bzxy)+lmb12*(0.5*ay*az*(bxxy**2-byxy**2)+ax*az*bxxy*byxy)+lmb21*(0.5*(ax**2-ay**2)*byxy*bzxy+ax*ay*bxxy*bzxy)
    xxy = d*(1.0+s*fxy)*bxxy
    yxy = d*(1.0+s*fxy)*byxy
    zxy = d*(1.0+s*fxy)*bzxy

    vxz = np.mgrid[0:2*np.pi:200j]
    bxxz = np.sin(vxz)
    byxz = 0.0
    bzxz = np.cos(vxz)
    fxz = lmb01*(bxxz**2+byxz**2)+lmb02*(bzxz**2)+lmba1*((az**2)-(1/3))*(bxxz**2+byxz**2)+lmba2*((az**2)-(1/3))*(bzxz**2)+lmbg1*(0.5*(ax**2-ay**2)*(bxxz**2-byxz**2)+ax*ay*bxxz*byxz)+lmbg2*(ax*az*bxxz*bzxz+ay*az*byxz*bzxz)+lmb12*(0.5*ay*az*(bxxz**2-byxz**2)+ax*az*bxxz*byxz)+lmb21*(0.5*(ax**2-ay**2)*byxz*bzxz+ax*ay*bxxz*bzxz)
    xxz = d*(1.0+s*fxz)*bxxz
    yxz = d*(1.0+s*fxz)*byxz
    zxz = d*(1.0+s*fxz)*bzxz

    vyz = np.mgrid[0:2*np.pi:200j]
    bxyz = 0.0
    byyz = np.sin(vyz)
    bzyz = np.cos(vyz)
    fyz = lmb01*(bxyz**2+byyz**2)+lmb02*(bzyz**2)+lmba1*((az**2)-(1/3))*(bxyz**2+byyz**2)+lmba2*((az**2)-(1/3))*(bzyz**2)+lmbg1*(0.5*(ax**2-ay**2)*(bxyz**2-byyz**2)+ax*ay*bxyz*byyz)+lmbg2*(ax*az*bxyz*bzyz+ay*az*byyz*bzyz)+lmb12*(0.5*ay*az*(bxyz**2-byyz**2)+ax*az*bxyz*byyz)+lmb21*(0.5*(ax**2-ay**2)*byyz*bzyz+ax*ay*bxyz*bzyz)
    xyz = d*(1.0+s*fyz)*bxyz
    yyz = d*(1.0+s*fyz)*byyz
    zyz = d*(1.0+s*fyz)*bzyz

    figtr2d = make_subplots(rows=1, cols=3,
                    specs=[[{'type': 'xy'}, {'type': 'xy'}, {'type': 'xy'}]],
                    subplot_titles=['Plot XY plane (z=0)', 'Plot XZ plane (y=0)', 'Plot YZ plane (x=0)'],)

    figtr2d.add_trace(go.Scatter(x=xxy, y=yxy, mode='lines',name='(lx,ly)'), 1, 1)
    figtr2d.add_trace(go.Scatter(x=bxxy, y=byxy, mode='lines',name='(lox,loy)'), 1, 1)

    figtr2d.add_trace(go.Scatter(x=xxz, y=zxz, mode='lines',name='(lx,lz)'), 1, 2)
    figtr2d.add_trace(go.Scatter(x=bxxz, y=bzxz, mode='lines',name='(lox,loz)'), 1, 2)

    figtr2d.add_trace(go.Scatter(x=yyz, y=zyz, mode='lines',name='(ly,lz)'), 1, 3)
    figtr2d.add_trace(go.Scatter(x=byyz, y=bzyz, mode='lines',name='(loy,loz)'), 1, 3)
    figtr2d.update_layout(xaxis_title='X-axis')
    figtr2d.update_layout(yaxis_title='Y-axis')
    figtr2d.update_layout(transition_duration=500)

    figtr2d.update_yaxes(automargin=True)
    figtr2d.update_xaxes(automargin=True)

    return figtr2d


############## tet-3D

@app.callback(
    Output('tet_3D', 'figure'),
    [Input(component_id='tefieldx', component_property='value'),
     Input(component_id='tefieldy', component_property='value'),
     Input(component_id='tefieldz', component_property='value'),
     Input(component_id='teL01', component_property='value'),
     Input(component_id='teL02', component_property='value'),
     Input(component_id='teLa1', component_property='value'),
     Input(component_id='teLa2', component_property='value'),
     Input(component_id='teLg', component_property='value'),
     Input(component_id='teLd', component_property='value'),
     Input(component_id='teLe', component_property='value'),
     Input(component_id='tescale', component_property='value'),
    ]
)


def update_tefig(hx,hy,hz,lmb01,lmb02,lmba1,lmba2,lmbg,lmbd,lmbe,s):

    d=1.0
    h=np.sqrt(hx*hx+hy*hy+hz*hz)
    ax=hx/h
    ay=hy/h
    az=hz/h
    u, v = np.mgrid[0:np.pi:100j, 0:2*np.pi:100j]
    bx = np.sin(u)*np.cos(v)
    by = np.sin(u)*np.sin(v)
    bz = np.cos(u)
    f = lmb01*(bx**2+by**2)+lmb02*(bz**2)+lmba1*((az**2)-(1/3))*(bx**2+by**2)+lmba2*((az**2)-(1/3))*(bz**2)+lmbg*0.5*(ax**2-ay**2)*(bx**2-by**2)+lmbd*2*ax*ay*bx*by+2*lmbe*(ax*az*bx*bz+az*ay*bz*by)
    x = d*(1.0+s*f)*bx
    y = d*(1.0+s*f)*by
    z = d*(1.0+s*f)*bz
    dl_l=(np.sqrt(x**2+y**2+z**2)-d)/d

    tefig = make_subplots(rows=1, cols=2,
                    specs=[[{'is_3d': True}, {'is_3d': True}]],
                    subplot_titles=['   Magnetic field (H), magnetization (\u03B1) and unit cell lattice vectors (a,b,c)           ', '        Color corresponds to (\u0394l/lo)*scaling_factor along direction \u03B2'],)

    tefig.add_trace(go.Cone(x=[1], y=[1], z=[1], u=[ax], v=[ay], w=[az],name="H",colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(255,0,0)']]),1, 1)
    tefig.add_trace(go.Cone(x=[1], y=[1], z=[0], u=[ax], v=[ay], w=[az],name="\u03B1",colorscale=[[0, 'rgb(255,127,14)'], [1, 'rgb(255,127,16)']]),1, 1)
    tefig.add_trace(go.Cone(x=[0.4], y=[0], z=[0], u=[2], v=[0], w=[0],name="a",colorscale=[[0, 'rgb(255,0,0)'], [1, 'rgb(240,0,0)']]),1, 1)
    tefig.add_trace(go.Cone(x=[0], y=[0.4], z=[0], u=[0], v=[2], w=[0],name="b",colorscale=[[0, 'rgb(0,255,0)'], [1, 'rgb(0,240,0)']]),1, 1)
    tefig.add_trace(go.Cone(x=[0], y=[0], z=[0.4], u=[0], v=[0], w=[2],name="c",colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(0,0,240)']]),1, 1)

    tefig.update_traces(hoverinfo="name", showscale=False)
    tefig.add_trace(go.Surface(x=x, y=y, z=z, surfacecolor=dl_l, name="l"), 1, 2)

    tefig.update_layout(transition_duration=500)

    tefig.update_yaxes(automargin=True)
    tefig.update_xaxes(automargin=True)

    return tefig

######################################## Tet 2D

@app.callback(
    Output('tet_2D', 'figure'),
    [Input(component_id='tefieldx', component_property='value'),
     Input(component_id='tefieldy', component_property='value'),
     Input(component_id='tefieldz', component_property='value'),
     Input(component_id='teL01', component_property='value'),
     Input(component_id='teL02', component_property='value'),
     Input(component_id='teLa1', component_property='value'),
     Input(component_id='teLa2', component_property='value'),
     Input(component_id='teLg', component_property='value'),
     Input(component_id='teLd', component_property='value'),
     Input(component_id='teLe', component_property='value'),
     Input(component_id='tescale', component_property='value'),
    ]
)

def update_figte2d(hx,hy,hz,lmb01,lmb02,lmba1,lmba2,lmbg,lmbd,lmbe,s):

    d=1.0
    h=np.sqrt(hx*hx+hy*hy+hz*hz)
    ax=hx/h
    ay=hy/h
    az=hz/h
    vxy = np.mgrid[0:2*np.pi:200j]
    bxxy = np.cos(vxy)
    byxy = np.sin(vxy)
    bzxy = 0.0
    fxy = lmb01*(bxxy**2+byxy**2)+lmb02*(bzxy**2)+lmba1*((az**2)-(1/3))*(bxxy**2+byxy**2)+lmba2*((az**2)-(1/3))*(bzxy**2)+lmbg*0.5*(ax**2-ay**2)*(bxxy**2-byxy**2)+lmbd*2*ax*ay*bxxy*byxy+2*lmbe*(ax*az*bxxy*bzxy+az*ay*bzxy*byxy)
    xxy = d*(1.0+s*fxy)*bxxy
    yxy = d*(1.0+s*fxy)*byxy
    zxy = d*(1.0+s*fxy)*bzxy

    vxz = np.mgrid[0:2*np.pi:200j]
    bxxz = np.sin(vxz)
    byxz = 0.0
    bzxz = np.cos(vxz)
    fxz = lmb01*(bxxz**2+byxz**2)+lmb02*(bzxz**2)+lmba1*((az**2)-(1/3))*(bxxz**2+byxz**2)+lmba2*((az**2)-(1/3))*(bzxz**2)+lmbg*0.5*(ax**2-ay**2)*(bxxz**2-byxz**2)+lmbd*2*ax*ay*bxxz*byxz+2*lmbe*(ax*az*bxxz*bzxz+az*ay*bzxz*byxz)
    xxz = d*(1.0+s*fxz)*bxxz
    yxz = d*(1.0+s*fxz)*byxz
    zxz = d*(1.0+s*fxz)*bzxz

    vyz = np.mgrid[0:2*np.pi:200j]
    bxyz = 0.0
    byyz = np.sin(vyz)
    bzyz = np.cos(vyz)
    fyz = lmb01*(bxyz**2+byyz**2)+lmb02*(bzyz**2)+lmba1*((az**2)-(1/3))*(bxyz**2+byyz**2)+lmba2*((az**2)-(1/3))*(bzyz**2)+lmbg*0.5*(ax**2-ay**2)*(bxyz**2-byyz**2)+lmbd*2*ax*ay*bxyz*byyz+2*lmbe*(ax*az*bxyz*bzyz+az*ay*bzyz*byyz)
    xyz = d*(1.0+s*fyz)*bxyz
    yyz = d*(1.0+s*fyz)*byyz
    zyz = d*(1.0+s*fyz)*bzyz

    figte2d = make_subplots(rows=1, cols=3,
                    specs=[[{'type': 'xy'}, {'type': 'xy'}, {'type': 'xy'}]],
                    subplot_titles=['Plot XY plane (z=0)', 'Plot XZ plane (y=0)', 'Plot YZ plane (x=0)'],)

    figte2d.add_trace(go.Scatter(x=xxy, y=yxy, mode='lines',name='(lx,ly)'), 1, 1)
    figte2d.add_trace(go.Scatter(x=bxxy, y=byxy, mode='lines',name='(lox,loy)'), 1, 1)

    figte2d.add_trace(go.Scatter(x=xxz, y=zxz, mode='lines',name='(lx,lz)'), 1, 2)
    figte2d.add_trace(go.Scatter(x=bxxz, y=bzxz, mode='lines',name='(lox,loz)'), 1, 2)

    figte2d.add_trace(go.Scatter(x=yyz, y=zyz, mode='lines',name='(ly,lz)'), 1, 3)
    figte2d.add_trace(go.Scatter(x=byyz, y=bzyz, mode='lines',name='(loy,loz)'), 1, 3)
    figte2d.update_layout(xaxis_title='X-axis')
    figte2d.update_layout(yaxis_title='Y-axis')
    figte2d.update_layout(transition_duration=500)

    figte2d.update_yaxes(automargin=True)
    figte2d.update_xaxes(automargin=True)

    return figte2d



############## Ort-3D

@app.callback(
    Output('ort_3D', 'figure'),
    [Input(component_id='ofieldx', component_property='value'),
     Input(component_id='ofieldy', component_property='value'),
     Input(component_id='ofieldz', component_property='value'),
     Input(component_id='oL01', component_property='value'),
     Input(component_id='oL02', component_property='value'),
     Input(component_id='oL03', component_property='value'),
     Input(component_id='oL1', component_property='value'),
     Input(component_id='oL2', component_property='value'),
     Input(component_id='oL3', component_property='value'),
     Input(component_id='oL4', component_property='value'),
     Input(component_id='oL5', component_property='value'),
     Input(component_id='oL6', component_property='value'),
     Input(component_id='oL7', component_property='value'),
     Input(component_id='oL8', component_property='value'),
     Input(component_id='oL9', component_property='value'),
     Input(component_id='oscale', component_property='value'),
    ]
)


def update_ofig(hx,hy,hz,lmb01,lmb02,lmb03,lmb1,lmb2,lmb3,lmb4,lmb5,lmb6,lmb7,lmb8,lmb9,s):

    d=1.0
    h=np.sqrt(hx*hx+hy*hy+hz*hz)
    ax=hx/h
    ay=hy/h
    az=hz/h
    u, v = np.mgrid[0:np.pi:100j, 0:2*np.pi:100j]
    bx = np.sin(u)*np.cos(v)
    by = np.sin(u)*np.sin(v)
    bz = np.cos(u)
    f = lmb01*bx**2+lmb02*by**2+lmb03*bz**2+lmb1*(ax**2*bx**2-ax*ay*bx*by-ax*az*bx*bz)+lmb2*(ay**2*bx**2-ax*ay*bx*by)+lmb3*(ax**2*by**2-ax*ay*bx*by)+lmb4*(ay**2*by**2-ax*ay*bx*by-ay*az*by*bz)+lmb5*(ax**2*bz**2-ax*az*bx*bz)+lmb6*(ay**2*bz**2-ay*az*by*bz)+lmb7*4*ax*ay*bx*by+lmb8*4*ax*az*bx*bz+lmb9*4*ay*az*by*bz
    x = d*(1.0+s*f)*bx
    y = d*(1.0+s*f)*by
    z = d*(1.0+s*f)*bz
    dl_l=(np.sqrt(x**2+y**2+z**2)-d)/d

    ofig = make_subplots(rows=1, cols=2,
                    specs=[[{'is_3d': True}, {'is_3d': True}]],
                    subplot_titles=['      Magnetic field (H), magnetization (\u03B1) and unit cell lattice vectors (a,b,c)           ', '        Color corresponds to (\u0394l/lo)*scaling_factor along direction \u03B2'],)

    ofig.add_trace(go.Cone(x=[1], y=[1], z=[1], u=[ax], v=[ay], w=[az],name="H",colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(255,0,0)']]),1, 1)
    ofig.add_trace(go.Cone(x=[1], y=[1], z=[0], u=[ax], v=[ay], w=[az],name="\u03B1",colorscale=[[0, 'rgb(255,127,14)'], [1, 'rgb(255,127,16)']]),1, 1)
    ofig.add_trace(go.Cone(x=[0.4], y=[0], z=[0], u=[2], v=[0], w=[0],name="a",colorscale=[[0, 'rgb(255,0,0)'], [1, 'rgb(240,0,0)']]),1, 1)
    ofig.add_trace(go.Cone(x=[0], y=[0.4], z=[0], u=[0], v=[2], w=[0],name="b",colorscale=[[0, 'rgb(0,255,0)'], [1, 'rgb(0,240,0)']]),1, 1)
    ofig.add_trace(go.Cone(x=[0], y=[0], z=[0.4], u=[0], v=[0], w=[2],name="c",colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(0,0,240)']]),1, 1)

    ofig.update_traces(hoverinfo="name", showscale=False)
    ofig.add_trace(go.Surface(x=x, y=y, z=z, surfacecolor=dl_l, name="l"), 1, 2)

    ofig.update_layout(transition_duration=500)

    ofig.update_yaxes(automargin=True)
    ofig.update_xaxes(automargin=True)

    return ofig

######################################## Ort 2D

@app.callback(
    Output('ort_2D', 'figure'),
    [Input(component_id='ofieldx', component_property='value'),
     Input(component_id='ofieldy', component_property='value'),
     Input(component_id='ofieldz', component_property='value'),
     Input(component_id='oL01', component_property='value'),
     Input(component_id='oL02', component_property='value'),
     Input(component_id='oL03', component_property='value'),
     Input(component_id='oL1', component_property='value'),
     Input(component_id='oL2', component_property='value'),
     Input(component_id='oL3', component_property='value'),
     Input(component_id='oL4', component_property='value'),
     Input(component_id='oL5', component_property='value'),
     Input(component_id='oL6', component_property='value'),
     Input(component_id='oL7', component_property='value'),
     Input(component_id='oL8', component_property='value'),
     Input(component_id='oL9', component_property='value'),
     Input(component_id='oscale', component_property='value'),
    ]
)


def update_figo2d(hx,hy,hz,lmb01,lmb02,lmb03,lmb1,lmb2,lmb3,lmb4,lmb5,lmb6,lmb7,lmb8,lmb9,s):

    d=1.0
    h=np.sqrt(hx*hx+hy*hy+hz*hz)
    ax=hx/h
    ay=hy/h
    az=hz/h
    vxy = np.mgrid[0:2*np.pi:200j]
    bxxy = np.cos(vxy)
    byxy = np.sin(vxy)
    bzxy = 0.0
    fxy = lmb01*bxxy**2+lmb02*byxy**2+lmb03*bzxy**2+lmb1*(ax**2*bxxy**2-ax*ay*bxxy*byxy-ax*az*bxxy*bzxy)+lmb2*(ay**2*bxxy**2-ax*ay*bxxy*byxy)+lmb3*(ax**2*byxy**2-ax*ay*bxxy*byxy)+lmb4*(ay**2*byxy**2-ax*ay*bxxy*byxy-ay*az*byxy*bzxy)+lmb5*(ax**2*bzxy**2-ax*az*bxxy*bzxy)+lmb6*(ay**2*bzxy**2-ay*az*byxy*bzxy)+lmb7*4*ax*ay*bxxy*byxy+lmb8*4*ax*az*bxxy*bzxy+lmb9*4*ay*az*byxy*bzxy
    xxy = d*(1.0+s*fxy)*bxxy
    yxy = d*(1.0+s*fxy)*byxy
    zxy = d*(1.0+s*fxy)*bzxy

    vxz = np.mgrid[0:2*np.pi:200j]
    bxxz = np.sin(vxz)
    byxz = 0.0
    bzxz = np.cos(vxz)
    fxz = lmb01*bxxz**2+lmb02*byxz**2+lmb03*bzxz**2+lmb1*(ax**2*bxxz**2-ax*ay*bxxz*byxz-ax*az*bxxz*bzxz)+lmb2*(ay**2*bxxz**2-ax*ay*bxxz*byxz)+lmb3*(ax**2*byxz**2-ax*ay*bxxz*byxz)+lmb4*(ay**2*byxz**2-ax*ay*bxxz*byxz-ay*az*byxz*bzxz)+lmb5*(ax**2*bzxz**2-ax*az*bxxz*bzxz)+lmb6*(ay**2*bzxz**2-ay*az*byxz*bzxz)+lmb7*4*ax*ay*bxxz*byxz+lmb8*4*ax*az*bxxz*bzxz+lmb9*4*ay*az*byxz*bzxz
    xxz = d*(1.0+s*fxz)*bxxz
    yxz = d*(1.0+s*fxz)*byxz
    zxz = d*(1.0+s*fxz)*bzxz

    vyz = np.mgrid[0:2*np.pi:200j]
    bxyz = 0.0
    byyz = np.sin(vyz)
    bzyz = np.cos(vyz)
    fyz = lmb01*bxyz**2+lmb02*byyz**2+lmb03*bzyz**2+lmb1*(ax**2*bxyz**2-ax*ay*bxyz*byyz-ax*az*bxyz*bzyz)+lmb2*(ay**2*bxyz**2-ax*ay*bxyz*byyz)+lmb3*(ax**2*byyz**2-ax*ay*bxyz*byyz)+lmb4*(ay**2*byyz**2-ax*ay*bxyz*byyz-ay*az*byyz*bzyz)+lmb5*(ax**2*bzyz**2-ax*az*bxyz*bzyz)+lmb6*(ay**2*bzyz**2-ay*az*byyz*bzyz)+lmb7*4*ax*ay*bxyz*byyz+lmb8*4*ax*az*bxyz*bzyz+lmb9*4*ay*az*byyz*bzyz
    xyz = d*(1.0+s*fyz)*bxyz
    yyz = d*(1.0+s*fyz)*byyz
    zyz = d*(1.0+s*fyz)*bzyz

    figo2d = make_subplots(rows=1, cols=3,
                    specs=[[{'type': 'xy'}, {'type': 'xy'}, {'type': 'xy'}]],
                    subplot_titles=['Plot XY plane (z=0)', 'Plot XZ plane (y=0)', 'Plot YZ plane (x=0)'],)

    figo2d.add_trace(go.Scatter(x=xxy, y=yxy, mode='lines',name='(lx,ly)'), 1, 1)
    figo2d.add_trace(go.Scatter(x=bxxy, y=byxy, mode='lines',name='(lox,loy)'), 1, 1)

    figo2d.add_trace(go.Scatter(x=xxz, y=zxz, mode='lines',name='(lx,lz)'), 1, 2)
    figo2d.add_trace(go.Scatter(x=bxxz, y=bzxz, mode='lines',name='(lox,loz)'), 1, 2)

    figo2d.add_trace(go.Scatter(x=yyz, y=zyz, mode='lines',name='(ly,lz)'), 1, 3)
    figo2d.add_trace(go.Scatter(x=byyz, y=bzyz, mode='lines',name='(loy,loz)'), 1, 3)
    figo2d.update_layout(xaxis_title='X-axis')
    figo2d.update_layout(yaxis_title='Y-axis')
    figo2d.update_layout(transition_duration=500)

    figo2d.update_yaxes(automargin=True)
    figo2d.update_xaxes(automargin=True)

    return figo2d


############## rod-3D

@app.callback(
    Output('rod', 'figure'),
    [Input(component_id='RR', component_property='value'),
     Input(component_id='LL', component_property='value'),
     Input(component_id='hlong', component_property='value'),
     Input(component_id='lmbs', component_property='value'),
     Input(component_id='imax', component_property='value'),
    ]
)


def update_rod(rad,length,hlong,llmbs,iimax):

    ss = 4.0*np.pi*rad**2.0
    phi = (3.0*llmbs*length*iimax)/(ss*hlong*2.0)
 
    s = np.linspace(0, 2 * np.pi, 100)
    t = np.linspace(0, 1, 100)
    sGrid, tGrid = np.meshgrid(s, t)
    x = rad * np.cos(sGrid) 
    y = rad * np.sin(sGrid)  
    z = length*tGrid                   
    
    ht=iimax/(2.0*np.pi*np.sqrt((4.0*rad)**2.0+(0.5*length)**2.0))
    
    hln=(0.15*length*hlong)/(np.sqrt(hlong**2+ht**2))
    htn=(0.15*length*ht)/(np.sqrt(hlong**2+ht**2))
    
    #hln=hlong
    #htn=ht
    
    

    figrod = make_subplots(rows=1, cols=1,
                    specs=[[{'is_3d': True}]],
                    subplot_titles=['Electric current I=Imax. Color of the rod corresponds to the twisted angle ϕ(z) in radian.'],)

    figrod.add_trace(go.Cone(x=[4*rad], y=[0], z=[0.5*length], u=[0], v=[htn], w=[hln],name="H∥+H⟂",sizemode="absolute",showscale=False,colorscale=[[0, 'rgb(0,255,0)'], [1, 'rgb(255,0,0)']]),1, 1)
    figrod.add_trace(go.Cone(x=[-4*rad], y=[0], z=[0.5*length], u=[0], v=[-htn], w=[hln],name="H∥+H⟂",sizemode="absolute",showscale=False,colorscale=[[0, 'rgb(0,255,0)'], [1, 'rgb(255,0,0)']]),1, 1)
    figrod.add_trace(go.Cone(x=[0], y=[4*rad], z=[0.5*length], u=[-htn], v=[0], w=[hln],name="H∥+H⟂",sizemode="absolute",showscale=False,colorscale=[[0, 'rgb(0,255,0)'], [1, 'rgb(255,0,0)']]),1, 1)
    figrod.add_trace(go.Cone(x=[0], y=[-4*rad], z=[0.5*length], u=[htn], v=[0], w=[hln],name="H∥+H⟂",sizemode="absolute",showscale=False,colorscale=[[0, 'rgb(0,255,0)'], [1, 'rgb(255,0,0)']]),1, 1)
    figrod.update_traces(hoverinfo="name")
    figrod.add_trace(go.Surface(x=x, y=y, z=z, surfacecolor=phi*(z/length), name="magnetic rod"), 1, 1)

    #figrod.update_traces(hoverinfo="name")
    figrod.update_layout(transition_duration=500)

    figrod.update_yaxes(automargin=True)
    figrod.update_xaxes(automargin=True)

    return figrod


############################################### twisted angle


@app.callback(
    Output('twist_angle', 'figure'),
    [Input(component_id='RR', component_property='value'),
     Input(component_id='LL', component_property='value'),
     Input(component_id='hlong', component_property='value'),
     Input(component_id='lmbs', component_property='value'),
     Input(component_id='imin', component_property='value'),
     Input(component_id='imax', component_property='value'),
    ]
)


def update_twist(rad,leng,hl,llmbs,iimin,iimax):

    u = np.mgrid[0:1:100j]   
    i = iimin+(iimax-iimin)*u
    s = 4.0*np.pi*rad**2.0
    
    phi = (3.0*llmbs*leng*i)/(s*hl*2.0)

    tfig = make_subplots(rows=1, cols=1,
                    specs=[[{'type': 'xy'}]],
                    subplot_titles=['Twisted angle at z=L vs electric current'],)

    tfig.add_trace(go.Scatter(x=i, y=phi, mode='lines', name='(I,ϕ)'), 1, 1)


    tfig.update_layout(xaxis_title='Electric current I (A)')
    tfig.update_layout(yaxis_title='Twisted angle ϕ(z=L) (rad)')
    tfig.update_layout(transition_duration=500)

    tfig.update_yaxes(automargin=True)
    tfig.update_xaxes(automargin=True)

    return tfig


############################################### torque


@app.callback(
    Output('torque', 'figure'),
    [Input(component_id='RR', component_property='value'),
     Input(component_id='LL', component_property='value'),
     Input(component_id='hlong', component_property='value'),
     Input(component_id='lmbs', component_property='value'),
     Input(component_id='imin', component_property='value'),
     Input(component_id='imax', component_property='value'),
     Input(component_id='y', component_property='value'),
     Input(component_id='sigma', component_property='value'),
    ]
)


def update_torque(rad,leng,hl,llmbs,iimin,iimax,yy,sig):

    u = np.mgrid[0:1:100j]   
    i = iimin+(iimax-iimin)*u
    s = 4.0*np.pi*rad**2.0
    
    phi = (3.0*llmbs*leng*i)/(s*hl*2.0)
    
    mu = (yy*10.0**9)/(2.0*(1.0+sig))
    
    tor = (mu*np.pi*rad**4.0*phi)/(2.0*leng)


    tofig = make_subplots(rows=1, cols=1,
                    specs=[[{'type': 'xy'}]],
                    subplot_titles=['Helical field-induced torque vs electric current'],)

    tofig.add_trace(go.Scatter(x=i, y=tor, mode='lines', name='(I,τ)'), 1, 1)


    tofig.update_layout(xaxis_title='Electric current I (A)')
    tofig.update_layout(yaxis_title='Helical field-induced torque τ (N*m)')
    tofig.update_layout(transition_duration=500)

    tofig.update_yaxes(automargin=True)
    tofig.update_xaxes(automargin=True)

    return tofig


if __name__ == '__main__':
        app.run_server(debug=True)
