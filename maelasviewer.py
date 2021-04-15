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
    html.H3("Source files"),
    dcc.Markdown('''[https://github.com/pnieves2019/MAELASviewer](https://github.com/pnieves2019/MAELASviewer)'''),
    html.Hr(),


])


################## Select system to show

@app.callback(
    dash.dependencies.Output('effect_output', 'children'),
    [dash.dependencies.Input('effect', 'value')])

def update_output(effect):

    if effect == 'Joule':

        return html.Div(children=[

            html.H3("The Joule effect"),
            html.H4("Mapping the Joule effect to a unit sphere"),
            dcc.Markdown('''This interactive tool maps the Joule effect to a unit sphere. You can visualize the relative length change \u0394l/lo=(l-lo)/lo of the material along an arbitrary measuring direction **β** as a function of the external magnetic field **H**, magnetocrystalline anisotropy (K) and magnetostrictive coefficients (λ). The quantity lo corresponds to the length of the magnetic material in the demagnetized state along the measuring direction **β**=(sinθ*cosφ, sinθ*sinφ, cosθ), where θ and φ are the polar and azimuthal angles, respectively. In the simulation, the length lo of the material in any measuring direction **β** is represented by a sphere with radius 1 (lo=1), so that it does not depend on the geometry of the material. In the magnetized state, the unit sphere is distorted, and the distance between a point on the surface and the origin (0, 0, 0) describes the simulated length (lsim) along the measuring direction **β**. The magnetization direction is given by **α**, and it is assumed that the magnetization is saturated (M=Ms) along the effective field **Heff**=**H**+**Ha**, where **Ha** is the magnetocrystalline anisotropy field.'''),
            html.Div([html.Img(src=app.get_asset_url('diagram_online.png'))]),
            dcc.Markdown('''The real length along the measuring direction **β** of a material with any shape in the magnetized state (lexp) can be easily obtained with this interactive applet by setting the magnetic field, the corresponding values of the magnetocrystalline anisotropy constants and magnetostrictive coefficients, and using the following equation: '''),
            html.Div([html.Img(src=app.get_asset_url('convert.png'))]),
            dcc.Markdown('''where lo,exp is the real length along the measuring direction **β** of a material in the demagnetized state, lsim is the simulated length value printed by the interactive app at the surface of the distorted unit sphere along direction **β**, and s is the scale factor set by the user to facilitate the visualization. When the user click on the surface of the distorted sphere, it prints the components of vector **l**sim, so that it gives both the length lsim=|**l**sim|, and direction **β** since it is aligned to **l**sim (**β**‖**l**sim) and |**β**|=1. '''),
            html.Div([html.Img(src=app.get_asset_url('convert_fig.png'))]),
            dcc.Markdown(''' It is important to note that this method gives only the length change in the measuring direction **β** but not the overall real shape deformation of a material induced by the Joule effect, not even for a sphere. The real shape deformation of a sphere due to Joule effect corresponds to an ellipsoid. This shape is different to the shape of the distorted sphere generated by mapping the Joule effect to a unit sphere that is implemented in MAELASviewer. The reason for that is the lack of transverse component of the displacement vector in this procedure (**u**⊥=0, **u**‖=**u**).'''),
            html.Div([html.Img(src=app.get_asset_url('comparison_sphere.png'))]),


            html.H4("Example"),


            dcc.Markdown(''' Let's consdier that we have a single crystal FCC Ni with lattice vectors are **a**=(a0,0,0), **b**=(0,a0,0) and **c**=(0,0,a0), and its length along the measuring direction **β**=(1,0,0) in the demagnetized state is 1µm (lo,exp=1µm) at room temperature (T=300K). The magnetization is saturated along an applied magnetic field of 2 Tesla in the direction (0,0,1). How long is the length of the magnetized material in the measuring direction **β**=(1,0,0)?  '''),
            dcc.Markdown('''**Solution** '''),
            dcc.Markdown(''' Firstly, we select the crystal system `Single crystal: Cubic I (space group numbers  207-230)` since FCC Ni is a cubic crystal with space group 225. Next, we set the components of the external magnetic field μ0**H**=(0,0,2)T, and the corresponding values of the magnetostrictive coefficients, magnetocrystalline anisotropy constants and saturation mangetization at room temperature of FCC Ni (λ001=-0.000046, λ111=-0.000024, K1=-0.005MJ/m^3, K2=-0.002MJ/m^3, μ0Ms=0.61T). In this example, we don't include the volume magnetostriction, so we set λα=0. In order to facilitate the visualization of magnetostriction we set the scale factor s=10000. Now, either in the 3D surface plot or in the 2D cross section plot in the plane XY (z=0) we can click on the measuring direction **β**=(1,0,0)  to read the value of the simulated length along **β**. Doing so, we get lsim=1.23. Finally, inserting lsim=1.23, s=10000 and lo,exp=1µm into Eq.(1), we obtain that the final length along **β**=(1,0,0) is equal to lexp=1.000023µm. In this case, we see that the magnetic field induced a very small length change in the scale of interatomic distances. The fractional change length is \u0394l/lo=0.000023 which is also provided by the interactive figures. In the 3D figure of the magnetocrystalline anisotropy energy  we can see the magnetic easy direction (111) and hard direction (100). '''),

            html.Hr(),
            html.H3("Available systems"),
            html.H6("1. Single crystal: Cubic I (space group numbers  207-230)"),
            html.H6("2. Polycrystal: Cubic I (space group numbers  207-230)"),
            html.H6("3. Single crystal: Hexagonal I (space group numbers  177-194) - Clark definition"),
            html.H6("4. Single crystal: Hexagonal I (space group numbers  177-194) - Mason definition"),
            html.H6("5. Single crystal: Trigonal I (space group numbers  149-167)"),
            html.H6("6. Single crystal: Tetragonal I (space group numbers  89-142)"),
            html.H6("7. Single crystal: Orthorhombic (space group numbers  16-74)"),
            html.Hr(),
            html.H3('Select a crystal system:'),
            dcc.Dropdown(
                id='joule_system',
                options=[
                    {'label': 'Single crystal: Cubic I (space group numbers  207-230)', 'value': 'cub'},
                    {'label': 'Polycrystal: Cubic I (space group numbers  207-230)', 'value': 'cubpol'},
                    {'label': 'Single crystal: Hexagonal I (space group numbers  177-194) - Clark definition', 'value': 'hex'},
                    {'label': 'Single crystal: Hexagonal I (space group numbers  177-194) - Mason definition', 'value': 'hexm'},
                    {'label': 'Single crystal: Trigonal I (space group numbers  149-167)', 'value': 'trig'},
                    {'label': 'Single crystal: Tetragonal I (space group numbers  89-142)', 'value': 'tet'},
                    {'label': 'Single crystal: Orthorhombic (space group numbers  16-74)', 'value': 'ort'},
                ],
                placeholder="Select a crystal system"
            ),
            html.Div(id='joule_system_output')

        ])

    elif effect == 'Wied':

        return html.Div(children=[


            html.H3("The Wiedemann effect"),
            html.H4("Theory"),


            html.H6("This interactive applet simulates the Wiedemann effect for a magnetostrictive isotropic cylindrical rod (the twist of a rod induced by helical magnetic field). The helical magnetic field that twists the magnetic rod is achieved by combining a magnetic field along the rod height axis (longitudinal field H‖), and an electric current (I) through the rod which induces a circular magnetic field due to Ampère’s circuital law (perpendicular field H⊥)."),
            html.Div([html.Img(src=app.get_asset_url('torsion_rod.png'))]),
            html.H6("For an isotropic magnetic cylindrical rod aligned to the z-axis with height L and radius R, the twisted angle ϕ induced by a helical magnetic field is given by"),
            html.Div([html.Img(src=app.get_asset_url('eq_twist.png'))]),
            html.H6("where \u03BBs is the isotropic magnetostrictive coefficient and A = 4πR² is the area of the cross section of the rod. The helical field-induced torque can be calculated as"),
            html.Div([html.Img(src=app.get_asset_url('eq_torque.png'))]),
            html.H6("where Y is the Young’s modulus and σ is the Poisson’s ratio. In the 3D visualization of the rod, the perpendicular and longitudinal magnetic fields in the exterior of the rod at z = L/2 are also plotted. The perpendicular field is calculated applying the Biot-Savart law for a finite wire. Note that the magnetic field generated by the magnetization of the rod is not plotted here."),

            html.H4("Parameters of the simulation"),
            html.H6("(Press Enter after changing any input to update the figures)"),
            html.Hr(),
            dcc.Markdown(''' **Geometry of the magnetic rod:**'''),
            html.Div(['L(m) = ',
              dcc.Input(id='LL', value=0.0005, type='number', debounce=True, step=0.000000001)]),
            html.Div(['R(m) = ',
              dcc.Input(id='RR', value=0.00001, type='number', debounce=True, step=0.000000001)]),
            dcc.Markdown(''' **Longitudinal external magnetic field:**'''),
            html.Div(['H',html.Sub('||'),"(A/m) = ",
              dcc.Input(id='hlong', value=0.000001, type='number', debounce=True, step=0.000000001)]),
            dcc.Markdown(''' **Isotropic magnetostrictive coefficient:**'''),
            html.Div(['\u03BB',html.Sub('S')," = ",
              dcc.Input(id='lmbs', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            dcc.Markdown(''' **Range of applied electric current:**'''),
            html.Div(['I',html.Sub('min'),"(A) = ",
              dcc.Input(id='imin', value=-0.00000001, type='number', debounce=True, step=0.000000001)]),
            html.Div(['I',html.Sub('max'),"(A) = ",
              dcc.Input(id='imax', value=0.00000001, type='number', debounce=True, step=0.000000001)]),
            dcc.Markdown(''' **Isotropic elastic properties:**'''),
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



        ])
    else:

        return 'No magnetostriction effect selected'


@app.callback(
    dash.dependencies.Output('joule_system_output', 'children'),
    [dash.dependencies.Input('joule_system', 'value')])

def update_output(system):

    if system == 'cub':

        return html.Div(children=[




            html.H3("Single crystal: Cubic I (space group numbers  207-230)"),
            html.H4("Theory"),
            html.H6("The relative length change for cubic (I) systems is given by:"),
            html.Div([html.Img(src=app.get_asset_url('eq_cub.png'))]),
            html.H6("where αi and βi (i=x,y,z) are the direction of magnetization and the measured length direction, respectively."),
            html.H6("Magentostriction is a small effect that is hard to visualize. To facilitate its visualization in the simulation, we multiply the right hand side of this equation by a scale factor parameter (s) that can be modified by the user"),
            html.Div([html.Img(src=app.get_asset_url('eq_cub_s1.png'))]),
            html.H6("Note that this scale preserve the ratio between the magnetostrictive coefficients. Obviously, the case with s=1 corresponds to the real situation. The printed length lsim given at the surface of the distorted sphere is calculated as"),
            html.Div([html.Img(src=app.get_asset_url('eq_cub_s2.png'))]),
            html.H6("where we took into account that lo=1. Similar procedure is applied to the other supported crystal systems."),

            html.H6("The magnetocrystalline anisotropy energy for cubic systems is"),
            html.Div([html.Img(src=app.get_asset_url('mae_cub.png'))]),
            html.H6("where K0, K1 and K2 are the magnetocrystalline anisotropy constants. K0 can be used to shift the minimum energy reference in the 3D visualization of the magnetocrystalline anisotropy energy. The magnetocrystalline anisotropy field Ha is"),
            html.Div([html.Img(src=app.get_asset_url('mae_field.png'))]),
            html.H6("where μ0 is the vacuum permeability and Ms is the saturation magnetization. In the simulation, it is assumed that the magnetization is saturated along the effective field Heff"),
            html.Div([html.Img(src=app.get_asset_url('eff_field.png'))]),
            html.H6("where H is the external magnetic field. The direction of the equilibrium magnetization α is calculated by the Landau-Lifshitz-Gilbert (LLG) equation"),
            html.Div([html.Img(src=app.get_asset_url('llg.png'))]),
            html.H6("where γ is the electron gyromagnetic ratio, and η is the damping parameter. The torque |α x μ0Heff| is used as the criterion for the numerical convergence of the equilibrium state, it is recommended to use a tolerance for the torque lower than 0.0001 Tesla. Default values of the damping parameter, time step (dt) and total number of iteration steps usually work fine, but the user can change them in case the tolerance for the torque is not achieved."),
            html.Div([html.Img(src=app.get_asset_url('heffc.png'))]),


            html.H4("Parameters of the simulation"),
            html.H6("(Press Enter after changing any input to update the figures)"),
            html.Hr(),
            dcc.Markdown(''' **External magnetic field:**'''),
            html.Div(['μ0H',html.Sub('x'),'(Tesla)'" = ",
              dcc.Input(id='fieldx', value=5.000, type='number', debounce=True, step=0.0000001)]),
            html.Div(['μ0H',html.Sub('y'),'(Tesla)'" = ",
              dcc.Input(id='fieldy', value=0.500, type='number', debounce=True, step=0.0000001)]),
            html.Div(['μ0H',html.Sub('z'),'(Tesla)'" = ",
              dcc.Input(id='fieldz', value=0.005, type='number', debounce=True, step=0.0000001)]),
            dcc.Markdown(''' **Magnetostrictive coefficients:**'''),
            html.Div(['\u03BB',html.Sup("\u03B1")," = ",
              dcc.Input(id='L0', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sub("001")," = ",
              dcc.Input(id='L1', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sub("111")," = ",
              dcc.Input(id='L2', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(["Scale factor = ",
              dcc.Input(id='scale', value=1.0, type='number', debounce=True)]),
            dcc.Markdown(''' **Magnetocrystalline anisotropy constants:**'''),
            html.Div(['K',html.Sub("0"), "(MJ/m",html.Sup("3"),")"," = ",
              dcc.Input(id='k0c', value=0.001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['K',html.Sub("1"), "(MJ/m",html.Sup("3"),")"," = ",
              dcc.Input(id='k1c', value=0.001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['K',html.Sub("2"),"(MJ/m",html.Sup("3"),")"," = ",
              dcc.Input(id='k2c', value=0.001, type='number', debounce=True, step=0.0000001)]),
            dcc.Markdown(''' **Saturation magnetization:**'''),
            html.Div(['μ0Ms (Tesla)'," = ",
              dcc.Input(id='ms', value=0.01, type='number', debounce=True, step=0.0001)]),
            dcc.Markdown(''' **Landau-Lifshitz-Gilbert solver to find the direction of the equilibrium magnetization (α):**'''),
            html.Div(['Damping LLG '," = ",
              dcc.Input(id='alpha0', value=0.95, type='number', debounce=True, step=0.00001)]),
            html.Div(['Torque tolerance |α x μ0Heff| (Tesla)'," = ",
              dcc.Input(id='tol', value=0.00001, type='number', debounce=True, step=0.000000001)]),
            html.Div(['Time step LLG (ps)'," = ",
              dcc.Input(id='dt', value=0.01, type='number', debounce=True, step=0.00001)]),
            html.Div(['Total number of iteration steps'," = ",
              dcc.Input(id='nt', value=500000, type='number', debounce=True, step=1)]),


            html.Div(id='my-output'),
            html.Hr(),
            html.H4("Simulation"),

            dcc.Markdown(''' **LLG solver output:**'''),
            html.Table([
                html.Tr([html.Td(['Parameters']), html.Td(['Calculated values']),html.Td(['Parameters']), html.Td(['Calculated values']),html.Td(['Parameters']), html.Td(['Calculated values'])]),
                html.Tr([html.Td(['α', html.Sub('x')]), html.Td(id='meqxc'),html.Td(['μ0H', html.Sub('eff,x'),'(T)']), html.Td(id='hefxc'),html.Td(['μ0H', html.Sub('a,x'),'(T)']), html.Td(id='haxc')]),
                html.Tr([html.Td(['α', html.Sub('y')]), html.Td(id='meqyc'),html.Td(['μ0H', html.Sub('eff,y'),'(T)']), html.Td(id='hefyc'),html.Td(['μ0H', html.Sub('a,y'),'(T)']), html.Td(id='hayc')]),
                html.Tr([html.Td(['α', html.Sub('z')]), html.Td(id='meqzc'),html.Td(['μ0H', html.Sub('eff,z'),'(T)']), html.Td(id='hefzc'),html.Td(['μ0H', html.Sub('a,z'),'(T)']), html.Td(id='hazc')]),
                html.Tr([html.Td(['Torque |α x μ0Heff| (T)']), html.Td(id='torquec')]),
                html.Tr([html.Td(['Status']), html.Td(id='statc')]),

            ]),
            html.Hr(),
            dcc.Markdown(''' **Generated figures to visualize magnetostriction:** '''),
            html.H6("The distance between a point on the surface and the origin (0,0,0) describes the simulated length lsim in the direction =(sinθ*cosφ, sinθ*sinφ, cosθ), where θ and φ are the polar and azimuthal angles, respectively. The color of the surface corresponds to the relative length change multiplied by the scale factor in the measuring length direction \u03B2."),

            dcc.Graph(id='cub_3D'),
            dcc.Graph(id='cub_2D'),
            html.Hr(),
            dcc.Markdown(''' **Generated figure to visualize the magnetocrystalline anisotropy energy:**'''),
            html.H6("In this interactive plot you can visualize the magnetocrystalline anisotropy energy (Ea) as a function of the magnetization direction α=(sinθ*cosφ, sinθ*sinφ, cosθ)."),

            dcc.Graph(id='mae_cub_3D'),
            html.Hr(),
            html.H4("Magnetic properties of some cubic crystals"),

            html.Table([
                html.Tr([html.Td('Material'), html.Td('Space group'),html.Td('Temperature (K)'),html.Td(html.Div(['    \u03BB',html.Sup("\u03B1"),'    '])), html.Td(html.Div(['\u03BB',html.Sub("001")])),html.Td(html.Div(['\u03BB',html.Sub("111")])),html.Td(html.Div(['K',html.Sub("1"),'(MJ/m',html.Sub("3"),')'])),html.Td(html.Div(['K',html.Sub("2"),'(MJ/m',html.Sub("3"),')'])),html.Td(html.Div(['μ0Ms(Tesla)']))]),
                html.Tr([html.Td('BCC Fe'), html.Td('229'),html.Td('0'),html.Td('    -   '),html.Td('0.000026'), html.Td('-0.00003'),html.Td('0.052'),html.Td('-0.018'),html.Td('2.19')]),
                html.Tr([html.Td('BCC Fe'), html.Td('229'),html.Td('300'),html.Td('   -   '),html.Td('0.000021'), html.Td('-0.000021'),html.Td('0.048'),html.Td('-0.010'),html.Td('2.14')]),
                html.Tr([html.Td('FCC Ni'), html.Td('225'),html.Td('0'),html.Td('   -   '),html.Td('-0.00006'), html.Td('-0.000035'),html.Td('-0.120'),html.Td('0.030'),html.Td('0.64')]),
                html.Tr([html.Td('FCC Ni'), html.Td('225'),html.Td('300'),html.Td('   -   '),html.Td('-0.000046'), html.Td('-0.000024'),html.Td('-0.005'),html.Td('-0.002'),html.Td('0.61')]),
                html.Tr([html.Td(html.Div(['Laves phase SmFe',html.Sub("2")])), html.Td('227'),html.Td('0'),html.Td('   -   '),html.Td('0.00003'), html.Td('-0.0041'),html.Td('   -   '),html.Td('   -   '),html.Td('   -   ')]),
                html.Tr([html.Td(html.Div(['Laves phase DyFe',html.Sub("2")])), html.Td('227'),html.Td('0'),html.Td('   -   '),html.Td('-0.00007'), html.Td('0.003'),html.Td('   -   '),html.Td('   -   '),html.Td('   -   ')]),
                html.Tr([html.Td(html.Div(['Laves phase TbCo',html.Sub("2")])), html.Td('227'),html.Td('0'),html.Td('   -   '),html.Td('-0.0012'), html.Td('0.0045'),html.Td('   -   '),html.Td('   -   '),html.Td('   -   ')]),
                html.Tr([html.Td(html.Div(['Laves phase ErCo',html.Sub("2")])), html.Td('227'),html.Td('0'),html.Td('   -   '),html.Td('-0.001'), html.Td('-0.0025'),html.Td('   -   '),html.Td('   -   '),html.Td('   -   ')]),
            ])
        ])


    elif system == 'cubpol':

        return html.Div(children=[

            html.H3("Polycrystal: Cubic I (space group numbers  207-230)"),
            html.H4("Theory"),
            html.H6("The theory of magnetostriction for polycrystalline materials is more complex. A widely used approximation is to assume that the stress distribution is uniform through the material. In this case the relative change in length may be put into the form:"),
            html.Div([html.Img(src=app.get_asset_url('eq_cub_poly.png'))]),
            html.H6("where"),
            html.Div([html.Img(src=app.get_asset_url('eq_cub_poly_lmb_s.png'))]),
            html.H6("The magnetization is considered to be saturated in the direction of the effective field Heff"),
            html.Div([html.Img(src=app.get_asset_url('eff_field.png'))]),
            html.H6("The magentocrystalline anisotropy field for a polycrystalline material is more complicated than for single crystals. Hence, in the simulation the user can set directly the final direction of equilibrium magnetization via the effective field Heff."),
            html.Div([html.Img(src=app.get_asset_url('poly.png'))]),

            html.H4("Parameters of the simulation"),
            html.H6("(Press Enter after changing any input to update the figures)"),
            html.Hr(),
            html.H6("Effective magnetic field:"),
            html.Div(['μ0H',html.Sub('eff,x'),"(T) = ",
              dcc.Input(id='pfieldx', value=1.0, type='number', debounce=True, step=0.1)]),
            html.Div(['μ0H',html.Sub('eff,y'),"(T) = ",
              dcc.Input(id='pfieldy', value=0.0, type='number', debounce=True, step=0.1)]),
            html.Div(['μ0H',html.Sub('eff,z'),"(T) = ",
              dcc.Input(id='pfieldz', value=0.0, type='number', debounce=True, step=0.1)]),
            html.H6("Magnetostrictive coefficients:"),
            html.Div(['\u03BB',html.Sub("S")," = ",
              dcc.Input(id='LS', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(["Scale factor = ",
              dcc.Input(id='pscale', value=1.0, type='number', debounce=True)]),

            html.Div(id='my-output-p'),
            html.Hr(),
            html.H4("Simulation"),
            dcc.Markdown(''' **Generated figures to visualize magnetostriction:** '''),
            html.H6("The distance between a point on the surface and the origin (0,0,0) describes the simulated length lsim in the direction \u03B2=(sinθ*cosφ, sinθ*sinφ, cosθ), where θ and φ are the polar and azimuthal angles, respectively. The color of the surface corresponds to the relative length change multiplied by the scale factor in the measuring length direction \u03B2."),


            dcc.Graph(id='graph-poly'),
            dcc.Graph(id='graph-poly2D'),
            html.H4("Magnetostrictive coefficients for some polycrystal materials (cubic I)"),
            html.Table([
                html.Tr([html.Td('Material'), html.Td('Space group'),html.Td('Temperature (K)'),html.Td(html.Div(['\u03BB',html.Sub("S")]))]),
                html.Tr([html.Td('BCC Fe'), html.Td('229'),html.Td('300'),html.Td('-0.000007')]),
                html.Tr([html.Td('FCC Ni'), html.Td('225'),html.Td('300'),html.Td('-0.000034')]),
                html.Tr([html.Td(html.Div(['Laves phase TbFe',html.Sub("2")])), html.Td('227'),html.Td('300'),html.Td('0.001753')]),
            ])
        ])

######## Hex I (Clark)

    elif system == 'hex':

        return html.Div(children=[

            html.H3(" Single crystal: Hexagonal I (space group numbers  177-194) - Clark definition"),
            html.H4("Theory"),
            html.H6("In Clark definition, the relative length change for hexagonal (I) systems is given by:"),
            html.Div([html.Img(src=app.get_asset_url('eq_hex.png'))]),
            html.H6("where αi and βi (i=x,y,z) are the direction of magnetization and the measuring length direction, respectively. Magentostriction is a small effect that is hard to visualize. To facilitate its visualization in the simulation, we multiply the right hand side of this equation by a scale factor parameter which can be modified by the user."),

            html.H6("The magnetocrystalline anisotropy energy for hexagonal systems is"),
            html.Div([html.Img(src=app.get_asset_url('mae_hex.png'))]),
            html.H6("where K0, K1 and K2 are the magnetocrystalline anisotropy constants. K0 can be used to shift the minimum energy reference in the 3D visualization of the magnetocrystalline anisotropy energy. The magnetocrystalline anisotropy field Ha is"),
            html.Div([html.Img(src=app.get_asset_url('mae_field.png'))]),
            html.H6("where μ0 is the vacuum permeability and Ms is the saturation magnetization. In the simulation, it is assumed that the magnetization is saturated along the effective field Heff"),
            html.Div([html.Img(src=app.get_asset_url('eff_field.png'))]),
            html.H6("where H is the external magnetic field. The direction of the equilibrium magnetization α is calculated by the Landau-Lifshitz-Gilbert equation"),
            html.Div([html.Img(src=app.get_asset_url('llg.png'))]),
            html.H6("where γ is the electron gyromagnetic ratio, and η is the damping parameter. The torque |α x μ0Heff| is used as the criterion for the numerical convergence of the equilibrium state, it is recommended to use a tolerance for the torque lower than 0.0001 Tesla. Default values of the damping parameter, time step (dt) and total number of iteration steps usually work fine, but the user can change them in case the tolerance for the torque is not achieved."),
            html.Div([html.Img(src=app.get_asset_url('heffh.png'))]),


            html.H4("Parameters of the simulation"),
            html.H6("(Press Enter after changing any input to update the figures)"),
            html.Hr(),
            dcc.Markdown(''' **External magnetic field:**'''),
            html.Div(['μ0H',html.Sub('x'),'(Tesla)'" = ",
              dcc.Input(id='hfieldx', value=5.000, type='number', debounce=True, step=0.0000001)]),
            html.Div(['μ0H',html.Sub('y'),'(Tesla)'" = ",
              dcc.Input(id='hfieldy', value=0.500, type='number', debounce=True, step=0.0000001)]),
            html.Div(['μ0H',html.Sub('z'),'(Tesla)'" = ",
              dcc.Input(id='hfieldz', value=0.005, type='number', debounce=True, step=0.0000001)]),
            dcc.Markdown(''' **Magnetostrictive coefficients:**'''),
            html.Div(['\u03BB',html.Sup("\u03B11,0")," = ",
              dcc.Input(id='hL01', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sup("\u03B12,0")," = ",
              dcc.Input(id='hL02', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sup("\u03B11,2")," = ",
              dcc.Input(id='hLa1', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sup("\u03B12,2")," = ",
              dcc.Input(id='hLa2', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sup("\u03B3,2")," = ",
              dcc.Input(id='hLg', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sup("\u03B5,2")," = ",
              dcc.Input(id='hLe', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(["Scale factor = ",
              dcc.Input(id='hscale', value=1.0, type='number', debounce=True)]),
            dcc.Markdown(''' **Magnetocrystalline anisotropy constants:**'''),
            html.Div(['K',html.Sub("0"), "(MJ/m",html.Sup("3"),")"," = ",
              dcc.Input(id='k0h', value=0.001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['K',html.Sub("1"), "(MJ/m",html.Sup("3"),")"," = ",
              dcc.Input(id='k1h', value=0.001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['K',html.Sub("2"),"(MJ/m",html.Sup("3"),")"," = ",
              dcc.Input(id='k2h', value=0.001, type='number', debounce=True, step=0.0000001)]),
            dcc.Markdown(''' **Saturation magnetization:**'''),
            html.Div(['μ0Ms (Tesla)'," = ",
              dcc.Input(id='msh', value=0.01, type='number', debounce=True, step=0.0001)]),
            dcc.Markdown(''' **Landau-Lifshitz-Gilbert solver to find the direction of the equilibrium magnetization (α):**'''),
            html.Div(['Damping LLG '," = ",
              dcc.Input(id='alpha0h', value=0.95, type='number', debounce=True, step=0.00001)]),
            html.Div(['Torque tolerance |α x μ0Heff| (Tesla)'," = ",
              dcc.Input(id='tolh', value=0.00001, type='number', debounce=True, step=0.000000001)]),
            html.Div(['Time step LLG (ps)'," = ",
              dcc.Input(id='dth', value=0.01, type='number', debounce=True, step=0.00001)]),
            html.Div(['Total number of iteration steps'," = ",
              dcc.Input(id='nth', value=500000, type='number', debounce=True, step=1)]),

            html.Div(id='my-output-h'),
            html.Hr(),

            html.H4("Simulation"),

            dcc.Markdown(''' **LLG solver output:**'''),
            html.Table([
                html.Tr([html.Td(['Parameters']), html.Td(['Calculated values']),html.Td(['Parameters']), html.Td(['Calculated values']),html.Td(['Parameters']), html.Td(['Calculated values'])]),
                html.Tr([html.Td(['α', html.Sub('x')]), html.Td(id='meqxh'),html.Td(['μ0H', html.Sub('eff,x'),'(T)']), html.Td(id='hefxh'),html.Td(['μ0H', html.Sub('a,x'),'(T)']), html.Td(id='haxh')]),
                html.Tr([html.Td(['α', html.Sub('y')]), html.Td(id='meqyh'),html.Td(['μ0H', html.Sub('eff,y'),'(T)']), html.Td(id='hefyh'),html.Td(['μ0H', html.Sub('a,y'),'(T)']), html.Td(id='hayh')]),
                html.Tr([html.Td(['α', html.Sub('z')]), html.Td(id='meqzh'),html.Td(['μ0H', html.Sub('eff,z'),'(T)']), html.Td(id='hefzh'),html.Td(['μ0H', html.Sub('a,z'),'(T)']), html.Td(id='hazh')]),
                html.Tr([html.Td(['Torque |α x μ0Heff| (T)']), html.Td(id='torqueh')]),
                html.Tr([html.Td(['Status']), html.Td(id='stath')]),

            ]),
            html.Hr(),
            dcc.Markdown(''' **Generated figures to visualize magnetostriction:** '''),
            html.H6("The distance between a point on the surface and the origin (0,0,0) describes the simulated length lsim in the direction \u03B2=(sinθ*cosφ, sinθ*sinφ, cosθ), where θ and φ are the polar and azimuthal angles, respectively. The color of the surface corresponds to the relative length change multiplied by the scale factor in the measuring length direction \u03B2."),

            dcc.Graph(id='hex_3D'),
            dcc.Graph(id='hex_2D'),
            html.Hr(),
            dcc.Markdown(''' **Generated figure to visualize the magnetocrystalline anisotropy energy:**'''),
            html.H6("In this interactive plot you can visualize the magnetocrystalline anisotropy energy (Ea) as a function of the magnetization direction α=(sinθ*cosφ, sinθ*sinφ, cosθ)."),

            dcc.Graph(id='mae_hex_3D'),
            html.Hr(),


            html.H4("Magnetic properties of some hexagonal crystals"),
            html.Table([
                html.Tr([html.Td('Material'), html.Td('Space group'),html.Td('Temperature (K)'),html.Td(html.Div(['    \u03BB',html.Sup('\u03B11,0'),'    '])), html.Td(html.Div(['    \u03BB',html.Sup('\u03B12,0'),'    '])), html.Td(html.Div(['\u03BB',html.Sup('\u03B11,2')])), html.Td(html.Div(['\u03BB',html.Sup('\u03B12,2')])), html.Td(html.Div(['\u03BB',html.Sup('\u03B3,2')])),html.Td(html.Div(['\u03BB',html.Sup('\u03B5,2')])),html.Td(html.Div(['K',html.Sub("1"),'(MJ/m',html.Sub("3"),')'])),html.Td(html.Div(['K',html.Sub("2"),'(MJ/m',html.Sub("3"),')'])),html.Td(html.Div(['μ0Ms(Tesla)']))]),
                html.Tr([html.Td('HCP Co'), html.Td('194'),html.Td('0'),html.Td('    -   '),html.Td('    -   '), html.Td('0.000095'), html.Td('-0.000126'),html.Td('0.000057'),html.Td('-0.000286'),html.Td('0.7'),html.Td('0.18'),html.Td('1.81')]),
                html.Tr([html.Td('HCP Gd'), html.Td('194'),html.Td('0'),html.Td('    -   '),html.Td('    -   '), html.Td('0.00014'), html.Td('-0.00013'),html.Td('0.00011'),html.Td('0.00002'),html.Td('-0.12'),html.Td('0.08'),html.Td('2.59')]),
                html.Tr([html.Td('HCP Tb'), html.Td('194'),html.Td('0'),html.Td('    -   '),html.Td('    -   '), html.Td('-0.0026'), html.Td('0.009'),html.Td('0.0087'),html.Td('0.015'),html.Td(' -56.5 '),html.Td(' -4.6 '),html.Td(' - ')]),
            ])

        ])


######### Hex I (Mason)


    elif system == 'hexm':

        return html.Div(children=[

            html.H3(" Single crystal: Hexagonal I (space group numbers  177-194) - Mason definition"),
            html.H4("Theory"),
            html.H6("In Mason definition, the relative length change for hexagonal (I) systems is given by:"),
            html.Div([html.Img(src=app.get_asset_url('eq_hexm.png'))]),
            html.H6("where αi and βi (i=x,y,z) are the direction of magnetization and the measuring length direction, respectively. Magentostriction is a small effect that is hard to visualize. To facilitate its visualization in the simulation, we multiply the right hand side of this equation by a scale factor parameter which can be modified by the user."),
            html.H6("These magnetostrictive coefficients are related to those defined by Clark through the following equations:"),
            html.Div([html.Img(src=app.get_asset_url('eq_hexm_conv.png'))]),
            html.H6("The magnetocrystalline anisotropy energy for hexagonal systems is"),
            html.Div([html.Img(src=app.get_asset_url('mae_hex.png'))]),
            html.H6("where K0, K1 and K2 are the magnetocrystalline anisotropy constants. K0 can be used to shift the minimum energy reference in the 3D visualization of the magnetocrystalline anisotropy energy. The magnetocrystalline anisotropy field Ha is"),
            html.Div([html.Img(src=app.get_asset_url('mae_field.png'))]),
            html.H6("where μ0 is the vacuum permeability and Ms is the saturation magnetization. In the simulation, it is assumed that the magnetization is saturated along the effective field Heff"),
            html.Div([html.Img(src=app.get_asset_url('eff_field.png'))]),
            html.H6("where H is the external magnetic field. The direction of the equilibrium magnetization α is calculated by the Landau-Lifshitz-Gilbert equation"),
            html.Div([html.Img(src=app.get_asset_url('llg.png'))]),
            html.H6("where γ is the electron gyromagnetic ratio, and η is the damping parameter. The torque |α x μ0Heff| is used as the criterion for the numerical convergence of the equilibrium state, it is recommended to use a tolerance for the torque lower than 0.0001 Tesla. Default values of the damping parameter, time step (dt) and total number of iteration steps usually work fine, but the user can change them in case the tolerance for the torque is not achieved."),
            html.Div([html.Img(src=app.get_asset_url('heffh.png'))]),


            html.H4("Parameters of the simulation"),
            html.H6("(Press Enter after changing any input to update the figures)"),
            html.Hr(),
            dcc.Markdown(''' **External magnetic field:**'''),
            html.Div(['μ0H',html.Sub('x'),'(Tesla)'" = ",
              dcc.Input(id='hfieldxm', value=5.000, type='number', debounce=True, step=0.0000001)]),
            html.Div(['μ0H',html.Sub('y'),'(Tesla)'" = ",
              dcc.Input(id='hfieldym', value=0.500, type='number', debounce=True, step=0.0000001)]),
            html.Div(['μ0H',html.Sub('z'),'(Tesla)'" = ",
              dcc.Input(id='hfieldzm', value=0.005, type='number', debounce=True, step=0.0000001)]),
            dcc.Markdown(''' **Magnetostrictive coefficients:**'''),
            html.Div(['\u03BB',html.Sup("\u03B11,0"),html.Sub("Mason")," = ",
              dcc.Input(id='hL01m', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sup("\u03B12,0"),html.Sub("Mason")," = ",
              dcc.Input(id='hL02m', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sub("A")," = ",
              dcc.Input(id='hLa', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sub("B")," = ",
              dcc.Input(id='hLb', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sub("C")," = ",
              dcc.Input(id='hLc', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sub("D")," = ",
              dcc.Input(id='hLd', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(["Scale factor = ",
              dcc.Input(id='hscalem', value=1.0, type='number', debounce=True)]),
            dcc.Markdown(''' **Magnetocrystalline anisotropy constants:**'''),
            html.Div(['K',html.Sub("0"), "(MJ/m",html.Sup("3"),")"," = ",
              dcc.Input(id='k0hm', value=0.001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['K',html.Sub("1"), "(MJ/m",html.Sup("3"),")"," = ",
              dcc.Input(id='k1hm', value=0.001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['K',html.Sub("2"),"(MJ/m",html.Sup("3"),")"," = ",
              dcc.Input(id='k2hm', value=0.001, type='number', debounce=True, step=0.0000001)]),
            dcc.Markdown(''' **Saturation magnetization:**'''),
            html.Div(['μ0Ms (Tesla)'," = ",
              dcc.Input(id='mshm', value=0.01, type='number', debounce=True, step=0.0001)]),
            dcc.Markdown(''' **Landau-Lifshitz-Gilbert solver to find the direction of the equilibrium magnetization (α):**'''),
            html.Div(['Damping LLG '," = ",
              dcc.Input(id='alpha0hm', value=0.95, type='number', debounce=True, step=0.00001)]),
            html.Div(['Torque tolerance |α x μ0Heff| (Tesla)'," = ",
              dcc.Input(id='tolhm', value=0.00001, type='number', debounce=True, step=0.000000001)]),
            html.Div(['Time step LLG (ps)'," = ",
              dcc.Input(id='dthm', value=0.01, type='number', debounce=True, step=0.00001)]),
            html.Div(['Total number of iteration steps'," = ",
              dcc.Input(id='nthm', value=500000, type='number', debounce=True, step=1)]),

            html.Div(id='my-output-hm'),
            html.Hr(),

           html.H4("Simulation"),

            dcc.Markdown(''' **LLG solver output:**'''),
            html.Table([
                html.Tr([html.Td(['Parameters']), html.Td(['Calculated values']),html.Td(['Parameters']), html.Td(['Calculated values']),html.Td(['Parameters']), html.Td(['Calculated values'])]),
                html.Tr([html.Td(['α', html.Sub('x')]), html.Td(id='meqxhm'),html.Td(['μ0H', html.Sub('eff,x'),'(T)']), html.Td(id='hefxhm'),html.Td(['μ0H', html.Sub('a,x'),'(T)']), html.Td(id='haxhm')]),
                html.Tr([html.Td(['α', html.Sub('y')]), html.Td(id='meqyhm'),html.Td(['μ0H', html.Sub('eff,y'),'(T)']), html.Td(id='hefyhm'),html.Td(['μ0H', html.Sub('a,y'),'(T)']), html.Td(id='hayhm')]),
                html.Tr([html.Td(['α', html.Sub('z')]), html.Td(id='meqzhm'),html.Td(['μ0H', html.Sub('eff,z'),'(T)']), html.Td(id='hefzhm'),html.Td(['μ0H', html.Sub('a,z'),'(T)']), html.Td(id='hazhm')]),
                html.Tr([html.Td(['Torque |α x μ0Heff| (T)']), html.Td(id='torquehm')]),
                html.Tr([html.Td(['Status']), html.Td(id='stathm')]),

            ]),
            html.Hr(),
            dcc.Markdown(''' **Generated figures to visualize magnetostriction:** '''),
            html.H6("The distance between a point on the surface and the origin (0,0,0) describes the simulated length lsim in the direction \u03B2=(sinθ*cosφ, sinθ*sinφ, cosθ), where θ and φ are the polar and azimuthal angles, respectively. The color of the surface corresponds to the relative length change multiplied by the scale factor in the measuring length direction \u03B2."),

            dcc.Graph(id='hex_3Dm'),
            dcc.Graph(id='hex_2Dm'),
            html.Hr(),
            dcc.Markdown(''' **Generated figure to visualize the magnetocrystalline anisotropy energy:**'''),
            html.H6("In this interactive plot you can visualize the magnetocrystalline anisotropy energy (Ea) as a function of the magnetization direction α=(sinθ*cosφ, sinθ*sinφ, cosθ)."),

            dcc.Graph(id='mae_hex_3Dm'),
            html.Hr(),


            html.H4("Magnetic properties of some hexagonal crystals"),
            html.Table([
                html.Tr([html.Td('Material'), html.Td('Space group'),html.Td('Temperature (K)'),html.Td(html.Div(['    \u03BB',html.Sup('\u03B11,0'),html.Sub('Mason'),'    '])), html.Td(html.Div(['    \u03BB',html.Sup('\u03B12,0'),html.Sub('Mason'),'    '])), html.Td(html.Div(['\u03BB',html.Sub('A')])), html.Td(html.Div(['\u03BB',html.Sub('B')])), html.Td(html.Div(['\u03BB',html.Sub('C')])),html.Td(html.Div(['\u03BB',html.Sub('D')])),html.Td(html.Div(['K',html.Sub("1"),'(MJ/m',html.Sub("3"),')'])),html.Td(html.Div(['K',html.Sub("2"),'(MJ/m',html.Sub("3"),')'])),html.Td(html.Div(['μ0Ms(Tesla)']))]),
                html.Tr([html.Td('HCP Co'), html.Td('194'),html.Td('0'),html.Td('    -   '),html.Td('    -   '), html.Td('-0.000066'), html.Td('-0.000123'),html.Td('0.000126'),html.Td('-0.000128'),html.Td('0.7'),html.Td('0.18'),html.Td('1.81')]),
            ])

        ])



######## Trig I

    elif system == 'trig':

        return html.Div(children=[

            html.H3("Single crystal: Trigonal I (space group numbers  149-167)"),
            html.H4("Theory"),
            html.H6("The relative length change for trigonal (I) systems is given by:"),
            html.Div([html.Img(src=app.get_asset_url('eq_trig.png'))]),
            html.H6("where αi and βi (i=x,y,z) are the direction of magnetization and the measuring length direction, respectively. Magentostriction is a small effect that is hard to visualize. To facilitate its visualization in the simulation, we multiply the right hand side of this equation by a scale factor parameter which can be modified by the user."),

            html.H6("The magnetocrystalline anisotropy energy for trigonal systems is"),
            html.Div([html.Img(src=app.get_asset_url('mae_hex.png'))]),
            html.H6("where K0, K1 and K2 are the magnetocrystalline anisotropy constants. K0 can be used to shift the minimum energy reference in the 3D visualization of the magnetocrystalline anisotropy energy. The magnetocrystalline anisotropy field Ha is"),
            html.Div([html.Img(src=app.get_asset_url('mae_field.png'))]),
            html.H6("where μ0 is the vacuum permeability and Ms is the saturation magnetization. In the simulation, it is assumed that the magnetization is saturated along the effective field Heff"),
            html.Div([html.Img(src=app.get_asset_url('eff_field.png'))]),
            html.H6("where H is the external magnetic field. The direction of the equilibrium magnetization α is calculated by the Landau-Lifshitz-Gilbert equation"),
            html.Div([html.Img(src=app.get_asset_url('llg.png'))]),
            html.H6("where γ is the electron gyromagnetic ratio, and η is the damping parameter. The torque |α x μ0Heff| is used as the criterion for the numerical convergence of the equilibrium state, it is recommended to use a tolerance for the torque lower than 0.0001 Tesla. Default values of the damping parameter, time step (dt) and total number of iteration steps usually work fine, but the user can change them in case the tolerance for the torque is not achieved."),
            html.Div([html.Img(src=app.get_asset_url('hefftr.png'))]),


            html.H4("Parameters of the simulation"),
            html.H6("(Press Enter after changing any input to update the figures)"),
            html.Hr(),

            dcc.Markdown(''' **External magnetic field:**'''),
            html.Div(['μ0H',html.Sub('x'),'(Tesla)'" = ",
              dcc.Input(id='trfieldx', value=5.000, type='number', debounce=True, step=0.0000001)]),
            html.Div(['μ0H',html.Sub('y'),'(Tesla)'" = ",
              dcc.Input(id='trfieldy', value=0.500, type='number', debounce=True, step=0.0000001)]),
            html.Div(['μ0H',html.Sub('z'),'(Tesla)'" = ",
              dcc.Input(id='trfieldz', value=0.005, type='number', debounce=True, step=0.0000001)]),
            dcc.Markdown(''' **Magnetostrictive coefficients:**'''),
            html.Div(['\u03BB',html.Sup("\u03B11,0")," = ",
              dcc.Input(id='trL01', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sup("\u03B12,0")," = ",
              dcc.Input(id='trL02', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sup("\u03B11,2")," = ",
              dcc.Input(id='trLa1', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sup("\u03B12,2")," = ",
              dcc.Input(id='trLa2', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sup("\u03B3,1")," = ",
              dcc.Input(id='trLg1', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sup("\u03B3,2")," = ",
              dcc.Input(id='trLg2', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sub("12")," = ",
              dcc.Input(id='trL12', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sub("21")," = ",
              dcc.Input(id='trL21', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(["Scale factor = ",
              dcc.Input(id='trscale', value=1.0, type='number', debounce=True)]),
            dcc.Markdown(''' **Magnetocrystalline anisotropy constants:**'''),
            html.Div(['K',html.Sub("0"), "(MJ/m",html.Sup("3"),")"," = ",
              dcc.Input(id='k0tr', value=0.001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['K',html.Sub("1"), "(MJ/m",html.Sup("3"),")"," = ",
              dcc.Input(id='k1tr', value=0.001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['K',html.Sub("2"),"(MJ/m",html.Sup("3"),")"," = ",
              dcc.Input(id='k2tr', value=0.001, type='number', debounce=True, step=0.0000001)]),
            dcc.Markdown(''' **Saturation magnetization:**'''),
            html.Div(['μ0Ms (Tesla)'," = ",
              dcc.Input(id='mstr', value=0.01, type='number', debounce=True, step=0.0001)]),
            dcc.Markdown(''' **Landau-Lifshitz-Gilbert solver to find the direction of the equilibrium magnetization (α):**'''),
            html.Div(['Damping LLG '," = ",
              dcc.Input(id='alpha0tr', value=0.95, type='number', debounce=True, step=0.00001)]),
            html.Div(['Torque tolerance |α x μ0Heff| (Tesla)'," = ",
              dcc.Input(id='toltr', value=0.00001, type='number', debounce=True, step=0.000000001)]),
            html.Div(['Time step LLG (ps)'," = ",
              dcc.Input(id='dttr', value=0.01, type='number', debounce=True, step=0.00001)]),
            html.Div(['Total number of iteration steps'," = ",
              dcc.Input(id='nttr', value=500000, type='number', debounce=True, step=1)]),


            html.Div(id='my-output-tr'),
            html.Hr(),
            html.H4("Simulation"),
            dcc.Markdown(''' **LLG solver output:**'''),
            html.Table([
                html.Tr([html.Td(['Parameters']), html.Td(['Calculated values']),html.Td(['Parameters']), html.Td(['Calculated values']),html.Td(['Parameters']), html.Td(['Calculated values'])]),
                html.Tr([html.Td(['α', html.Sub('x')]), html.Td(id='meqxtr'),html.Td(['μ0H', html.Sub('eff,x'),'(T)']), html.Td(id='hefxtr'),html.Td(['μ0H', html.Sub('a,x'),'(T)']), html.Td(id='haxtr')]),
                html.Tr([html.Td(['α', html.Sub('y')]), html.Td(id='meqytr'),html.Td(['μ0H', html.Sub('eff,y'),'(T)']), html.Td(id='hefytr'),html.Td(['μ0H', html.Sub('a,y'),'(T)']), html.Td(id='haytr')]),
                html.Tr([html.Td(['α', html.Sub('z')]), html.Td(id='meqztr'),html.Td(['μ0H', html.Sub('eff,z'),'(T)']), html.Td(id='hefztr'),html.Td(['μ0H', html.Sub('a,z'),'(T)']), html.Td(id='haztr')]),
                html.Tr([html.Td(['Torque |α x μ0Heff| (T)']), html.Td(id='torquetr')]),
                html.Tr([html.Td(['Status']), html.Td(id='stattr')]),
            ]),
            html.Hr(),
            dcc.Markdown(''' **Generated figures to visualize magnetostriction:** '''),
            html.H6("The distance between a point on the surface and the origin (0,0,0) describes the simulated length lsim in the direction \u03B2=(sinθ*cosφ, sinθ*sinφ, cosθ), where θ and φ are the polar and azimuthal angles, respectively. The color of the surface corresponds to the relative length change multiplied by the scale factor in the measuring length direction \u03B2."),

            dcc.Graph(id='tri_3D'),
            dcc.Graph(id='tri_2D'),
            html.Hr(),
            dcc.Markdown(''' **Generated figure to visualize the magnetocrystalline anisotropy energy:**'''),
            html.H6("In this interactive plot you can visualize the magnetocrystalline anisotropy energy (Ea) as a function of the magnetization direction α=(sinθ*cosφ, sinθ*sinφ, cosθ)."),

            dcc.Graph(id='mae_tri_3D'),

        ])

######## Tet I

    elif system == 'tet':

        return html.Div(children=[

            html.H3("Single crystal: Tetragonal I (space group numbers  89-142)"),
            html.H4("Theory"),
            html.H6("The relative length change for tetragoanl (I) systems is given by:"),
            html.Div([html.Img(src=app.get_asset_url('eq_tet.png'))]),
            html.H6("where αi and βi (i=x,y,z) are the direction of magnetization and the measuring length direction, respectively. Magentostriction is a small effect that is hard to visualize. To facilitate its visualization in the simulation, we multiply the right hand side of this equation by a scale factor parameter which can be modified by the user."),

            html.H6("The magnetocrystalline anisotropy energy for tetragonal systems is"),
            html.Div([html.Img(src=app.get_asset_url('mae_hex.png'))]),
            html.H6("where K0, K1 and K2 are the magnetocrystalline anisotropy constants. K0 can be used to shift the minimum energy reference in the 3D visualization of the magnetocrystalline anisotropy energy. The magnetocrystalline anisotropy field Ha is"),
            html.Div([html.Img(src=app.get_asset_url('mae_field.png'))]),
            html.H6("where μ0 is the vacuum permeability and Ms is the saturation magnetization. In the simulation, it is assumed that the magnetization is saturated along the effective field Heff"),
            html.Div([html.Img(src=app.get_asset_url('eff_field.png'))]),
            html.H6("where H is the external magnetic field. The direction of the equilibrium magnetization α is calculated by the Landau-Lifshitz-Gilbert equation"),
            html.Div([html.Img(src=app.get_asset_url('llg.png'))]),
            html.H6("where γ is the electron gyromagnetic ratio, and η is the damping parameter. The torque |α x μ0Heff| is used as the criterion for the numerical convergence of the equilibrium state, it is recommended to use a tolerance for the torque lower than 0.0001 Tesla. Default values of the damping parameter, time step (dt) and total number of iteration steps usually work fine, but the user can change them in case the tolerance for the torque is not achieved."),
            html.Div([html.Img(src=app.get_asset_url('heffte.png'))]),


            html.H4("Parameters of the simulation"),
            html.H6("(Press Enter after changing any input to update the figures)"),
            html.Hr(),

            dcc.Markdown(''' **External magnetic field:**'''),
            html.Div(['μ0H',html.Sub('x'),'(Tesla)'" = ",
              dcc.Input(id='tefieldx', value=5.000, type='number', debounce=True, step=0.0000001)]),
            html.Div(['μ0H',html.Sub('y'),'(Tesla)'" = ",
              dcc.Input(id='tefieldy', value=0.500, type='number', debounce=True, step=0.0000001)]),
            html.Div(['μ0H',html.Sub('z'),'(Tesla)'" = ",
              dcc.Input(id='tefieldz', value=0.005, type='number', debounce=True, step=0.0000001)]),
            dcc.Markdown(''' **Magnetostrictive coefficients:**'''),
            html.Div(['\u03BB',html.Sup("\u03B11,0")," = ",
              dcc.Input(id='teL01', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sup("\u03B12,0")," = ",
              dcc.Input(id='teL02', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sup("\u03B11,2")," = ",
              dcc.Input(id='teLa1', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sup("\u03B12,2")," = ",
              dcc.Input(id='teLa2', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sup("\u03B3,2")," = ",
              dcc.Input(id='teLg', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sup("\u03B4,2")," = ",
              dcc.Input(id='teLd', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sup("\u03B5,2")," = ",
              dcc.Input(id='teLe', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(["Scale factor = ",
              dcc.Input(id='tescale', value=1.0, type='number', debounce=True)]),
            dcc.Markdown(''' **Magnetocrystalline anisotropy constants:**'''),
            html.Div(['K',html.Sub("0"), "(MJ/m",html.Sup("3"),")"," = ",
              dcc.Input(id='k0te', value=0.001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['K',html.Sub("1"), "(MJ/m",html.Sup("3"),")"," = ",
              dcc.Input(id='k1te', value=0.001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['K',html.Sub("2"),"(MJ/m",html.Sup("3"),")"," = ",
              dcc.Input(id='k2te', value=0.001, type='number', debounce=True, step=0.0000001)]),
            dcc.Markdown(''' **Saturation magnetization:**'''),
            html.Div(['μ0Ms (Tesla)'," = ",
              dcc.Input(id='mste', value=0.01, type='number', debounce=True, step=0.0001)]),
            dcc.Markdown(''' **Landau-Lifshitz-Gilbert solver to find the direction of the equilibrium magnetization (α):**'''),
            html.Div(['Damping LLG '," = ",
              dcc.Input(id='alpha0te', value=0.95, type='number', debounce=True, step=0.00001)]),
            html.Div(['Torque tolerance |α x μ0Heff| (Tesla)'," = ",
              dcc.Input(id='tolte', value=0.00001, type='number', debounce=True, step=0.000000001)]),
            html.Div(['Time step LLG (ps)'," = ",
              dcc.Input(id='dtte', value=0.01, type='number', debounce=True, step=0.00001)]),
            html.Div(['Total number of iteration steps'," = ",
              dcc.Input(id='ntte', value=500000, type='number', debounce=True, step=1)]),


            html.Div(id='my-output-te'),
            html.Hr(),
            html.H4("Simulation"),
            dcc.Markdown(''' **LLG solver output:**'''),
            html.Table([
                html.Tr([html.Td(['Parameters']), html.Td(['Calculated values']),html.Td(['Parameters']), html.Td(['Calculated values']),html.Td(['Parameters']), html.Td(['Calculated values'])]),
                html.Tr([html.Td(['α', html.Sub('x')]), html.Td(id='meqxte'),html.Td(['μ0H', html.Sub('eff,x'),'(T)']), html.Td(id='hefxte'),html.Td(['μ0H', html.Sub('a,x'),'(T)']), html.Td(id='haxte')]),
                html.Tr([html.Td(['α', html.Sub('y')]), html.Td(id='meqyte'),html.Td(['μ0H', html.Sub('eff,y'),'(T)']), html.Td(id='hefyte'),html.Td(['μ0H', html.Sub('a,y'),'(T)']), html.Td(id='hayte')]),
                html.Tr([html.Td(['α', html.Sub('z')]), html.Td(id='meqzte'),html.Td(['μ0H', html.Sub('eff,z'),'(T)']), html.Td(id='hefzte'),html.Td(['μ0H', html.Sub('a,z'),'(T)']), html.Td(id='hazte')]),
                html.Tr([html.Td(['Torque |α x μ0Heff| (T)']), html.Td(id='torquete')]),
                html.Tr([html.Td(['Status']), html.Td(id='statte')]),
            ]),
            html.Hr(),
            dcc.Markdown(''' **Generated figures to visualize magnetostriction:** '''),
            html.H6("The distance between a point on the surface and the origin (0,0,0) describes the simulated length lsim in the direction \u03B2=(sinθ*cosφ, sinθ*sinφ, cosθ), where θ and φ are the polar and azimuthal angles, respectively. The color of the surface corresponds to the relative length change multiplied by the scale factor in the measuring length direction \u03B2."),
            dcc.Graph(id='tet_3D'),
            dcc.Graph(id='tet_2D'),
            html.Hr(),
            dcc.Markdown(''' **Generated figure to visualize the magnetocrystalline anisotropy energy:**'''),
            html.H6("In this interactive plot you can visualize the magnetocrystalline anisotropy energy (Ea) as a function of the magnetization direction α=(sinθ*cosφ, sinθ*sinφ, cosθ)."),

            dcc.Graph(id='mae_tet_3D'),
        ])


######## Orth

    elif system == 'ort':

        return html.Div(children=[

            html.H3("Single crystal: Orthorhombic (space group numbers  16-74)"),
            html.H4("Theory"),
            html.H6("This interactive applet shows the magnetostriction due to Joule effect for some crystal systems."),
            html.Div([html.Img(src=app.get_asset_url('eq_ort.png'))]),
            html.H6("where αi and βi (i=x,y,z) are the direction of magnetization and the measuring length direction, respectively. Magentostriction is a small effect that is hard to visualize. To facilitate its visualization in the simulation, we multiply the right hand side of this equation by a scale factor parameter which can be modified by the user."),

            html.H6("The magnetocrystalline anisotropy energy for orthorhombic systems is"),
            html.Div([html.Img(src=app.get_asset_url('mae_ort.png'))]),
            html.H6("where K0, K1 and K2 are the magnetocrystalline anisotropy constants. K0 can be used to shift the minimum energy reference in the 3D visualization of the magnetocrystalline anisotropy energy. The magnetocrystalline anisotropy field Ha is"),
            html.Div([html.Img(src=app.get_asset_url('mae_field.png'))]),
            html.H6("where μ0 is the vacuum permeability and Ms is the saturation magnetization. In the simulation, it is assumed that the magnetization is saturated along the effective field Heff"),
            html.Div([html.Img(src=app.get_asset_url('eff_field.png'))]),
            html.H6("where H is the external magnetic field. The direction of the equilibrium magnetization α is calculated by the Landau-Lifshitz-Gilbert equation"),
            html.Div([html.Img(src=app.get_asset_url('llg.png'))]),
            html.H6("where γ is the electron gyromagnetic ratio, and η is the damping parameter. The torque |α x μ0Heff| is used as the criterion for the numerical convergence of the equilibrium state, it is recommended to use a tolerance for the torque lower than 0.0001 Tesla. Default values of the damping parameter, time step (dt) and total number of iteration steps usually work fine, but the user can change them in case the tolerance for the torque is not achieved."),
            html.Div([html.Img(src=app.get_asset_url('heffo.png'))]),

            html.H4("Parameters of the simulation"),
            html.H6("(Press Enter after changing any input to update the figures)"),
            html.Hr(),

            dcc.Markdown(''' **External magnetic field:**'''),
            html.Div(['μ0H',html.Sub('x'),'(Tesla)'" = ",
              dcc.Input(id='ofieldx', value=5.000, type='number', debounce=True, step=0.0000001)]),
            html.Div(['μ0H',html.Sub('y'),'(Tesla)'" = ",
              dcc.Input(id='ofieldy', value=0.500, type='number', debounce=True, step=0.0000001)]),
            html.Div(['μ0H',html.Sub('z'),'(Tesla)'" = ",
              dcc.Input(id='ofieldz', value=0.005, type='number', debounce=True, step=0.0000001)]),
            dcc.Markdown(''' **Magnetostrictive coefficients:**'''),
            html.Div(['\u03BB',html.Sup("\u03B11,0")," = ",
              dcc.Input(id='oL01', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sup("\u03B12,0")," = ",
              dcc.Input(id='oL02', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sup("\u03B13,0")," = ",
              dcc.Input(id='oL03', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sub("1")," = ",
              dcc.Input(id='oL1', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sub("2")," = ",
              dcc.Input(id='oL2', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sub("3")," = ",
              dcc.Input(id='oL3', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sub("4")," = ",
              dcc.Input(id='oL4', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sub("5")," = ",
              dcc.Input(id='oL5', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sub("6")," = ",
              dcc.Input(id='oL6', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sub("7")," = ",
              dcc.Input(id='oL7', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sub("8")," = ",
              dcc.Input(id='oL8', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['\u03BB',html.Sub("9")," = ",
              dcc.Input(id='oL9', value=0.0000001, type='number', debounce=True, step=0.0000001)]),
            html.Div(["Scale factor = ",
              dcc.Input(id='oscale', value=1.0, type='number', debounce=True)]),
            dcc.Markdown(''' **Magnetocrystalline anisotropy constants:**'''),
            html.Div(['K',html.Sub("0"), "(MJ/m",html.Sup("3"),")"," = ",
              dcc.Input(id='k0o', value=0.001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['K',html.Sub("1"), "(MJ/m",html.Sup("3"),")"," = ",
              dcc.Input(id='k1o', value=0.001, type='number', debounce=True, step=0.0000001)]),
            html.Div(['K',html.Sub("2"),"(MJ/m",html.Sup("3"),")"," = ",
              dcc.Input(id='k2o', value=0.001, type='number', debounce=True, step=0.0000001)]),
            dcc.Markdown(''' **Saturation magnetization:**'''),
            html.Div(['μ0Ms (Tesla)'," = ",
              dcc.Input(id='mso', value=0.01, type='number', debounce=True, step=0.0001)]),
            dcc.Markdown(''' **Landau-Lifshitz-Gilbert solver to find the direction of the equilibrium magnetization (α):**'''),
            html.Div(['Damping LLG '," = ",
              dcc.Input(id='alpha0o', value=0.95, type='number', debounce=True, step=0.00001)]),
            html.Div(['Torque tolerance |α x μ0Heff| (Tesla)'," = ",
              dcc.Input(id='tolo', value=0.00001, type='number', debounce=True, step=0.000000001)]),
            html.Div(['Time step LLG (ps)'," = ",
              dcc.Input(id='dto', value=0.01, type='number', debounce=True, step=0.00001)]),
            html.Div(['Total number of iteration steps'," = ",
              dcc.Input(id='nto', value=500000, type='number', debounce=True, step=1)]),

            html.Div(id='my-output-o'),
            html.Hr(),
            html.H4("Simulation"),
            dcc.Markdown(''' **LLG solver output:**'''),
            html.Table([
                html.Tr([html.Td(['Parameters']), html.Td(['Calculated values']),html.Td(['Parameters']), html.Td(['Calculated values']),html.Td(['Parameters']), html.Td(['Calculated values'])]),
                html.Tr([html.Td(['α', html.Sub('x')]), html.Td(id='meqxo'),html.Td(['μ0H', html.Sub('eff,x'),'(T)']), html.Td(id='hefxo'),html.Td(['μ0H', html.Sub('a,x'),'(T)']), html.Td(id='haxo')]),
                html.Tr([html.Td(['α', html.Sub('y')]), html.Td(id='meqyo'),html.Td(['μ0H', html.Sub('eff,y'),'(T)']), html.Td(id='hefyo'),html.Td(['μ0H', html.Sub('a,y'),'(T)']), html.Td(id='hayo')]),
                html.Tr([html.Td(['α', html.Sub('z')]), html.Td(id='meqzo'),html.Td(['μ0H', html.Sub('eff,z'),'(T)']), html.Td(id='hefzo'),html.Td(['μ0H', html.Sub('a,z'),'(T)']), html.Td(id='hazo')]),
                html.Tr([html.Td(['Torque |α x μ0Heff| (T)']), html.Td(id='torqueo')]),
                html.Tr([html.Td(['Status']), html.Td(id='stato')]),
            ]),
            html.Hr(),
            dcc.Markdown(''' **Generated figures to visualize magnetostriction:** '''),
            html.H6("The distance between a point on the surface and the origin (0,0,0) describes the simulated length lsim in the direction \u03B2=(sinθ*cosφ, sinθ*sinφ, cosθ), where θ and φ are the polar and azimuthal angles, respectively. The color of the surface corresponds to the relative length change multiplied by the scale factor in the measuring length direction \u03B2."),
            dcc.Graph(id='ort_3D'),
            dcc.Graph(id='ort_2D'),
            html.Hr(),
            dcc.Markdown(''' **Generated figure to visualize the magnetocrystalline anisotropy energy:**'''),
            html.H6("In this interactive plot you can visualize the magnetocrystalline anisotropy energy (Ea) as a function of the magnetization direction α=(sinθ*cosφ, sinθ*sinφ, cosθ)."),

            dcc.Graph(id='mae_ort_3D'),
        ])



    else:

        return 'No crystal system selected'


def energy(crystal,sx,sy,sz,ms0,efx,efy,efz,k01,k02):

    mu0=4.0*np.pi*10**(-7)
    ene_ef=-(ms0/mu0)*(sx*efx+sy*efy+sz*efz)
    k01=k01*10.0**6
    k02=k02*10.0**6

    if crystal=='cub':
        ene_af=k01*(sx**2.0*sy**2.0+sx**2.0*sz**2.0+sy**2.0*sz**2.0)+k02*sx**2.0*sy**2.0*sz**2.0
    elif crystal=='uni':
        ene_af=k01*(1.0-sz**2.0)+k02*(1.0-sz**2.0)**2.0
    elif crystal=='ort':
        ene_af=k01*sx**2.0+k02*sy**2.0

    ene_tot=ene_ef+ene_af
    return ene_tot

def field(crystal,sx,sy,sz,ms0,efx,efy,efz,k01,k02):

    mu0=4.0*np.pi*10**(-7)

    k01=k01*10.0**6
    k02=k02*10.0**6

    if crystal=='cub':
        hax=-(2.0*k01*mu0/ms0)*sx*(sy**2.0+sz**2.0)-(2.0*k02*mu0/ms0)*sx*(sy**2.0)*(sz**2.0)
        hay=-(2.0*k01*mu0/ms0)*sy*(sx**2.0+sz**2.0)-(2.0*k02*mu0/ms0)*sy*(sx**2.0)*(sz**2.0)
        haz=-(2.0*k01*mu0/ms0)*sz*(sy**2.0+sx**2.0)-(2.0*k02*mu0/ms0)*sz*(sy**2.0)*(sx**2.0)
    elif crystal=='uni':
       # hax=-(2.0*k01*mu0/ms0)*sx-(4.0*k02*mu0/ms0)*sx*(sx**2.0+sy**2.0)
       # hay=-(2.0*k01*mu0/ms0)*sy-(4.0*k02*mu0/ms0)*sy*(sx**2.0+sy**2.0)
        hax=0.0
        hay=0.0
        haz=(2.0*k01*mu0/ms0)*sz+(4.0*k02*mu0/ms0)*sz*(1.0-sz**2.0)
    elif crystal=='ort':
        hax=-(2.0*k01*mu0/ms0)*sx
        hay=-(2.0*k02*mu0/ms0)*sy
        haz=0

    fieldx=hax+efx
    fieldy=hay+efy
    fieldz=haz+efz

    return fieldx,fieldy,fieldz

def llg(crystal0,mms0,effx,effy,effz,kk01,kk02,alpha,tol0,dt,ntot):

    hmod=np.sqrt(effx**2.0+effy**2.0+effz**2.0)

    if hmod < 10.0**(-6):
        uu=np.random.random_sample()*np.pi
        vv=np.random.random_sample()*2.0*np.pi

        a0x = np.sin(uu)*np.cos(vv)
        a0y = np.sin(uu)*np.sin(vv)
        a0z = np.cos(uu)

    else:

        a0x = effx/hmod
        a0y = effy/hmod
        a0z = effz/hmod

        a0x=a0x/np.sqrt(a0x**2+a0y**2+a0z**2.0)
        a0y=a0y/np.sqrt(a0x**2+a0y**2+a0z**2.0)
        a0z=a0z/np.sqrt(a0x**2+a0y**2+a0z**2.0)

        a0x=a0x+np.random.uniform(-1,1)*10.0**(-3)
        a0y=a0y+np.random.uniform(-1,1)*10.0**(-3)
        a0z=a0z+np.random.uniform(-1,1)*10.0**(-3)

        a0x=a0x/np.sqrt(a0x**2+a0y**2+a0z**2.0)
        a0y=a0y/np.sqrt(a0x**2+a0y**2+a0z**2.0)
        a0z=a0z/np.sqrt(a0x**2+a0y**2+a0z**2.0)


    dt=dt*10.0**(-12)

    mu0=4.0*np.pi*10**(-7)
    gamma=1.760859*10**11

    torque=1.0

    count=0

    while count < ntot:

        heffx,heffy,heffz=field(crystal0,a0x,a0y,a0z,mms0,effx,effy,effz,kk01,kk02)

        dax=((alpha*gamma)/(1.0+alpha**2.0))*(heffx*(a0y**2.0+a0z**2.0)-a0x*(heffy*a0y+heffz*a0z))
        day=((alpha*gamma)/(1.0+alpha**2.0))*(heffy*(a0x**2.0+a0z**2.0)-a0y*(heffx*a0x+heffz*a0z))
        daz=((alpha*gamma)/(1.0+alpha**2.0))*(heffz*(a0x**2.0+a0y**2.0)-a0z*(heffx*a0x+heffy*a0y))

        dax=dax+(gamma/(1.0+alpha**2.0))*(heffy*a0z-heffz*a0y)
        day=day+(gamma/(1.0+alpha**2.0))*(heffz*a0x-heffx*a0z)
        daz=daz+(gamma/(1.0+alpha**2.0))*(heffx*a0y-heffy*a0x)

        a0x=a0x+dt*dax
        a0y=a0y+dt*day
        a0z=a0z+dt*daz

        a0x=a0x/np.sqrt(a0x**2+a0y**2+a0z**2.0)
        a0y=a0y/np.sqrt(a0x**2+a0y**2+a0z**2.0)
        a0z=a0z/np.sqrt(a0x**2+a0y**2+a0z**2.0)

        torque=np.sqrt((heffy*a0z-heffz*a0y)**2.0+(heffz*a0x-heffx*a0z)**2.0+(heffx*a0y-heffy*a0x)**2.0)

        if torque < tol0:

            stat = 'The calculation converged successfully. The torque is lower than the tolerance.'
            return a0x,a0y,a0z,stat


        count=count+1

    stat = 'The calculation did not converge. The torque is greater than the tolerance. Try to increase the total number of iteration steps and press Enter.'

    return a0x,a0y,a0z,stat

def mc(crystal0,mms0,effx,effy,effz,kk01,kk02,tol0,ntot):

    mu0=4.0*np.pi*10**(-7)
    gamma=1.760859*10**11
    kb=1.3806503*10**(-23)
    muB=9.274009994*10**(-24)
    temp=1.0
    sigma=(2.0/25.0)*((kb*temp)/muB)**(1.0/5.0)

    torque=1.0
    count=0

    uu=np.random.random_sample()*np.pi
    vv=np.random.random_sample()*2.0*np.pi

    ssx = np.sin(uu)*np.cos(vv)
    ssy = np.sin(uu)*np.sin(vv)
    ssz = np.cos(uu)

    ene0=energy(crystal0,ssx,ssy,ssz,mms0,effx,effy,effz,kk01,kk02)

    while count < ntot:

            ssx1=ssx+sigma*np.random.normal()
            ssy1=ssy+sigma*np.random.normal()
            ssz1=ssz+sigma*np.random.normal()

            ssx1=ssx1/np.sqrt(ssx1**2.0+ssy1**2.0+ssz1**2.0)
            ssy1=ssy1/np.sqrt(ssx1**2.0+ssy1**2.0+ssz1**2.0)
            ssz1=ssz1/np.sqrt(ssx1**2.0+ssy1**2.0+ssz1**2.0)

            ene=energy(crystal0,ssx1,ssy1,ssz1,mms0,effx,effy,effz,kk01,kk02)

            dene=ene-ene0

            p=np.random.random_sample()

            pp=np.exp(-dene/(kb*temp))

            if p < pp:
                ene0=ene
                ssx=ssx1
                ssy=ssy1
                ssz=ssz1

                heffx,heffy,heffz=field(crystal0,ssx,ssy,ssz,mms0,effx,effy,effz,kk01,kk02)
                torque=np.sqrt((heffy*ssz-heffz*ssy)**2.0+(heffz*ssx-heffx*ssz)**2.0+(heffx*ssy-heffy*ssx)**2.0)

                if torque < tol0:
                    return ssx,ssy,ssz

            count=count+1

    return ssx,ssy,ssz



##############Cubic I -3D

@app.callback(
    [Output('cub_3D', 'figure'),
     Output('meqxc', 'children'),
     Output('meqyc', 'children'),
     Output('meqzc', 'children'),
     Output('hefxc', 'children'),
     Output('hefyc', 'children'),
     Output('hefzc', 'children'),
     Output('torquec', 'children'),
     Output('haxc', 'children'),
     Output('hayc', 'children'),
     Output('hazc', 'children'),
     Output('statc', 'children'),
     ],
    [Input(component_id='fieldx', component_property='value'),
     Input(component_id='fieldy', component_property='value'),
     Input(component_id='fieldz', component_property='value'),
     Input(component_id='L0', component_property='value'),
     Input(component_id='L1', component_property='value'),
     Input(component_id='L2', component_property='value'),
     Input(component_id='scale', component_property='value'),
     Input(component_id='k1c', component_property='value'),
     Input(component_id='k2c', component_property='value'),
     Input(component_id='ms', component_property='value'),
     Input(component_id='alpha0', component_property='value'),
     Input(component_id='tol', component_property='value'),
     Input(component_id='dt', component_property='value'),
     Input(component_id='nt', component_property='value'),

    ]
)


def update_fig(hx,hy,hz,lmb0,lmb1,lmb2,s,kk1,kk2,mms,alph,tol00,dtt,ntt):

    d=1.0

    crys = 'cub'


    ax,ay,az,stat0=llg(crys,mms,hx,hy,hz,kk1,kk2,alph,tol00,dtt,ntt)
    heffx,heffy,heffz=field(crys,ax,ay,az,mms,hx,hy,hz,kk1,kk2)
    torque=np.sqrt((heffy*az-heffz*ay)**2.0+(heffz*ax-heffx*az)**2.0+(heffx*ay-heffy*ax)**2.0)

    u, v = np.mgrid[0:np.pi:100j, 0:2.0*np.pi:100j]

    bx = np.sin(u)*np.cos(v)
    by = np.sin(u)*np.sin(v)
    bz = np.cos(u)
    f = lmb0+lmb1*1.5*(ax*ax*bx*bx+ay*ay*by*by+az*az*bz*bz-(1/3))+3*lmb2*(ax*ay*bx*by+ay*az*by*bz+ax*az*bx*bz)
    x = d*(1.0+s*f)*bx
    y = d*(1.0+s*f)*by
    z = d*(1.0+s*f)*bz
    dl_l=(np.sqrt(x**2+y**2+z**2)-d)/d
    hmod=np.sqrt(hx**2.0+hy**2.0+hz**2.0)+10.0**(-12)

    hhx=hx/hmod
    hhy=hy/hmod
    hhz=hz/hmod

    lsim = np.sqrt(x**2+y**2+z**2)
    list00 = np.stack((np.transpose(x),np.transpose(y),np.transpose(z),np.transpose(dl_l)),axis=0)
    list0 = np.stack((lsim,dl_l,dl_l/s,bx,by,bz,u*(180.0/np.pi),v*(180.0/np.pi)),axis=-1)
    hax=heffx-hx
    hay=heffy-hy
    haz=heffz-hz
    hamod=np.sqrt(hax**2.0+hay**2.0+haz**2.0)+10.0**(-12)
    hhax=hax/hamod
    hhay=hay/hamod
    hhaz=haz/hamod
    heffmod=np.sqrt(heffx**2.0+heffy**2.0+heffz**2.0)+10.0**(-12)
    hheffx=heffx/heffmod
    hheffy=heffy/heffmod
    hheffz=heffz/heffmod


    fig = make_subplots(rows=1, cols=2,
                    specs=[[{'is_3d': True}, {'is_3d': True}]],
                    subplot_titles=['  Fields (H,Ha,Heff), magnetization (\u03B1) and unit cell lattice vectors (a,b,c)           ',
                                    '        Color corresponds to (\u0394l/lo)*scale_factor along direction \u03B2(θ,φ)'],)

    fig.add_trace(go.Cone(x=[0], y=[1], z=[1], u=[hhx], v=[hhy], w=[hhz],name="""<b>H<b>/|<b>H<b>|""",
                          hovertemplate = """Hx/H = %{u:.6g}<br>Hy/H = %{v:.6g}<br>Hz/H = %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(255,0,0)']]),1, 1)

    fig.add_trace(go.Cone(x=[1], y=[0], z=[1], u=[hhax], v=[hhay], w=[hhaz],name="""<b>Ha<b>/|<b>Ha<b>|""",
                          hovertemplate = """Hax/Ha = %{u:.6g}<br>Hay/Ha = %{v:.6g}<br>Haz/Ha = %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(225,0,0)']]),1, 1)

    fig.add_trace(go.Cone(x=[1], y=[1], z=[1], u=[hheffx], v=[hheffy], w=[hheffz],name="""<b>Heff<b>/|<b>Heff<b>|""",
                          hovertemplate = """Heffx/Heff = %{u:.6g}<br>Heffy/Heff = %{v:.6g}<br>Heffz/Heff = %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(200,0,0)']]),1, 1)


    fig.add_trace(go.Cone(x=[1], y=[1], z=[0], u=[ax], v=[ay], w=[az],name="""<b>\u03B1<b>""",
                          hovertemplate = """\u03B1x = %{u:.6g}<br>\u03B1y= %{v:.6g}<br>\u03B1z= %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(255,127,14)'], [1, 'rgb(255,127,16)']]),1, 1)
    fig.add_trace(go.Cone(x=[0.4], y=[0], z=[0], u=[2], v=[0], w=[0],name="a",hoverinfo="name",
                          colorscale=[[0, 'rgb(255,0,0)'], [1, 'rgb(240,0,0)']]),1, 1)
    fig.add_trace(go.Cone(x=[0], y=[0.4], z=[0], u=[0], v=[2], w=[0],name="b",hoverinfo="name",
                          colorscale=[[0, 'rgb(0,255,0)'], [1, 'rgb(0,240,0)']]),1, 1)
    fig.add_trace(go.Cone(x=[0], y=[0], z=[0.4], u=[0], v=[0], w=[2],name="c",hoverinfo="name",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(0,0,240)']]),1, 1)

    fig.update_traces(showscale=False)

    fig.add_trace(go.Surface(x=list00[0], y=list00[1], z=list00[2],
                             surfacecolor=list00[3], customdata=list0, name=" ",hoverinfo="name",
                             hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 2)


    fig.update_layout(transition_duration=500)



    fig.update_yaxes(automargin=True)
    fig.update_xaxes(automargin=True)

    return fig,ax,ay,az,heffx,heffy,heffz,torque,hax,hay,haz,stat0

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
     Input(component_id='k1c', component_property='value'),
     Input(component_id='k2c', component_property='value'),
     Input(component_id='ms', component_property='value'),
     Input(component_id='alpha0', component_property='value'),
     Input(component_id='tol', component_property='value'),
     Input(component_id='dt', component_property='value'),
     Input(component_id='nt', component_property='value'),

    ]
)


def update_figc2d(hx,hy,hz,lmb0,lmb1,lmb2,s,kk1,kk2,mms,alph,tol00,dtt,ntt):

    d=1.0

    crys = 'cub'

    ax,ay,az,stat0=llg(crys,mms,hx,hy,hz,kk1,kk2,alph,tol00,dtt,ntt)


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



    dl_lxy = (np.sqrt(xxy**2+yxy**2+zxy**2)-d)/d
    lsimxy = np.sqrt(xxy**2+yxy**2+zxy**2)
    list00xy = np.stack((np.transpose(xxy),np.transpose(yxy),np.transpose(zxy)),axis=0)
    list0xy = np.stack((lsimxy,dl_lxy,dl_lxy/s,xxy/lsimxy,yxy/lsimxy,zxy/lsimxy, vxy*0.0,vxy*(180.0/np.pi)),axis=-1)

    dl_lxz = (np.sqrt(xxz**2+yxz**2+zxz**2)-d)/d
    lsimxz = np.sqrt(xxz**2+yxz**2+zxz**2)
    list00xz = np.stack((np.transpose(xxz),np.transpose(yxz),np.transpose(zxz)),axis=0)
    list0xz = np.stack((lsimxz,dl_lxz,dl_lxz/s,xxz/lsimxz,yxz/lsimxz,zxz/lsimxz,vxz*(180.0/np.pi),vxz*0.0),axis=-1)

    dl_lyz = (np.sqrt(xyz**2+yyz**2+zyz**2)-d)/d
    lsimyz = np.sqrt(xyz**2+yyz**2+zyz**2)
    uyz = np.mgrid[90.0:90.0:200j]
    list00yz = np.stack((np.transpose(xyz),np.transpose(yyz),np.transpose(zyz)),axis=0)
    list0yz = np.stack((lsimyz,dl_lyz,dl_lyz/s,xyz/lsimyz,yyz/lsimyz,zyz/lsimyz,vyz*(180.0/np.pi),uyz),axis=-1)



    figc2d = make_subplots(rows=1, cols=3,
                    specs=[[{'type': 'xy'}, {'type': 'xy'}, {'type': 'xy'}]],
                    subplot_titles=['Plot XY plane (z=0)', 'Plot XZ plane (y=0)', 'Plot YZ plane (x=0)'],)


    figc2d.add_trace(go.Scatter(x=bxxy, y=byxy, mode='lines',name='(lox,loy)'), 1, 1)
    figc2d.add_trace(go.Scatter(x=list00xy[0], y=list00xy[1], mode='lines', customdata=list0xy, name="(lsimx,lsimy)", hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 1)


    figc2d.add_trace(go.Scatter(x=bxxz, y=bzxz, mode='lines',name='(lox,loz)'), 1, 2)
    figc2d.add_trace(go.Scatter(x=list00xz[0], y=list00xz[2], mode='lines', customdata=list0xz, name="(lsimx,lsimz)", hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 2)


    figc2d.add_trace(go.Scatter(x=byyz, y=bzyz, mode='lines',name='(loy,loz)'), 1, 3)
    figc2d.add_trace(go.Scatter(x=list00yz[1], y=list00yz[2], mode='lines', customdata=list0yz, name="(lsimy,lsimz)", hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 3)


    figc2d.update_layout(xaxis_title='X-axis')
    figc2d.update_layout(yaxis_title='Y-axis')
    figc2d.update_layout(transition_duration=500)

    figc2d.update_yaxes(automargin=True)
    figc2d.update_xaxes(automargin=True)

    return figc2d



##############MAE Cubic I -3D

@app.callback(
    Output('mae_cub_3D', 'figure'),
    [Input(component_id='k0c', component_property='value'),
     Input(component_id='k1c', component_property='value'),
     Input(component_id='k2c', component_property='value'),
    ]
)


def update_figmaec(kk0,kk1,kk2):


    kk0=kk0*10.0**3
    kk1=kk1*10.0**3
    kk2=kk2*10.0**3


    u, v = np.mgrid[0:np.pi:100j, 0:2*np.pi:100j]
    bx = np.sin(u)*np.cos(v)
    by = np.sin(u)*np.sin(v)
    bz = np.cos(u)
    f = kk0+kk1*(bx**2.0*by**2.0+by**2.0*bz**2.0+bx**2.0*bz**2.0)+kk2*bx**2.0*by**2.0*bz**2.0
    f0 = f-kk0
    f0min = np.amin(f0)
    f0max = np.amax(f0)
    ff0 = f0 + np.absolute(f0min) + (f0max-f0min)*0.45

    x = ff0*bx
    y = ff0*by
    z = ff0*bz



    list00 = np.stack((np.transpose(x),np.transpose(y),np.transpose(z),np.transpose(f)),axis=0)
    list0 = np.stack((f,bx,by,bz,u*(180.0/np.pi),v*(180.0/np.pi)),axis=-1)



    figmaec = make_subplots(rows=1, cols=1,
                    specs=[[{'is_3d': True}]],
                    subplot_titles=['  Color corresponds to Magnetocrystalline Anisotropy Energy (KJ/m^3) '])


    figmaec.add_trace(go.Surface(x=list00[0], y=list00[1], z=list00[2], surfacecolor=list00[3], customdata=list0, name=" ", hoverinfo="name",
                                 hovertemplate = """Ea = %{customdata[0]:.6g} KJ/m^3<br>\u03B1x = %{customdata[1]:.6g}<br>\u03B1y = %{customdata[2]:.6g}<br>\u03B1z = %{customdata[3]:.6g}<br>θ = %{customdata[4]:.6g}°<br>φ = %{customdata[5]:.6g}°<br>"""),1, 1)


    figmaec.update_layout(transition_duration=500)

    figmaec.update_layout(
        scene = {
            "xaxis": {"showticklabels": False},
            "yaxis": {"showticklabels": False},
            "zaxis": {"showticklabels": False}
        })



    return figmaec




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
    h=np.sqrt(hx*hx+hy*hy+hz*hz)+1.0e-12
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

    dl_l=(np.sqrt(x**2+y**2+z**2)-d)/d
    lsim = np.sqrt(x**2+y**2+z**2)
    list00 = np.stack((np.transpose(x),np.transpose(y),np.transpose(z),np.transpose(dl_l)),axis=0)
    list0 = np.stack((lsim,dl_l,dl_l/s,bx,by,bz,u*(180.0/np.pi),v*(180.0/np.pi)),axis=-1)

    figp = make_subplots(rows=1, cols=2,
                    specs=[[{'is_3d': True}, {'is_3d': True}]],
                    subplot_titles=['       Effective field (Heff), magnetization (\u03B1) and axis (x,y,z)           ', '        Color corresponds to (\u0394l/lo)*scale_factor along direction \u03B2'],)

    figp.add_trace(go.Cone(x=[1], y=[1], z=[1], u=[ax], v=[ay], w=[az],name="""<b>Heff<b>/|<b>Heff<b>|""",
                          hovertemplate = """Heffx/Heff = %{u:.6g}<br>Heffy/Heff = %{v:.6g}<br>Heffz/Heff = %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(200,0,0)']]),1, 1)

    figp.add_trace(go.Cone(x=[1], y=[1], z=[0], u=[ax], v=[ay], w=[az],name="""<b>\u03B1<b>""",
                          hovertemplate = """\u03B1x = %{u:.6g}<br>\u03B1y= %{v:.6g}<br>\u03B1z= %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(255,127,14)'], [1, 'rgb(255,127,16)']]),1, 1)
    figp.add_trace(go.Cone(x=[0.4], y=[0], z=[0], u=[2], v=[0], w=[0],name="X",hoverinfo="name",
                          colorscale=[[0, 'rgb(255,0,0)'], [1, 'rgb(240,0,0)']]),1, 1)
    figp.add_trace(go.Cone(x=[0], y=[0.4], z=[0], u=[0], v=[2], w=[0],name="Y",hoverinfo="name",
                          colorscale=[[0, 'rgb(0,255,0)'], [1, 'rgb(0,240,0)']]),1, 1)
    figp.add_trace(go.Cone(x=[0], y=[0], z=[0.4], u=[0], v=[0], w=[2],name="Z",hoverinfo="name",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(0,0,240)']]),1, 1)

    figp.update_traces(showscale=False)

    figp.add_trace(go.Surface(x=list00[0], y=list00[1], z=list00[2],
                             surfacecolor=list00[3], customdata=list0, name=" ",hoverinfo="name",
                             hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 2)

    figp.update_layout(transition_duration=500)


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


    dl_lxy = (np.sqrt(xxy**2+yxy**2+zxy**2)-d)/d
    lsimxy = np.sqrt(xxy**2+yxy**2+zxy**2)
    list00xy = np.stack((np.transpose(xxy),np.transpose(yxy),np.transpose(zxy)),axis=0)
    list0xy = np.stack((lsimxy,dl_lxy,dl_lxy/s,xxy/lsimxy,yxy/lsimxy,zxy/lsimxy, vxy*0.0,vxy*(180.0/np.pi)),axis=-1)

    dl_lxz = (np.sqrt(xxz**2+yxz**2+zxz**2)-d)/d
    lsimxz = np.sqrt(xxz**2+yxz**2+zxz**2)
    list00xz = np.stack((np.transpose(xxz),np.transpose(yxz),np.transpose(zxz)),axis=0)
    list0xz = np.stack((lsimxz,dl_lxz,dl_lxz/s,xxz/lsimxz,yxz/lsimxz,zxz/lsimxz,vxz*(180.0/np.pi),vxz*0.0),axis=-1)

    dl_lyz = (np.sqrt(xyz**2+yyz**2+zyz**2)-d)/d
    lsimyz = np.sqrt(xyz**2+yyz**2+zyz**2)
    uyz = np.mgrid[90.0:90.0:200j]
    list00yz = np.stack((np.transpose(xyz),np.transpose(yyz),np.transpose(zyz)),axis=0)
    list0yz = np.stack((lsimyz,dl_lyz,dl_lyz/s,xyz/lsimyz,yyz/lsimyz,zyz/lsimyz,vyz*(180.0/np.pi),uyz),axis=-1)



    figcp2d = make_subplots(rows=1, cols=3,
                    specs=[[{'type': 'xy'}, {'type': 'xy'}, {'type': 'xy'}]],
                    subplot_titles=['Plot XY plane (z=0)', 'Plot XZ plane (y=0)', 'Plot YZ plane (x=0)'],)


    figcp2d.add_trace(go.Scatter(x=bxxy, y=byxy, mode='lines',name='(lox,loy)'), 1, 1)
    figcp2d.add_trace(go.Scatter(x=list00xy[0], y=list00xy[1], mode='lines', customdata=list0xy, name="(lsimx,lsimy)", hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 1)


    figcp2d.add_trace(go.Scatter(x=bxxz, y=bzxz, mode='lines',name='(lox,loz)'), 1, 2)
    figcp2d.add_trace(go.Scatter(x=list00xz[0], y=list00xz[2], mode='lines', customdata=list0xz, name="(lsimx,lsimz)", hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 2)


    figcp2d.add_trace(go.Scatter(x=byyz, y=bzyz, mode='lines',name='(loy,loz)'), 1, 3)
    figcp2d.add_trace(go.Scatter(x=list00yz[1], y=list00yz[2], mode='lines', customdata=list0yz, name="(lsimy,lsimz)", hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 3)


    figcp2d.update_layout(xaxis_title='X-axis')
    figcp2d.update_layout(yaxis_title='Y-axis')
    figcp2d.update_layout(transition_duration=500)

    figcp2d.update_yaxes(automargin=True)
    figcp2d.update_xaxes(automargin=True)



    return figcp2d



############## hex-3D-clark


@app.callback(
    [Output('hex_3D', 'figure'),
     Output('meqxh', 'children'),
     Output('meqyh', 'children'),
     Output('meqzh', 'children'),
     Output('hefxh', 'children'),
     Output('hefyh', 'children'),
     Output('hefzh', 'children'),
     Output('torqueh', 'children'),
     Output('haxh', 'children'),
     Output('hayh', 'children'),
     Output('hazh', 'children'),
     Output('stath', 'children'),
     ],
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
     Input(component_id='k1h', component_property='value'),
     Input(component_id='k2h', component_property='value'),
     Input(component_id='msh', component_property='value'),
     Input(component_id='alpha0h', component_property='value'),
     Input(component_id='tolh', component_property='value'),
     Input(component_id='dth', component_property='value'),
     Input(component_id='nth', component_property='value'),

    ]
)



def update_hfig(hx,hy,hz,lmb01,lmb02,lmba1,lmba2,lmbg,lmbe,s,kk1,kk2,mms,alph,tol00,dtt,ntt):

    d=1.0

    crys = 'uni'

    ax,ay,az,stat0=llg(crys,mms,hx,hy,hz,kk1,kk2,alph,tol00,dtt,ntt)
    heffx,heffy,heffz=field(crys,ax,ay,az,mms,hx,hy,hz,kk1,kk2)
    torque=np.sqrt((heffy*az-heffz*ay)**2.0+(heffz*ax-heffx*az)**2.0+(heffx*ay-heffy*ax)**2.0)



    u, v = np.mgrid[0:np.pi:100j, 0:2*np.pi:100j]
    bx = np.sin(u)*np.cos(v)
    by = np.sin(u)*np.sin(v)
    bz = np.cos(u)
    f = lmb01*(bx**2+by**2)+lmb02*(bz**2)+lmba1*((az**2)-(1/3))*(bx**2+by**2)+lmba2*((az**2)-(1/3))*(bz**2)+lmbg*(0.5*(ax**2-ay**2)*(bx**2-by**2)+2*ax*ay*bx*by)+2*lmbe*(ax*az*bx*bz+az*ay*bz*by)
    x = d*(1.0+s*f)*bx
    y = d*(1.0+s*f)*by
    z = d*(1.0+s*f)*bz
    dl_l=(np.sqrt(x**2+y**2+z**2)-d)/d

    hmod=np.sqrt(hx**2.0+hy**2.0+hz**2.0)+10.0**(-8)

    hhx=hx/hmod
    hhy=hy/hmod
    hhz=hz/hmod


    lsim = np.sqrt(x**2+y**2+z**2)
    list00 = np.stack((np.transpose(x),np.transpose(y),np.transpose(z),np.transpose(dl_l)),axis=0)
    list0 = np.stack((lsim,dl_l,dl_l/s,bx,by,bz,u*(180.0/np.pi),v*(180.0/np.pi)),axis=-1)
    hax=heffx-hx
    hay=heffy-hy
    haz=heffz-hz
    hamod=np.sqrt(hax**2.0+hay**2.0+haz**2.0)+10.0**(-12)
    hhax=hax/hamod
    hhay=hay/hamod
    hhaz=haz/hamod
    heffmod=np.sqrt(heffx**2.0+heffy**2.0+heffz**2.0)+10.0**(-12)
    hheffx=heffx/heffmod
    hheffy=heffy/heffmod
    hheffz=heffz/heffmod


    hfig = make_subplots(rows=1, cols=2,
                    specs=[[{'is_3d': True}, {'is_3d': True}]],
                    subplot_titles=['  Fields (H,Ha,Heff), magnetization (\u03B1) and unit cell lattice vectors (a,b,c)           ',
                                    '        Color corresponds to (\u0394l/lo)*scale_factor along direction \u03B2(θ,φ)'],)

    hfig.add_trace(go.Cone(x=[0], y=[1], z=[1], u=[hhx], v=[hhy], w=[hhz],name="""<b>H<b>/|<b>H<b>|""",
                          hovertemplate = """Hx/H = %{u:.6g}<br>Hy/H = %{v:.6g}<br>Hz/H = %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(255,0,0)']]),1, 1)

    hfig.add_trace(go.Cone(x=[1], y=[0], z=[1], u=[hhax], v=[hhay], w=[hhaz],name="""<b>Ha<b>/|<b>Ha<b>|""",
                          hovertemplate = """Hax/Ha = %{u:.6g}<br>Hay/Ha = %{v:.6g}<br>Haz/Ha = %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(225,0,0)']]),1, 1)

    hfig.add_trace(go.Cone(x=[1], y=[1], z=[1], u=[hheffx], v=[hheffy], w=[hheffz],name="""<b>Heff<b>/|<b>Heff<b>|""",
                          hovertemplate = """Heffx/Heff = %{u:.6g}<br>Heffy/Heff = %{v:.6g}<br>Heffz/Heff = %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(200,0,0)']]),1, 1)


    hfig.add_trace(go.Cone(x=[1], y=[1], z=[0], u=[ax], v=[ay], w=[az],name="""<b>\u03B1<b>""",
                          hovertemplate = """\u03B1x = %{u:.6g}<br>\u03B1y= %{v:.6g}<br>\u03B1z= %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(255,127,14)'], [1, 'rgb(255,127,16)']]),1, 1)
    hfig.add_trace(go.Cone(x=[0.4], y=[0], z=[0], u=[2], v=[0], w=[0],name="a",hoverinfo="name",
                          colorscale=[[0, 'rgb(255,0,0)'], [1, 'rgb(240,0,0)']]),1, 1)
    hfig.add_trace(go.Cone(x=[-0.2], y=[0.34641016], z=[0], u=[-1], v=[1.7320508],w=[0],name="b",hoverinfo="name",
                          colorscale=[[0, 'rgb(0,255,0)'], [1, 'rgb(0,240,0)']]),1, 1)
    hfig.add_trace(go.Cone(x=[0], y=[0], z=[0.4], u=[0], v=[0], w=[2],name="c",hoverinfo="name",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(0,0,240)']]),1, 1)

    hfig.update_traces(showscale=False)

    hfig.add_trace(go.Surface(x=list00[0], y=list00[1], z=list00[2],
                             surfacecolor=list00[3], customdata=list0, name=" ",hoverinfo="name",
                             hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 2)


    hfig.update_layout(transition_duration=500)



    hfig.update_yaxes(automargin=True)
    hfig.update_xaxes(automargin=True)

    return hfig,ax,ay,az,heffx,heffy,heffz,torque,hax,hay,haz,stat0


######################################## Hex 2D - clark


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
     Input(component_id='k1h', component_property='value'),
     Input(component_id='k2h', component_property='value'),
     Input(component_id='msh', component_property='value'),
     Input(component_id='alpha0h', component_property='value'),
     Input(component_id='tolh', component_property='value'),
     Input(component_id='dth', component_property='value'),
     Input(component_id='nth', component_property='value'),
    ]
)


def update_figh2d(hx,hy,hz,lmb01,lmb02,lmba1,lmba2,lmbg,lmbe,s,kkk1,kkk2,mmms,alph,tol00,dtt,ntt):

    d=1.0

    crys = 'uni'

    ax,ay,az,stat0=llg(crys,mmms,hx,hy,hz,kkk1,kkk2,alph,tol00,dtt,ntt)

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

    dl_lxy = (np.sqrt(xxy**2+yxy**2+zxy**2)-d)/d
    lsimxy = np.sqrt(xxy**2+yxy**2+zxy**2)
    list00xy = np.stack((np.transpose(xxy),np.transpose(yxy),np.transpose(zxy)),axis=0)
    list0xy = np.stack((lsimxy,dl_lxy,dl_lxy/s,xxy/lsimxy,yxy/lsimxy,zxy/lsimxy, vxy*0.0,vxy*(180.0/np.pi)),axis=-1)

    dl_lxz = (np.sqrt(xxz**2+yxz**2+zxz**2)-d)/d
    lsimxz = np.sqrt(xxz**2+yxz**2+zxz**2)
    list00xz = np.stack((np.transpose(xxz),np.transpose(yxz),np.transpose(zxz)),axis=0)
    list0xz = np.stack((lsimxz,dl_lxz,dl_lxz/s,xxz/lsimxz,yxz/lsimxz,zxz/lsimxz,vxz*(180.0/np.pi),vxz*0.0),axis=-1)

    dl_lyz = (np.sqrt(xyz**2+yyz**2+zyz**2)-d)/d
    lsimyz = np.sqrt(xyz**2+yyz**2+zyz**2)
    uyz = np.mgrid[90.0:90.0:200j]
    list00yz = np.stack((np.transpose(xyz),np.transpose(yyz),np.transpose(zyz)),axis=0)
    list0yz = np.stack((lsimyz,dl_lyz,dl_lyz/s,xyz/lsimyz,yyz/lsimyz,zyz/lsimyz,vyz*(180.0/np.pi),uyz),axis=-1)

    figh2d = make_subplots(rows=1, cols=3,
                    specs=[[{'type': 'xy'}, {'type': 'xy'}, {'type': 'xy'}]],
                    subplot_titles=['Plot XY plane (z=0)', 'Plot XZ plane (y=0)', 'Plot YZ plane (x=0)'],)

    figh2d.add_trace(go.Scatter(x=bxxy, y=byxy, mode='lines',name='(lox,loy)'), 1, 1)
    figh2d.add_trace(go.Scatter(x=list00xy[0], y=list00xy[1], mode='lines', customdata=list0xy, name="(lsimx,lsimy)", hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 1)

    figh2d.add_trace(go.Scatter(x=bxxz, y=bzxz, mode='lines',name='(lox,loz)'), 1, 2)
    figh2d.add_trace(go.Scatter(x=list00xz[0], y=list00xz[2], mode='lines', customdata=list0xz, name="(lsimx,lsimz)", hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 2)

    figh2d.add_trace(go.Scatter(x=byyz, y=bzyz, mode='lines',name='(loy,loz)'), 1, 3)
    figh2d.add_trace(go.Scatter(x=list00yz[1], y=list00yz[2], mode='lines', customdata=list0yz, name="(lsimy,lsimz)", hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 3)

    figh2d.update_layout(xaxis_title='X-axis')
    figh2d.update_layout(yaxis_title='Y-axis')
    figh2d.update_layout(transition_duration=500)

    figh2d.update_yaxes(automargin=True)
    figh2d.update_xaxes(automargin=True)

    return figh2d


##############MAE hex -3D - clark

@app.callback(
    Output('mae_hex_3D', 'figure'),
    [Input(component_id='k0h', component_property='value'),
     Input(component_id='k1h', component_property='value'),
     Input(component_id='k2h', component_property='value'),
    ]
)


def update_figmaeh(kkkk0,kkkk1,kkkk2):


    kkkk0=kkkk0*10**3
    kkkk1=kkkk1*10**3
    kkkk2=kkkk2*10**3


    u, v = np.mgrid[0:np.pi:100j, 0:2*np.pi:100j]
    bx = np.sin(u)*np.cos(v)
    by = np.sin(u)*np.sin(v)
    bz = np.cos(u)
    f = kkkk0+kkkk1*(1.0-bz**2.0)+kkkk2*(1.0-bz**2.0)**2.0
    f0 = f-kkkk0
    f0min = np.amin(f0)
    f0max = np.amax(f0)
    ff0 = f0 + np.absolute(f0min) + (f0max-f0min)*0.45

    x = ff0*bx
    y = ff0*by
    z = ff0*bz


    list00 = np.stack((np.transpose(x),np.transpose(y),np.transpose(z),np.transpose(f)),axis=0)
    list0 = np.stack((f,bx,by,bz,u*(180.0/np.pi),v*(180.0/np.pi)),axis=-1)



    figmaeh = make_subplots(rows=1, cols=1,
                    specs=[[{'is_3d': True}]],
                    subplot_titles=['  Color corresponds to Magnetocrystalline Anisotropy Energy (KJ/m^3) '])


    figmaeh.add_trace(go.Surface(x=list00[0], y=list00[1], z=list00[2], surfacecolor=list00[3], customdata=list0, name=" ", hoverinfo="name",
                                 hovertemplate = """Ea = %{customdata[0]:.6g} KJ/m^3<br>\u03B1x = %{customdata[1]:.6g}<br>\u03B1y = %{customdata[2]:.6g}<br>\u03B1z = %{customdata[3]:.6g}<br>θ = %{customdata[4]:.6g}°<br>φ = %{customdata[5]:.6g}°<br>"""), 1, 1)


    figmaeh.update_layout(transition_duration=500)

    figmaeh.update_layout(
        scene = {
            "xaxis": {"showticklabels": False},
            "yaxis": {"showticklabels": False},
            "zaxis": {"showticklabels": False}
        })

    return figmaeh



############## hex-3D-mason


@app.callback(
    [Output('hex_3Dm', 'figure'),
     Output('meqxhm', 'children'),
     Output('meqyhm', 'children'),
     Output('meqzhm', 'children'),
     Output('hefxhm', 'children'),
     Output('hefyhm', 'children'),
     Output('hefzhm', 'children'),
     Output('torquehm', 'children'),
     Output('haxhm', 'children'),
     Output('hayhm', 'children'),
     Output('hazhm', 'children'),
     Output('stathm', 'children'),
     ],
    [Input(component_id='hfieldxm', component_property='value'),
     Input(component_id='hfieldym', component_property='value'),
     Input(component_id='hfieldzm', component_property='value'),
     Input(component_id='hL01m', component_property='value'),
     Input(component_id='hL02m', component_property='value'),
     Input(component_id='hLa', component_property='value'),
     Input(component_id='hLb', component_property='value'),
     Input(component_id='hLc', component_property='value'),
     Input(component_id='hLd', component_property='value'),
     Input(component_id='hscalem', component_property='value'),
     Input(component_id='k1hm', component_property='value'),
     Input(component_id='k2hm', component_property='value'),
     Input(component_id='mshm', component_property='value'),
     Input(component_id='alpha0hm', component_property='value'),
     Input(component_id='tolhm', component_property='value'),
     Input(component_id='dthm', component_property='value'),
     Input(component_id='nthm', component_property='value'),

    ]
)



def update_hfigm(hx,hy,hz,lmb01,lmb02,lmba,lmbb,lmbc,lmbd,s,kk1,kk2,mms,alph,tol00,dtt,ntt):

    d=1.0

    crys = 'uni'

    ax,ay,az,stat0=llg(crys,mms,hx,hy,hz,kk1,kk2,alph,tol00,dtt,ntt)
    heffx,heffy,heffz=field(crys,ax,ay,az,mms,hx,hy,hz,kk1,kk2)
    torque=np.sqrt((heffy*az-heffz*ay)**2.0+(heffz*ax-heffx*az)**2.0+(heffx*ay-heffy*ax)**2.0)



    u, v = np.mgrid[0:np.pi:100j, 0:2*np.pi:100j]
    bx = np.sin(u)*np.cos(v)
    by = np.sin(u)*np.sin(v)
    bz = np.cos(u)

    f = lmb01*(bx**2+by**2)+lmb02*(bz**2)+lmba*((ax*bx+ay*by)**2-(ax*bx+ay*by)*az*bz)+lmbb*((1.0-az**2)*(1.0-bz**2)-(ax*bx+ay*by)**2)+lmbc*((1.0-az**2)*bz**2-(ax*bx+ay*by)*az*bz)+4.0*lmbd*(ax*bx+ay*by)*az*bz

    x = d*(1.0+s*f)*bx
    y = d*(1.0+s*f)*by
    z = d*(1.0+s*f)*bz
    dl_l=(np.sqrt(x**2+y**2+z**2)-d)/d

    hmod=np.sqrt(hx**2.0+hy**2.0+hz**2.0)+10.0**(-8)

    hhx=hx/hmod
    hhy=hy/hmod
    hhz=hz/hmod


    lsim = np.sqrt(x**2+y**2+z**2)
    list00 = np.stack((np.transpose(x),np.transpose(y),np.transpose(z),np.transpose(dl_l)),axis=0)
    list0 = np.stack((lsim,dl_l,dl_l/s,bx,by,bz,u*(180.0/np.pi),v*(180.0/np.pi)),axis=-1)
    hax=heffx-hx
    hay=heffy-hy
    haz=heffz-hz
    hamod=np.sqrt(hax**2.0+hay**2.0+haz**2.0)+10.0**(-12)
    hhax=hax/hamod
    hhay=hay/hamod
    hhaz=haz/hamod
    heffmod=np.sqrt(heffx**2.0+heffy**2.0+heffz**2.0)+10.0**(-12)
    hheffx=heffx/heffmod
    hheffy=heffy/heffmod
    hheffz=heffz/heffmod


    hfigm = make_subplots(rows=1, cols=2,
                    specs=[[{'is_3d': True}, {'is_3d': True}]],
                    subplot_titles=['  Fields (H,Ha,Heff), magnetization (\u03B1) and unit cell lattice vectors (a,b,c)           ',
                                    '        Color corresponds to (\u0394l/lo)*scale_factor along direction \u03B2(θ,φ)'],)

    hfigm.add_trace(go.Cone(x=[0], y=[1], z=[1], u=[hhx], v=[hhy], w=[hhz],name="""<b>H<b>/|<b>H<b>|""",
                          hovertemplate = """Hx/H = %{u:.6g}<br>Hy/H = %{v:.6g}<br>Hz/H = %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(255,0,0)']]),1, 1)

    hfigm.add_trace(go.Cone(x=[1], y=[0], z=[1], u=[hhax], v=[hhay], w=[hhaz],name="""<b>Ha<b>/|<b>Ha<b>|""",
                          hovertemplate = """Hax/Ha = %{u:.6g}<br>Hay/Ha = %{v:.6g}<br>Haz/Ha = %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(225,0,0)']]),1, 1)

    hfigm.add_trace(go.Cone(x=[1], y=[1], z=[1], u=[hheffx], v=[hheffy], w=[hheffz],name="""<b>Heff<b>/|<b>Heff<b>|""",
                          hovertemplate = """Heffx/Heff = %{u:.6g}<br>Heffy/Heff = %{v:.6g}<br>Heffz/Heff = %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(200,0,0)']]),1, 1)


    hfigm.add_trace(go.Cone(x=[1], y=[1], z=[0], u=[ax], v=[ay], w=[az],name="""<b>\u03B1<b>""",
                          hovertemplate = """\u03B1x = %{u:.6g}<br>\u03B1y= %{v:.6g}<br>\u03B1z= %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(255,127,14)'], [1, 'rgb(255,127,16)']]),1, 1)
    hfigm.add_trace(go.Cone(x=[0.4], y=[0], z=[0], u=[2], v=[0], w=[0],name="a",hoverinfo="name",
                          colorscale=[[0, 'rgb(255,0,0)'], [1, 'rgb(240,0,0)']]),1, 1)
    hfigm.add_trace(go.Cone(x=[-0.2], y=[0.34641016], z=[0], u=[-1], v=[1.7320508],w=[0],name="b",hoverinfo="name",
                          colorscale=[[0, 'rgb(0,255,0)'], [1, 'rgb(0,240,0)']]),1, 1)
    hfigm.add_trace(go.Cone(x=[0], y=[0], z=[0.4], u=[0], v=[0], w=[2],name="c",hoverinfo="name",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(0,0,240)']]),1, 1)

    hfigm.update_traces(showscale=False)

    hfigm.add_trace(go.Surface(x=list00[0], y=list00[1], z=list00[2],
                             surfacecolor=list00[3], customdata=list0, name=" ",hoverinfo="name",
                             hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 2)


    hfigm.update_layout(transition_duration=500)



    hfigm.update_yaxes(automargin=True)
    hfigm.update_xaxes(automargin=True)

    return hfigm,ax,ay,az,heffx,heffy,heffz,torque,hax,hay,haz,stat0



######################################## Hex 2D - mason


@app.callback(
    Output('hex_2Dm', 'figure'),
    [Input(component_id='hfieldxm', component_property='value'),
     Input(component_id='hfieldym', component_property='value'),
     Input(component_id='hfieldzm', component_property='value'),
     Input(component_id='hL01m', component_property='value'),
     Input(component_id='hL02m', component_property='value'),
     Input(component_id='hLa', component_property='value'),
     Input(component_id='hLb', component_property='value'),
     Input(component_id='hLc', component_property='value'),
     Input(component_id='hLd', component_property='value'),
     Input(component_id='hscalem', component_property='value'),
     Input(component_id='k1hm', component_property='value'),
     Input(component_id='k2hm', component_property='value'),
     Input(component_id='mshm', component_property='value'),
     Input(component_id='alpha0hm', component_property='value'),
     Input(component_id='tolhm', component_property='value'),
     Input(component_id='dthm', component_property='value'),
     Input(component_id='nthm', component_property='value'),
    ]
)


def update_figh2dm(hx,hy,hz,lmb01,lmb02,lmba,lmbb,lmbc,lmbd,s,kkk1,kkk2,mmms,alph,tol00,dtt,ntt):

    d=1.0

    crys = 'uni'

    ax,ay,az,stat0=llg(crys,mmms,hx,hy,hz,kkk1,kkk2,alph,tol00,dtt,ntt)

    vxy = np.mgrid[0:2*np.pi:200j]
    bxxy = np.cos(vxy)
    byxy = np.sin(vxy)
    bzxy = 0.0
    fxy = lmb01*(bxxy**2+byxy**2)+lmb02*(bzxy**2)+lmba*((ax*bxxy+ay*byxy)**2-(ax*bxxy+ay*byxy)*az*bzxy)+lmbb*((1.0-az**2)*(1.0-bzxy**2)-(ax*bxxy+ay*byxy)**2)+lmbc*((1.0-az**2)*bzxy**2-(ax*bxxy+ay*byxy)*az*bzxy)+4.0*lmbd*(ax*bxxy+ay*byxy)*az*bzxy
    xxy = d*(1.0+s*fxy)*bxxy
    yxy = d*(1.0+s*fxy)*byxy
    zxy = d*(1.0+s*fxy)*bzxy

    vxz = np.mgrid[0:2*np.pi:200j]
    bxxz = np.sin(vxz)
    byxz = 0.0
    bzxz = np.cos(vxz)
    fxz = lmb01*(bxxz**2+byxz**2)+lmb02*(bzxz**2)+lmba*((ax*bxxz+ay*byxz)**2-(ax*bxxz+ay*byxz)*az*bzxz)+lmbb*((1.0-az**2)*(1.0-bzxz**2)-(ax*bxxz+ay*byxz)**2)+lmbc*((1.0-az**2)*bzxz**2-(ax*bxxz+ay*byxz)*az*bzxz)+4.0*lmbd*(ax*bxxz+ay*byxz)*az*bzxz
    xxz = d*(1.0+s*fxz)*bxxz
    yxz = d*(1.0+s*fxz)*byxz
    zxz = d*(1.0+s*fxz)*bzxz

    vyz = np.mgrid[0:2*np.pi:200j]
    bxyz = 0.0
    byyz = np.sin(vyz)
    bzyz = np.cos(vyz)
    fyz = lmb01*(bxyz**2+byyz**2)+lmb02*(bzyz**2)+lmba*((ax*bxyz+ay*byyz)**2-(ax*bxyz+ay*byyz)*az*bzyz)+lmbb*((1.0-az**2)*(1.0-bzyz**2)-(ax*bxyz+ay*byyz)**2)+lmbc*((1.0-az**2)*bzyz**2-(ax*bxyz+ay*byyz)*az*bzyz)+4.0*lmbd*(ax*bxyz+ay*byyz)*az*bzyz
    xyz = d*(1.0+s*fyz)*bxyz
    yyz = d*(1.0+s*fyz)*byyz
    zyz = d*(1.0+s*fyz)*bzyz

    dl_lxy = (np.sqrt(xxy**2+yxy**2+zxy**2)-d)/d
    lsimxy = np.sqrt(xxy**2+yxy**2+zxy**2)
    list00xy = np.stack((np.transpose(xxy),np.transpose(yxy),np.transpose(zxy)),axis=0)
    list0xy = np.stack((lsimxy,dl_lxy,dl_lxy/s,xxy/lsimxy,yxy/lsimxy,zxy/lsimxy, vxy*0.0,vxy*(180.0/np.pi)),axis=-1)

    dl_lxz = (np.sqrt(xxz**2+yxz**2+zxz**2)-d)/d
    lsimxz = np.sqrt(xxz**2+yxz**2+zxz**2)
    list00xz = np.stack((np.transpose(xxz),np.transpose(yxz),np.transpose(zxz)),axis=0)
    list0xz = np.stack((lsimxz,dl_lxz,dl_lxz/s,xxz/lsimxz,yxz/lsimxz,zxz/lsimxz,vxz*(180.0/np.pi),vxz*0.0),axis=-1)

    dl_lyz = (np.sqrt(xyz**2+yyz**2+zyz**2)-d)/d
    lsimyz = np.sqrt(xyz**2+yyz**2+zyz**2)
    uyz = np.mgrid[90.0:90.0:200j]
    list00yz = np.stack((np.transpose(xyz),np.transpose(yyz),np.transpose(zyz)),axis=0)
    list0yz = np.stack((lsimyz,dl_lyz,dl_lyz/s,xyz/lsimyz,yyz/lsimyz,zyz/lsimyz,vyz*(180.0/np.pi),uyz),axis=-1)

    figh2dm = make_subplots(rows=1, cols=3,
                    specs=[[{'type': 'xy'}, {'type': 'xy'}, {'type': 'xy'}]],
                    subplot_titles=['Plot XY plane (z=0)', 'Plot XZ plane (y=0)', 'Plot YZ plane (x=0)'],)

    figh2dm.add_trace(go.Scatter(x=bxxy, y=byxy, mode='lines',name='(lox,loy)'), 1, 1)
    figh2dm.add_trace(go.Scatter(x=list00xy[0], y=list00xy[1], mode='lines', customdata=list0xy, name="(lsimx,lsimy)", hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 1)

    figh2dm.add_trace(go.Scatter(x=bxxz, y=bzxz, mode='lines',name='(lox,loz)'), 1, 2)
    figh2dm.add_trace(go.Scatter(x=list00xz[0], y=list00xz[2], mode='lines', customdata=list0xz, name="(lsimx,lsimz)", hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 2)

    figh2dm.add_trace(go.Scatter(x=byyz, y=bzyz, mode='lines',name='(loy,loz)'), 1, 3)
    figh2dm.add_trace(go.Scatter(x=list00yz[1], y=list00yz[2], mode='lines', customdata=list0yz, name="(lsimy,lsimz)", hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 3)

    figh2dm.update_layout(xaxis_title='X-axis')
    figh2dm.update_layout(yaxis_title='Y-axis')
    figh2dm.update_layout(transition_duration=500)

    figh2dm.update_yaxes(automargin=True)
    figh2dm.update_xaxes(automargin=True)

    return figh2dm



##############MAE hex -3D - mason

@app.callback(
    Output('mae_hex_3Dm', 'figure'),
    [Input(component_id='k0hm', component_property='value'),
     Input(component_id='k1hm', component_property='value'),
     Input(component_id='k2hm', component_property='value'),
    ]
)


def update_figmaehm(kkkk0,kkkk1,kkkk2):


    kkkk0=kkkk0*10**3
    kkkk1=kkkk1*10**3
    kkkk2=kkkk2*10**3


    u, v = np.mgrid[0:np.pi:100j, 0:2*np.pi:100j]
    bx = np.sin(u)*np.cos(v)
    by = np.sin(u)*np.sin(v)
    bz = np.cos(u)
    f = kkkk0+kkkk1*(1.0-bz**2.0)+kkkk2*(1.0-bz**2.0)**2.0
    f0 = f-kkkk0
    f0min = np.amin(f0)
    f0max = np.amax(f0)
    ff0 = f0 + np.absolute(f0min) + (f0max-f0min)*0.45

    x = ff0*bx
    y = ff0*by
    z = ff0*bz


    list00 = np.stack((np.transpose(x),np.transpose(y),np.transpose(z),np.transpose(f)),axis=0)
    list0 = np.stack((f,bx,by,bz,u*(180.0/np.pi),v*(180.0/np.pi)),axis=-1)



    figmaehm = make_subplots(rows=1, cols=1,
                    specs=[[{'is_3d': True}]],
                    subplot_titles=['  Color corresponds to Magnetocrystalline Anisotropy Energy (KJ/m^3) '])


    figmaehm.add_trace(go.Surface(x=list00[0], y=list00[1], z=list00[2], surfacecolor=list00[3], customdata=list0, name=" ", hoverinfo="name",
                                 hovertemplate = """Ea = %{customdata[0]:.6g} KJ/m^3<br>\u03B1x = %{customdata[1]:.6g}<br>\u03B1y = %{customdata[2]:.6g}<br>\u03B1z = %{customdata[3]:.6g}<br>θ = %{customdata[4]:.6g}°<br>φ = %{customdata[5]:.6g}°<br>"""), 1, 1)


    figmaehm.update_layout(transition_duration=500)

    figmaehm.update_layout(
        scene = {
            "xaxis": {"showticklabels": False},
            "yaxis": {"showticklabels": False},
            "zaxis": {"showticklabels": False}
        })

    return figmaehm



############## tri-3D



@app.callback(
    [Output('tri_3D', 'figure'),
     Output('meqxtr', 'children'),
     Output('meqytr', 'children'),
     Output('meqztr', 'children'),
     Output('hefxtr', 'children'),
     Output('hefytr', 'children'),
     Output('hefztr', 'children'),
     Output('haxtr', 'children'),
     Output('haytr', 'children'),
     Output('haztr', 'children'),
     Output('torquetr', 'children'),
     Output('stattr', 'children'),
     ],
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
     Input(component_id='k1tr', component_property='value'),
     Input(component_id='k2tr', component_property='value'),
     Input(component_id='mstr', component_property='value'),
     Input(component_id='alpha0tr', component_property='value'),
     Input(component_id='toltr', component_property='value'),
     Input(component_id='dttr', component_property='value'),
     Input(component_id='nttr', component_property='value'),

    ]
)



def update_trfig(hx,hy,hz,lmb01,lmb02,lmba1,lmba2,lmbg1,lmbg2,lmb12,lmb21,s,kk1,kk2,mms,alph,tol00,dtt,ntt):

    d=1.0

    crys = 'uni'

    ax,ay,az,stat0=llg(crys,mms,hx,hy,hz,kk1,kk2,alph,tol00,dtt,ntt)
    heffx,heffy,heffz=field(crys,ax,ay,az,mms,hx,hy,hz,kk1,kk2)
    torque=np.sqrt((heffy*az-heffz*ay)**2.0+(heffz*ax-heffx*az)**2.0+(heffx*ay-heffy*ax)**2.0)


    u, v = np.mgrid[0:np.pi:100j, 0:2*np.pi:100j]
    bx = np.sin(u)*np.cos(v)
    by = np.sin(u)*np.sin(v)
    bz = np.cos(u)
    f = lmb01*(bx**2+by**2)+lmb02*(bz**2)+lmba1*((az**2)-(1/3))*(bx**2+by**2)+lmba2*((az**2)-(1/3))*(bz**2)+lmbg1*(0.5*(ax**2-ay**2)*(bx**2-by**2)+ax*ay*bx*by)+lmbg2*(ax*az*bx*bz+ay*az*by*bz)+lmb12*(0.5*ay*az*(bx**2-by**2)+ax*az*bx*by)+lmb21*(0.5*(ax**2-ay**2)*by*bz+ax*ay*bx*bz)
    x = d*(1.0+s*f)*bx
    y = d*(1.0+s*f)*by
    z = d*(1.0+s*f)*bz
    dl_l=(np.sqrt(x**2+y**2+z**2)-d)/d

    hmod=np.sqrt(hx**2.0+hy**2.0+hz**2.0)+10.0**(-8)

    hhx=hx/hmod
    hhy=hy/hmod
    hhz=hz/hmod

    lsim = np.sqrt(x**2+y**2+z**2)
    list00 = np.stack((np.transpose(x),np.transpose(y),np.transpose(z),np.transpose(dl_l)),axis=0)
    list0 = np.stack((lsim,dl_l,dl_l/s,bx,by,bz,u*(180.0/np.pi),v*(180.0/np.pi)),axis=-1)
    hax=heffx-hx
    hay=heffy-hy
    haz=heffz-hz
    hamod=np.sqrt(hax**2.0+hay**2.0+haz**2.0)+10.0**(-12)
    hhax=hax/hamod
    hhay=hay/hamod
    hhaz=haz/hamod
    heffmod=np.sqrt(heffx**2.0+heffy**2.0+heffz**2.0)+10.0**(-12)
    hheffx=heffx/heffmod
    hheffy=heffy/heffmod
    hheffz=heffz/heffmod


    trfig = make_subplots(rows=1, cols=2,
                    specs=[[{'is_3d': True}, {'is_3d': True}]],
                    subplot_titles=['  Fields (H,Ha,Heff), magnetization (\u03B1) and unit cell lattice vectors (a,b,c)           ',
                                    '        Color corresponds to (\u0394l/lo)*scale_factor along direction \u03B2(θ,φ)'],)

    trfig.add_trace(go.Cone(x=[0], y=[1], z=[1], u=[hhx], v=[hhy], w=[hhz],name="""<b>H<b>/|<b>H<b>|""",
                          hovertemplate = """Hx/H = %{u:.6g}<br>Hy/H = %{v:.6g}<br>Hz/H = %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(255,0,0)']]),1, 1)

    trfig.add_trace(go.Cone(x=[1], y=[0], z=[1], u=[hhax], v=[hhay], w=[hhaz],name="""<b>Ha<b>/|<b>Ha<b>|""",
                          hovertemplate = """Hax/Ha = %{u:.6g}<br>Hay/Ha = %{v:.6g}<br>Haz/Ha = %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(225,0,0)']]),1, 1)

    trfig.add_trace(go.Cone(x=[1], y=[1], z=[1], u=[hheffx], v=[hheffy], w=[hheffz],name="""<b>Heff<b>/|<b>Heff<b>|""",
                          hovertemplate = """Heffx/Heff = %{u:.6g}<br>Heffy/Heff = %{v:.6g}<br>Heffz/Heff = %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(200,0,0)']]),1, 1)


    trfig.add_trace(go.Cone(x=[1], y=[1], z=[0], u=[ax], v=[ay], w=[az],name="""<b>\u03B1<b>""",
                          hovertemplate = """\u03B1x = %{u:.6g}<br>\u03B1y= %{v:.6g}<br>\u03B1z= %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(255,127,14)'], [1, 'rgb(255,127,16)']]),1, 1)
    trfig.add_trace(go.Cone(x=[0.4], y=[0], z=[0], u=[2], v=[0], w=[0],name="a",hoverinfo="name",
                          colorscale=[[0, 'rgb(255,0,0)'], [1, 'rgb(240,0,0)']]),1, 1)
    trfig.add_trace(go.Cone(x=[-0.2], y=[0.34641016], z=[0], u=[-1], v=[1.7320508],w=[0],name="b",hoverinfo="name",
                          colorscale=[[0, 'rgb(0,255,0)'], [1, 'rgb(0,240,0)']]),1, 1)
    trfig.add_trace(go.Cone(x=[0], y=[0], z=[0.4], u=[0], v=[0], w=[2],name="c",hoverinfo="name",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(0,0,240)']]),1, 1)

    trfig.update_traces(showscale=False)

    trfig.add_trace(go.Surface(x=list00[0], y=list00[1], z=list00[2],
                             surfacecolor=list00[3], customdata=list0, name=" ",hoverinfo="name",
                             hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 2)


    trfig.update_layout(transition_duration=500)



    trfig.update_yaxes(automargin=True)
    trfig.update_xaxes(automargin=True)

    return trfig,ax,ay,az,heffx,heffy,heffz,torque,hax,hay,haz,stat0





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
     Input(component_id='k1tr', component_property='value'),
     Input(component_id='k2tr', component_property='value'),
     Input(component_id='mstr', component_property='value'),
     Input(component_id='alpha0tr', component_property='value'),
     Input(component_id='toltr', component_property='value'),
     Input(component_id='dttr', component_property='value'),
     Input(component_id='nttr', component_property='value'),
    ]
)


def update_figtr2d(hx,hy,hz,lmb01,lmb02,lmba1,lmba2,lmbg1,lmbg2,lmb12,lmb21,s,kkk1,kkk2,mmms,alph,tol00,dtt,ntt):

    d=1.0

    crys = 'uni'

    ax,ay,az,stat0=llg(crys,mmms,hx,hy,hz,kkk1,kkk2,alph,tol00,dtt,ntt)

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


    dl_lxy = (np.sqrt(xxy**2+yxy**2+zxy**2)-d)/d
    lsimxy = np.sqrt(xxy**2+yxy**2+zxy**2)
    list00xy = np.stack((np.transpose(xxy),np.transpose(yxy),np.transpose(zxy)),axis=0)
    list0xy = np.stack((lsimxy,dl_lxy,dl_lxy/s,xxy/lsimxy,yxy/lsimxy,zxy/lsimxy, vxy*0.0,vxy*(180.0/np.pi)),axis=-1)

    dl_lxz = (np.sqrt(xxz**2+yxz**2+zxz**2)-d)/d
    lsimxz = np.sqrt(xxz**2+yxz**2+zxz**2)
    list00xz = np.stack((np.transpose(xxz),np.transpose(yxz),np.transpose(zxz)),axis=0)
    list0xz = np.stack((lsimxz,dl_lxz,dl_lxz/s,xxz/lsimxz,yxz/lsimxz,zxz/lsimxz,vxz*(180.0/np.pi),vxz*0.0),axis=-1)

    dl_lyz = (np.sqrt(xyz**2+yyz**2+zyz**2)-d)/d
    lsimyz = np.sqrt(xyz**2+yyz**2+zyz**2)
    uyz = np.mgrid[90.0:90.0:200j]
    list00yz = np.stack((np.transpose(xyz),np.transpose(yyz),np.transpose(zyz)),axis=0)
    list0yz = np.stack((lsimyz,dl_lyz,dl_lyz/s,xyz/lsimyz,yyz/lsimyz,zyz/lsimyz,vyz*(180.0/np.pi),uyz),axis=-1)

    figtr2d = make_subplots(rows=1, cols=3,
                    specs=[[{'type': 'xy'}, {'type': 'xy'}, {'type': 'xy'}]],
                    subplot_titles=['Plot XY plane (z=0)', 'Plot XZ plane (y=0)', 'Plot YZ plane (x=0)'],)

    figtr2d.add_trace(go.Scatter(x=bxxy, y=byxy, mode='lines',name='(lox,loy)'), 1, 1)
    figtr2d.add_trace(go.Scatter(x=list00xy[0], y=list00xy[1], mode='lines', customdata=list0xy, name="(lsimx,lsimy)", hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 1)

    figtr2d.add_trace(go.Scatter(x=bxxz, y=bzxz, mode='lines',name='(lox,loz)'), 1, 2)
    figtr2d.add_trace(go.Scatter(x=list00xz[0], y=list00xz[2], mode='lines', customdata=list0xz, name="(lsimx,lsimz)", hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 2)

    figtr2d.add_trace(go.Scatter(x=byyz, y=bzyz, mode='lines',name='(loy,loz)'), 1, 3)
    figtr2d.add_trace(go.Scatter(x=list00yz[1], y=list00yz[2], mode='lines', customdata=list0yz, name="(lsimy,lsimz)", hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 3)

    figtr2d.update_layout(xaxis_title='X-axis')
    figtr2d.update_layout(yaxis_title='Y-axis')
    figtr2d.update_layout(transition_duration=500)

    figtr2d.update_yaxes(automargin=True)
    figtr2d.update_xaxes(automargin=True)


    return figtr2d

##############MAE tri -3D

@app.callback(
    Output('mae_tri_3D', 'figure'),
    [Input(component_id='k0tr', component_property='value'),
     Input(component_id='k1tr', component_property='value'),
     Input(component_id='k2tr', component_property='value'),
    ]
)


def update_figmaetr(kkkk0,kkkk1,kkkk2):


    kkkk0=kkkk0*10**3
    kkkk1=kkkk1*10**3
    kkkk2=kkkk2*10**3


    u, v = np.mgrid[0:np.pi:100j, 0:2*np.pi:100j]
    bx = np.sin(u)*np.cos(v)
    by = np.sin(u)*np.sin(v)
    bz = np.cos(u)
    f = kkkk0+kkkk1*(1.0-bz**2.0)+kkkk2*(1.0-bz**2.0)**2.0
    f0 = f-kkkk0
    f0min = np.amin(f0)
    f0max = np.amax(f0)
    ff0 = f0 + np.absolute(f0min) + (f0max-f0min)*0.45

    x = ff0*bx
    y = ff0*by
    z = ff0*bz


    list00 = np.stack((np.transpose(x),np.transpose(y),np.transpose(z),np.transpose(f)),axis=0)
    list0 = np.stack((f,bx,by,bz,u*(180.0/np.pi),v*(180.0/np.pi)),axis=-1)


    figmaetr = make_subplots(rows=1, cols=1,
                    specs=[[{'is_3d': True}]],
                    subplot_titles=['  Color corresponds to Magnetocrystalline Anisotropy Energy (KJ/m^3) '])


    figmaetr.add_trace(go.Surface(x=list00[0], y=list00[1], z=list00[2], surfacecolor=list00[3], customdata=list0, name=" ", hoverinfo="name",
                                 hovertemplate = """Ea = %{customdata[0]:.6g} KJ/m^3<br>\u03B1x = %{customdata[1]:.6g}<br>\u03B1y = %{customdata[2]:.6g}<br>\u03B1z = %{customdata[3]:.6g}<br>θ = %{customdata[4]:.6g}°<br>φ = %{customdata[5]:.6g}°<br>"""), 1, 1)


    figmaetr.update_layout(transition_duration=500)

    figmaetr.update_layout(
        scene = {
            "xaxis": {"showticklabels": False},
            "yaxis": {"showticklabels": False},
            "zaxis": {"showticklabels": False}
        })

    return figmaetr



############## tet-3D



@app.callback(
    [Output('tet_3D', 'figure'),
     Output('meqxte', 'children'),
     Output('meqyte', 'children'),
     Output('meqzte', 'children'),
     Output('hefxte', 'children'),
     Output('hefyte', 'children'),
     Output('hefzte', 'children'),
     Output('torquete', 'children'),
     Output('haxte', 'children'),
     Output('hayte', 'children'),
     Output('hazte', 'children'),
     Output('statte', 'children'),
     ],
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
     Input(component_id='k1te', component_property='value'),
     Input(component_id='k2te', component_property='value'),
     Input(component_id='mste', component_property='value'),
     Input(component_id='alpha0te', component_property='value'),
     Input(component_id='tolte', component_property='value'),
     Input(component_id='dtte', component_property='value'),
     Input(component_id='ntte', component_property='value'),

    ]
)



def update_tefig(hx,hy,hz,lmb01,lmb02,lmba1,lmba2,lmbg,lmbd,lmbe,s,kk1,kk2,mms,alph,tol00,dtt,ntt):

    d=1.0

    crys = 'uni'

    ax,ay,az,stat0=llg(crys,mms,hx,hy,hz,kk1,kk2,alph,tol00,dtt,ntt)
    heffx,heffy,heffz=field(crys,ax,ay,az,mms,hx,hy,hz,kk1,kk2)
    torque=np.sqrt((heffy*az-heffz*ay)**2.0+(heffz*ax-heffx*az)**2.0+(heffx*ay-heffy*ax)**2.0)



    u, v = np.mgrid[0:np.pi:100j, 0:2*np.pi:100j]
    bx = np.sin(u)*np.cos(v)
    by = np.sin(u)*np.sin(v)
    bz = np.cos(u)
    f = lmb01*(bx**2+by**2)+lmb02*(bz**2)+lmba1*((az**2)-(1/3))*(bx**2+by**2)+lmba2*((az**2)-(1/3))*(bz**2)+lmbg*0.5*(ax**2-ay**2)*(bx**2-by**2)+lmbd*2*ax*ay*bx*by+2*lmbe*(ax*az*bx*bz+az*ay*bz*by)
    x = d*(1.0+s*f)*bx
    y = d*(1.0+s*f)*by
    z = d*(1.0+s*f)*bz
    dl_l=(np.sqrt(x**2+y**2+z**2)-d)/d

    hmod=np.sqrt(hx**2.0+hy**2.0+hz**2.0)+10.0**(-8)

    hhx=hx/hmod
    hhy=hy/hmod
    hhz=hz/hmod


    lsim = np.sqrt(x**2+y**2+z**2)
    list00 = np.stack((np.transpose(x),np.transpose(y),np.transpose(z),np.transpose(dl_l)),axis=0)
    list0 = np.stack((lsim,dl_l,dl_l/s,bx,by,bz,u*(180.0/np.pi),v*(180.0/np.pi)),axis=-1)
    hax=heffx-hx
    hay=heffy-hy
    haz=heffz-hz
    hamod=np.sqrt(hax**2.0+hay**2.0+haz**2.0)+10.0**(-12)
    hhax=hax/hamod
    hhay=hay/hamod
    hhaz=haz/hamod
    heffmod=np.sqrt(heffx**2.0+heffy**2.0+heffz**2.0)+10.0**(-12)
    hheffx=heffx/heffmod
    hheffy=heffy/heffmod
    hheffz=heffz/heffmod


    tefig = make_subplots(rows=1, cols=2,
                    specs=[[{'is_3d': True}, {'is_3d': True}]],
                    subplot_titles=['  Fields (H,Ha,Heff), magnetization (\u03B1) and unit cell lattice vectors (a,b,c)           ',
                                    '        Color corresponds to (\u0394l/lo)*scale_factor along direction \u03B2(θ,φ)'],)

    tefig.add_trace(go.Cone(x=[0], y=[1], z=[1], u=[hhx], v=[hhy], w=[hhz],name="""<b>H<b>/|<b>H<b>|""",
                          hovertemplate = """Hx/H = %{u:.6g}<br>Hy/H = %{v:.6g}<br>Hz/H = %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(255,0,0)']]),1, 1)

    tefig.add_trace(go.Cone(x=[1], y=[0], z=[1], u=[hhax], v=[hhay], w=[hhaz],name="""<b>Ha<b>/|<b>Ha<b>|""",
                          hovertemplate = """Hax/Ha = %{u:.6g}<br>Hay/Ha = %{v:.6g}<br>Haz/Ha = %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(225,0,0)']]),1, 1)

    tefig.add_trace(go.Cone(x=[1], y=[1], z=[1], u=[hheffx], v=[hheffy], w=[hheffz],name="""<b>Heff<b>/|<b>Heff<b>|""",
                          hovertemplate = """Heffx/Heff = %{u:.6g}<br>Heffy/Heff = %{v:.6g}<br>Heffz/Heff = %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(200,0,0)']]),1, 1)


    tefig.add_trace(go.Cone(x=[1], y=[1], z=[0], u=[ax], v=[ay], w=[az],name="""<b>\u03B1<b>""",
                          hovertemplate = """\u03B1x = %{u:.6g}<br>\u03B1y= %{v:.6g}<br>\u03B1z= %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(255,127,14)'], [1, 'rgb(255,127,16)']]),1, 1)
    tefig.add_trace(go.Cone(x=[0.4], y=[0], z=[0], u=[2], v=[0], w=[0],name="a",hoverinfo="name",
                          colorscale=[[0, 'rgb(255,0,0)'], [1, 'rgb(240,0,0)']]),1, 1)
    tefig.add_trace(go.Cone(x=[0.0], y=[0.4], z=[0], u=[0], v=[2], w=[0],name="b",hoverinfo="name",
                          colorscale=[[0, 'rgb(0,255,0)'], [1, 'rgb(0,240,0)']]),1, 1)
    tefig.add_trace(go.Cone(x=[0], y=[0], z=[0.4], u=[0], v=[0], w=[2],name="c",hoverinfo="name",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(0,0,240)']]),1, 1)

    tefig.update_traces(showscale=False)

    tefig.add_trace(go.Surface(x=list00[0], y=list00[1], z=list00[2],
                             surfacecolor=list00[3], customdata=list0, name=" ",hoverinfo="name",
                             hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 2)


    tefig.update_layout(transition_duration=500)



    tefig.update_yaxes(automargin=True)
    tefig.update_xaxes(automargin=True)

    return tefig,ax,ay,az,heffx,heffy,heffz,torque,hax,hay,haz,stat0




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
     Input(component_id='k1te', component_property='value'),
     Input(component_id='k2te', component_property='value'),
     Input(component_id='mste', component_property='value'),
     Input(component_id='alpha0te', component_property='value'),
     Input(component_id='tolte', component_property='value'),
     Input(component_id='dtte', component_property='value'),
     Input(component_id='ntte', component_property='value'),
    ]
)


def update_figte2d(hx,hy,hz,lmb01,lmb02,lmba1,lmba2,lmbg,lmbd,lmbe,s,kkk1,kkk2,mmms,alph,tol00,dtt,ntt):

    d=1.0

    crys = 'uni'

    ax,ay,az,stat0=llg(crys,mmms,hx,hy,hz,kkk1,kkk2,alph,tol00,dtt,ntt)

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


    dl_lxy = (np.sqrt(xxy**2+yxy**2+zxy**2)-d)/d
    lsimxy = np.sqrt(xxy**2+yxy**2+zxy**2)
    list00xy = np.stack((np.transpose(xxy),np.transpose(yxy),np.transpose(zxy)),axis=0)
    list0xy = np.stack((lsimxy,dl_lxy,dl_lxy/s,xxy/lsimxy,yxy/lsimxy,zxy/lsimxy, vxy*0.0,vxy*(180.0/np.pi)),axis=-1)

    dl_lxz = (np.sqrt(xxz**2+yxz**2+zxz**2)-d)/d
    lsimxz = np.sqrt(xxz**2+yxz**2+zxz**2)
    list00xz = np.stack((np.transpose(xxz),np.transpose(yxz),np.transpose(zxz)),axis=0)
    list0xz = np.stack((lsimxz,dl_lxz,dl_lxz/s,xxz/lsimxz,yxz/lsimxz,zxz/lsimxz,vxz*(180.0/np.pi),vxz*0.0),axis=-1)

    dl_lyz = (np.sqrt(xyz**2+yyz**2+zyz**2)-d)/d
    lsimyz = np.sqrt(xyz**2+yyz**2+zyz**2)
    uyz = np.mgrid[90.0:90.0:200j]
    list00yz = np.stack((np.transpose(xyz),np.transpose(yyz),np.transpose(zyz)),axis=0)
    list0yz = np.stack((lsimyz,dl_lyz,dl_lyz/s,xyz/lsimyz,yyz/lsimyz,zyz/lsimyz,vyz*(180.0/np.pi),uyz),axis=-1)

    figte2d = make_subplots(rows=1, cols=3,
                    specs=[[{'type': 'xy'}, {'type': 'xy'}, {'type': 'xy'}]],
                    subplot_titles=['Plot XY plane (z=0)', 'Plot XZ plane (y=0)', 'Plot YZ plane (x=0)'],)

    figte2d.add_trace(go.Scatter(x=bxxy, y=byxy, mode='lines',name='(lox,loy)'), 1, 1)
    figte2d.add_trace(go.Scatter(x=list00xy[0], y=list00xy[1], mode='lines', customdata=list0xy, name="(lsimx,lsimy)", hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 1)

    figte2d.add_trace(go.Scatter(x=bxxz, y=bzxz, mode='lines',name='(lox,loz)'), 1, 2)
    figte2d.add_trace(go.Scatter(x=list00xz[0], y=list00xz[2], mode='lines', customdata=list0xz, name="(lsimx,lsimz)", hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 2)

    figte2d.add_trace(go.Scatter(x=byyz, y=bzyz, mode='lines',name='(loy,loz)'), 1, 3)
    figte2d.add_trace(go.Scatter(x=list00yz[1], y=list00yz[2], mode='lines', customdata=list0yz, name="(lsimy,lsimz)", hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 3)

    figte2d.update_layout(xaxis_title='X-axis')
    figte2d.update_layout(yaxis_title='Y-axis')
    figte2d.update_layout(transition_duration=500)

    figte2d.update_yaxes(automargin=True)
    figte2d.update_xaxes(automargin=True)

    return figte2d


##############MAE tet -3D

@app.callback(
    Output('mae_tet_3D', 'figure'),
    [Input(component_id='k0te', component_property='value'),
     Input(component_id='k1te', component_property='value'),
     Input(component_id='k2te', component_property='value'),
    ]
)


def update_figmaete(kkkk0,kkkk1,kkkk2):


    kkkk0=kkkk0*10**3
    kkkk1=kkkk1*10**3
    kkkk2=kkkk2*10**3


    u, v = np.mgrid[0:np.pi:100j, 0:2*np.pi:100j]
    bx = np.sin(u)*np.cos(v)
    by = np.sin(u)*np.sin(v)
    bz = np.cos(u)
    f = kkkk0+kkkk1*(1.0-bz**2.0)+kkkk2*(1.0-bz**2.0)**2.0
    f0 = f-kkkk0
    f0min = np.amin(f0)
    f0max = np.amax(f0)
    ff0 = f0 + np.absolute(f0min) + (f0max-f0min)*0.45

    x = ff0*bx
    y = ff0*by
    z = ff0*bz


    list00 = np.stack((np.transpose(x),np.transpose(y),np.transpose(z),np.transpose(f)),axis=0)
    list0 = np.stack((f,bx,by,bz,u*(180.0/np.pi),v*(180.0/np.pi)),axis=-1)


    figmaete = make_subplots(rows=1, cols=1,
                    specs=[[{'is_3d': True}]],
                    subplot_titles=['  Color corresponds to Magnetocrystalline Anisotropy Energy (KJ/m^3) '])


    figmaete.add_trace(go.Surface(x=list00[0], y=list00[1], z=list00[2], surfacecolor=list00[3], customdata=list0, name=" ", hoverinfo="name",
                                 hovertemplate = """Ea = %{customdata[0]:.6g} KJ/m^3<br>\u03B1x = %{customdata[1]:.6g}<br>\u03B1y = %{customdata[2]:.6g}<br>\u03B1z = %{customdata[3]:.6g}<br>θ = %{customdata[4]:.6g}°<br>φ = %{customdata[5]:.6g}°<br>"""), 1, 1)


    figmaete.update_layout(transition_duration=500)

    figmaete.update_layout(
        scene = {
            "xaxis": {"showticklabels": False},
            "yaxis": {"showticklabels": False},
            "zaxis": {"showticklabels": False}
        })


    return figmaete




############## Ort-3D


@app.callback(
    [Output('ort_3D', 'figure'),
     Output('meqxo', 'children'),
     Output('meqyo', 'children'),
     Output('meqzo', 'children'),
     Output('hefxo', 'children'),
     Output('hefyo', 'children'),
     Output('hefzo', 'children'),
     Output('torqueo', 'children'),
     Output('haxo', 'children'),
     Output('hayo', 'children'),
     Output('hazo', 'children'),
     Output('stato', 'children'),
     ],
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
     Input(component_id='k1o', component_property='value'),
     Input(component_id='k2o', component_property='value'),
     Input(component_id='mso', component_property='value'),
     Input(component_id='alpha0o', component_property='value'),
     Input(component_id='tolo', component_property='value'),
     Input(component_id='dto', component_property='value'),
     Input(component_id='nto', component_property='value'),

    ]
)



def update_ofig(hx,hy,hz,lmb01,lmb02,lmb03,lmb1,lmb2,lmb3,lmb4,lmb5,lmb6,lmb7,lmb8,lmb9,s,kk1,kk2,mms,alph,tol00,dtt,ntt):

    d=1.0

    crys = 'ort'

    ax,ay,az,stat0=llg(crys,mms,hx,hy,hz,kk1,kk2,alph,tol00,dtt,ntt)
    heffx,heffy,heffz=field(crys,ax,ay,az,mms,hx,hy,hz,kk1,kk2)
    torque=np.sqrt((heffy*az-heffz*ay)**2.0+(heffz*ax-heffx*az)**2.0+(heffx*ay-heffy*ax)**2.0)



    u, v = np.mgrid[0:np.pi:100j, 0:2*np.pi:100j]
    bx = np.sin(u)*np.cos(v)
    by = np.sin(u)*np.sin(v)
    bz = np.cos(u)
    f = lmb01*bx**2+lmb02*by**2+lmb03*bz**2+lmb1*(ax**2*bx**2-ax*ay*bx*by-ax*az*bx*bz)+lmb2*(ay**2*bx**2-ax*ay*bx*by)+lmb3*(ax**2*by**2-ax*ay*bx*by)+lmb4*(ay**2*by**2-ax*ay*bx*by-ay*az*by*bz)+lmb5*(ax**2*bz**2-ax*az*bx*bz)+lmb6*(ay**2*bz**2-ay*az*by*bz)+lmb7*4*ax*ay*bx*by+lmb8*4*ax*az*bx*bz+lmb9*4*ay*az*by*bz


    x = d*(1.0+s*f)*bx
    y = d*(1.0+s*f)*by
    z = d*(1.0+s*f)*bz
    dl_l=(np.sqrt(x**2+y**2+z**2)-d)/d

    hmod=np.sqrt(hx**2.0+hy**2.0+hz**2.0)+10.0**(-8)

    hhx=hx/hmod
    hhy=hy/hmod
    hhz=hz/hmod


    lsim = np.sqrt(x**2+y**2+z**2)
    list00 = np.stack((np.transpose(x),np.transpose(y),np.transpose(z),np.transpose(dl_l)),axis=0)
    list0 = np.stack((lsim,dl_l,dl_l/s,bx,by,bz,u*(180.0/np.pi),v*(180.0/np.pi)),axis=-1)
    hax=heffx-hx
    hay=heffy-hy
    haz=heffz-hz
    hamod=np.sqrt(hax**2.0+hay**2.0+haz**2.0)+10.0**(-12)
    hhax=hax/hamod
    hhay=hay/hamod
    hhaz=haz/hamod
    heffmod=np.sqrt(heffx**2.0+heffy**2.0+heffz**2.0)+10.0**(-12)
    hheffx=heffx/heffmod
    hheffy=heffy/heffmod
    hheffz=heffz/heffmod


    ofig = make_subplots(rows=1, cols=2,
                    specs=[[{'is_3d': True}, {'is_3d': True}]],
                    subplot_titles=['  Fields (H,Ha,Heff), magnetization (\u03B1) and unit cell lattice vectors (a,b,c)           ',
                                    '        Color corresponds to (\u0394l/lo)*scale_factor along direction \u03B2(θ,φ)'],)

    ofig.add_trace(go.Cone(x=[0], y=[1], z=[1], u=[hhx], v=[hhy], w=[hhz],name="""<b>H<b>/|<b>H<b>|""",
                          hovertemplate = """Hx/H = %{u:.6g}<br>Hy/H = %{v:.6g}<br>Hz/H = %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(255,0,0)']]),1, 1)

    ofig.add_trace(go.Cone(x=[1], y=[0], z=[1], u=[hhax], v=[hhay], w=[hhaz],name="""<b>Ha<b>/|<b>Ha<b>|""",
                          hovertemplate = """Hax/Ha = %{u:.6g}<br>Hay/Ha = %{v:.6g}<br>Haz/Ha = %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(225,0,0)']]),1, 1)

    ofig.add_trace(go.Cone(x=[1], y=[1], z=[1], u=[hheffx], v=[hheffy], w=[hheffz],name="""<b>Heff<b>/|<b>Heff<b>|""",
                          hovertemplate = """Heffx/Heff = %{u:.6g}<br>Heffy/Heff = %{v:.6g}<br>Heffz/Heff = %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(200,0,0)']]),1, 1)


    ofig.add_trace(go.Cone(x=[1], y=[1], z=[0], u=[ax], v=[ay], w=[az],name="""<b>\u03B1<b>""",
                          hovertemplate = """\u03B1x = %{u:.6g}<br>\u03B1y= %{v:.6g}<br>\u03B1z= %{w:.6g}<br>""",
                          colorscale=[[0, 'rgb(255,127,14)'], [1, 'rgb(255,127,16)']]),1, 1)
    ofig.add_trace(go.Cone(x=[0.4], y=[0], z=[0], u=[2], v=[0], w=[0],name="a",hoverinfo="name",
                          colorscale=[[0, 'rgb(255,0,0)'], [1, 'rgb(240,0,0)']]),1, 1)
    ofig.add_trace(go.Cone(x=[0.0], y=[0.4], z=[0], u=[0], v=[2], w=[0],name="b",hoverinfo="name",
                          colorscale=[[0, 'rgb(0,255,0)'], [1, 'rgb(0,240,0)']]),1, 1)
    ofig.add_trace(go.Cone(x=[0], y=[0], z=[0.4], u=[0], v=[0], w=[2],name="c",hoverinfo="name",
                          colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(0,0,240)']]),1, 1)

    ofig.update_traces(showscale=False)

    ofig.add_trace(go.Surface(x=list00[0], y=list00[1], z=list00[2],
                             surfacecolor=list00[3], customdata=list0, name=" ",hoverinfo="name",
                             hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 2)


    ofig.update_layout(transition_duration=500)



    ofig.update_yaxes(automargin=True)
    ofig.update_xaxes(automargin=True)

    return ofig,ax,ay,az,heffx,heffy,heffz,torque,hax,hay,haz,stat0




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
     Input(component_id='k1o', component_property='value'),
     Input(component_id='k2o', component_property='value'),
     Input(component_id='mso', component_property='value'),
     Input(component_id='alpha0o', component_property='value'),
     Input(component_id='tolo', component_property='value'),
     Input(component_id='dto', component_property='value'),
     Input(component_id='nto', component_property='value'),
    ]
)


def update_figo2d(hx,hy,hz,lmb01,lmb02,lmb03,lmb1,lmb2,lmb3,lmb4,lmb5,lmb6,lmb7,lmb8,lmb9,s,kkk1,kkk2,mmms,alph,tol00,dtt,ntt):

    d=1.0

    crys = 'ort'

    ax,ay,az,stat0=llg(crys,mmms,hx,hy,hz,kkk1,kkk2,alph,tol00,dtt,ntt)

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

    dl_lxy = (np.sqrt(xxy**2+yxy**2+zxy**2)-d)/d
    lsimxy = np.sqrt(xxy**2+yxy**2+zxy**2)
    list00xy = np.stack((np.transpose(xxy),np.transpose(yxy),np.transpose(zxy)),axis=0)
    list0xy = np.stack((lsimxy,dl_lxy,dl_lxy/s,xxy/lsimxy,yxy/lsimxy,zxy/lsimxy, vxy*0.0,vxy*(180.0/np.pi)),axis=-1)

    dl_lxz = (np.sqrt(xxz**2+yxz**2+zxz**2)-d)/d
    lsimxz = np.sqrt(xxz**2+yxz**2+zxz**2)
    list00xz = np.stack((np.transpose(xxz),np.transpose(yxz),np.transpose(zxz)),axis=0)
    list0xz = np.stack((lsimxz,dl_lxz,dl_lxz/s,xxz/lsimxz,yxz/lsimxz,zxz/lsimxz,vxz*(180.0/np.pi),vxz*0.0),axis=-1)

    dl_lyz = (np.sqrt(xyz**2+yyz**2+zyz**2)-d)/d
    lsimyz = np.sqrt(xyz**2+yyz**2+zyz**2)
    uyz = np.mgrid[90.0:90.0:200j]
    list00yz = np.stack((np.transpose(xyz),np.transpose(yyz),np.transpose(zyz)),axis=0)
    list0yz = np.stack((lsimyz,dl_lyz,dl_lyz/s,xyz/lsimyz,yyz/lsimyz,zyz/lsimyz,vyz*(180.0/np.pi),uyz),axis=-1)

    figo2d = make_subplots(rows=1, cols=3,
                    specs=[[{'type': 'xy'}, {'type': 'xy'}, {'type': 'xy'}]],
                    subplot_titles=['Plot XY plane (z=0)', 'Plot XZ plane (y=0)', 'Plot YZ plane (x=0)'],)

    figo2d.add_trace(go.Scatter(x=bxxy, y=byxy, mode='lines',name='(lox,loy)'), 1, 1)
    figo2d.add_trace(go.Scatter(x=list00xy[0], y=list00xy[1], mode='lines', customdata=list0xy, name="(lsimx,lsimy)", hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 1)

    figo2d.add_trace(go.Scatter(x=bxxz, y=bzxz, mode='lines',name='(lox,loz)'), 1, 2)
    figo2d.add_trace(go.Scatter(x=list00xz[0], y=list00xz[2], mode='lines', customdata=list0xz, name="(lsimx,lsimz)", hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 2)

    figo2d.add_trace(go.Scatter(x=byyz, y=bzyz, mode='lines',name='(loy,loz)'), 1, 3)
    figo2d.add_trace(go.Scatter(x=list00yz[1], y=list00yz[2], mode='lines', customdata=list0yz, name="(lsimy,lsimz)", hovertemplate = """lsim = %{customdata[0]:.6g} a.u.<br>(\u0394l/lo)*s = %{customdata[1]:.6g}<br>\u0394l/lo = %{customdata[2]:.6g}<br>\u03B2x = %{customdata[3]:.6g}<br>\u03B2y = %{customdata[4]:.6g}<br>\u03B2z = %{customdata[5]:.6g}<br>θ = %{customdata[6]:.6g}°<br>φ = %{customdata[7]:.6g}°<br>"""), 1, 3)

    figo2d.update_layout(xaxis_title='X-axis')
    figo2d.update_layout(yaxis_title='Y-axis')
    figo2d.update_layout(transition_duration=500)

    figo2d.update_yaxes(automargin=True)
    figo2d.update_xaxes(automargin=True)


    return figo2d


##############MAE ort -3D

@app.callback(
    Output('mae_ort_3D', 'figure'),
    [Input(component_id='k0o', component_property='value'),
     Input(component_id='k1o', component_property='value'),
     Input(component_id='k2o', component_property='value'),
    ]
)


def update_figmaeo(kkkk0,kkkk1,kkkk2):


    kkkk0=kkkk0*10**3
    kkkk1=kkkk1*10**3
    kkkk2=kkkk2*10**3


    u, v = np.mgrid[0:np.pi:100j, 0:2*np.pi:100j]
    bx = np.sin(u)*np.cos(v)
    by = np.sin(u)*np.sin(v)
    bz = np.cos(u)
    f = kkkk0+kkkk1*bx**2.0+kkkk2*by**2.0
    f0 = f-kkkk0
    f0min = np.amin(f0)
    f0max = np.amax(f0)
    ff0 = f0 + np.absolute(f0min) + (f0max-f0min)*0.45

    x = ff0*bx
    y = ff0*by
    z = ff0*bz


    list00 = np.stack((np.transpose(x),np.transpose(y),np.transpose(z),np.transpose(f)),axis=0)
    list0 = np.stack((f,bx,by,bz,u*(180.0/np.pi),v*(180.0/np.pi)),axis=-1)


    figmaeo = make_subplots(rows=1, cols=1,
                    specs=[[{'is_3d': True}]],
                    subplot_titles=['  Color corresponds to Magnetocrystalline Anisotropy Energy (KJ/m^3) '])


    figmaeo.add_trace(go.Surface(x=list00[0], y=list00[1], z=list00[2], surfacecolor=list00[3], customdata=list0, name=" ", hoverinfo="name",
                                 hovertemplate = """Ea = %{customdata[0]:.6g} KJ/m^3<br>\u03B1x = %{customdata[1]:.6g}<br>\u03B1y = %{customdata[2]:.6g}<br>\u03B1z = %{customdata[3]:.6g}<br>θ = %{customdata[4]:.6g}°<br>φ = %{customdata[5]:.6g}°<br>"""), 1, 1)


    figmaeo.update_layout(transition_duration=500)

    figmaeo.update_layout(
        scene = {
            "xaxis": {"showticklabels": False},
            "yaxis": {"showticklabels": False},
            "zaxis": {"showticklabels": False}
        })

    return figmaeo


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

    phix = 0.0*x
    phiy = 0.0*y
    phiz = phi*(z/length)
    phizd = phi*(z/length)*(180.0/np.pi)

    list00 = np.stack((np.transpose(x),np.transpose(y),np.transpose(z),np.transpose(phiz)),axis=0)
    list0 = np.stack((phix,phiy,phiz,phizd),axis=-1)


    figrod = make_subplots(rows=1, cols=1,
                    specs=[[{'is_3d': True}]],
                    subplot_titles=['Electric current I=Imax. Color of the rod corresponds to the twisted angle ϕ(z) in radian.'],)

    figrod.add_trace(go.Cone(x=[4*rad], y=[0], z=[0.5*length], u=[0], v=[htn], w=[hln],name="H∥+H⟂",sizemode="absolute",showscale=False,colorscale=[[0, 'rgb(0,255,0)'], [1, 'rgb(255,0,0)']]),1, 1)
    figrod.add_trace(go.Cone(x=[-4*rad], y=[0], z=[0.5*length], u=[0], v=[-htn], w=[hln],name="H∥+H⟂",sizemode="absolute",showscale=False,colorscale=[[0, 'rgb(0,255,0)'], [1, 'rgb(255,0,0)']]),1, 1)
    figrod.add_trace(go.Cone(x=[0], y=[4*rad], z=[0.5*length], u=[-htn], v=[0], w=[hln],name="H∥+H⟂",sizemode="absolute",showscale=False,colorscale=[[0, 'rgb(0,255,0)'], [1, 'rgb(255,0,0)']]),1, 1)
    figrod.add_trace(go.Cone(x=[0], y=[-4*rad], z=[0.5*length], u=[htn], v=[0], w=[hln],name="H∥+H⟂",sizemode="absolute",showscale=False,colorscale=[[0, 'rgb(0,255,0)'], [1, 'rgb(255,0,0)']]),1, 1)
    figrod.update_traces(hoverinfo="name")

    figrod.add_trace(go.Surface(x=list00[0], y=list00[1], z=list00[2], surfacecolor=list00[3], customdata=list0, name="magnetic rod",hoverinfo="name",
                               hovertemplate = """x = %{x:.6g} m<br>y = %{y:.6g} m<br>z = %{z:.6g} m<br>ϕ(z) = %{customdata[2]:.6g} rad<br>ϕ(z) = %{customdata[3]:.6g} °<br>"""), 1, 1)


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
