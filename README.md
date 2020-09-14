# MAELASviewer
Online visualization of magnetostriction


![GitHub Logo](assets/logo_maelasviewer.png)

MAELASviewer

Authors: P. Nieves, S. Arapan, S.H. Zhang, A.P. Kądzielawa, Lukáš Kývala, R.F. Zhang and D. Legut

-------------------------
WHAT IS MAELASviewer?
-------------------------

MAELASviewer is an online tool to visualize magnetostriction phenomena. This interactive applet shows the magnetostriction for some crystal systems. You can visualize the relative length change (Δl/lo=[l-lo]/lo) of the material along an arbitrary direction (β) as a function of the external magnetic field (H) and magnetostrictive coefficients (λ).

------------------
INSTALLATION
------------------

The MAELAS code requires to have ```Python3(>=3.6)```. For example, in Ubuntu Linux machine you can check the installed version of ```python3``` by opening a terminal and typing
```bash
python3 --version
```
In case you need to install python3 in your machine, you can type
```bash
sudo apt-get update
sudo apt-get install python3
```
Additionally, you need to install the following dependencies

```bash
dash(>=0.5.10.2)
plotly(>=4.7.1)
numpy(>=1.18.4)
```
you can easily install them with pip3
```bash
pip3 install dash
pip3 install plotly
pip3 isntall numpy
```
If you need to install pip3 on Ubuntu Linux, then type
```bash
sudo apt-get update
sudo apt-get install python3-pip
```
Next, download and extract the .zip file, go to the folder that contains the file ```maelasviewer.py``` and type
```bash
python3 maelasviewer.py
```
then visit http://127.0.0.1:8050/ in your web browser to use MAELASviewer.

