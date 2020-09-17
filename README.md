# MAELASviewer
Online visualization of magnetostriction


![GitHub Logo](assets/logo_maelasviewer.png)

MAELASviewer

Authors: P. Nieves, S. Arapan, A.P. Kądzielawa and D. Legut

-------------------------
WHAT IS MAELASviewer?
-------------------------

MAELASviewer is an online tool to visualize magnetostriction. This interactive applet shows the magnetostriction for some crystal systems. You can visualize the relative length change of the material along an arbitrary direction as a function of the external magnetic field and magnetostrictive coefficients.

------------------------
ONLINE APP
------------------------

MAELASviewer can be used in the following link

https://maelasviewer.herokuapp.com/

-------------------------------
RUNNING MAELASviewer FROM SOURCE FILES 
-------------------------------

Alternatively, it is also possible to use MAELASviewer running the source files. It requires to have ```Python3(>=3.6)```. 
For example, in Ubuntu Linux machine you can check the installed version of ```python3``` by opening a terminal and typing
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
pip3 install numpy
```
If you need to install pip3 on Ubuntu Linux, then type
```bash
sudo apt-get update
sudo apt-get install python3-pip
```
Once all depencencies are installed, download and extract the MAELASviewer-master.zip file, go to the folder that contains the file ```maelasviewer.py``` and type
```bash
python3 maelasviewer.py
```
then visit http://127.0.0.1:8050/ in your web browser to use MAELASviewer. This will run the app on ```localhost```, it means only on your local machine (development server).

-------------------------------
MANUSCRIPT
-------------------------------

P. Nieves, S. Arapan, A.P. Kądzielawa and D. Legut, MAELASviewer: an online tool to visualize magnetostriction, 2020, arXiv