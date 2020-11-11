# MAELASviewer
Online visualization of magnetostriction


![GitHub Logo](assets/logo_maelasviewer.png)

MAELASviewer

Authors: P. Nieves, S. Arapan, A.P. Kądzielawa and D. Legut

-------------------------
WHAT IS MAELASviewer?
-------------------------

MAELASviewer is an online tool to visualize magnetostriction. This interactive applet shows the magnetostriction for some crystal systems. You can simulate the Joule and Wiedemann effects.

-------------------------
ONLINE ACCESS
-------------------------

The online application is available at 

'''[https://maelasviewer.herokuapp.com](https://maelasviewer.herokuapp.com)'''


--------------------------
OFFLINE AND DEVELOPER MODE
--------------------------

MAELASviewer can also be executed in your local computer. It requires to have ```Python3(>=3.6)```. For example, in Ubuntu Linux machine you can check the installed version of ```python3``` by opening a terminal and typing
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
then visit http://127.0.0.1:8050/ in your web browser to use MAELASviewer.

------------------------------
DOCUMENTATION
------------------------------

More details of this application can be found in:

Nieves, P.; Arapan, S.; Kądzielawa, A.P.; Legut, D. MAELASviewer: An Online Tool to Visualize Magnetostriction. Sensors 2020, 20, 6436. 

'''[https://www.mdpi.com/1424-8220/20/22/6436/pdf](https://www.mdpi.com/1424-8220/20/22/6436/pdf)'''

