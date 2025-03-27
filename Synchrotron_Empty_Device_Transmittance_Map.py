# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 10:50:37 2023

@author: amart
"""

import numpy as np

import pyFAI

import pandas as pd

import seaborn as sns

import matplotlib.pyplot as plt
import matplotlib as mpl
import os 

import fabio

from scipy.integrate import simpson
from numpy import trapz
from sklearn.metrics import auc

from scipy.optimize import curve_fit

sns.set_theme()
sns.set_style("ticks")
sns.set_context("paper",font_scale=0.9)
sns.set_palette("tab10")


os.chdir(r"C:\Users\Alberto\OneDrive - UAM\Alberto\ALBA\ID2023027489_November_2023\DATA\Analysis\Absorcion_200um")


def mesh(positions):
    x = {}

    z = {}

    i = 0

    for sx in positions.sx:
        
        if str(round(sx,2)) not in x.keys() and sx>-5 and sx<5:
            x[str(round(sx,2))]=i
            i+=1
        
    i = 0
    for stz in positions.stz:
        
        if str(round(stz,2)) not in z.keys():
            z[str(round(stz,2))]=i
            i+=1


    files = [i for i in os.listdir(r"raw/mesh_snapdmesh_4065/data_00") if "pilatus" in i]


    mesh = np.zeros((len(z.keys()),len(x.keys())))

    for file in files:
        
        No = int(file[:-4].split("_")[-1])
        
        sx = positions.loc[positions.loc[:,"Pt_No"]==No,"sx"]
        
        sx = round(sx.iloc[0],2)
        
        if sx >-5 and sx<5:
            col = x[str(sx)]
            
            stz = positions.loc[positions.loc[:,"Pt_No"]==No,"stz"]
            
            stz = round(stz.iloc[0],2)
            
            row = z[str(stz)]
            
            
            
            img = fabio.open(f"raw/mesh_snapdmesh_4065/data_00/{file}")
            
            value = float(img.header["Photo"])/float(img.header["MonitorRing"])
            """
            if row!=0:
                for j in range(33*row,33*(row+1)):
                    mesh[j,col] = value
            
            else:
                for j in range(0,33*(row+1)):
                    mesh[j,col] = value
            """
            mesh[row,col] = value
        
    plt.imshow(mesh,aspect=3)    
    plt.colorbar()
    plt.axis("Off")
    plt.savefig("figures/map.tif",bbox_inches="tight",dpi=300)
    plt.show()
    plt.close()
    
    
positions = r"raw/mesh_snapdmesh_4065/mesh_snapdmesh_4065.dat"

positions = pd.read_table(positions,skiprows = 6,header= None,sep=" ",skipfooter=1,engine="python", names = ["Pt_No","sx",  "stz",  "uxtimer",  "pilatus_timer",  "rayonix_timer",  "pilatus_roi",  "rayonix_roi",  "tfg_timer",  "tfg_ch1",  "tfg_ch2",  "tfg_ch3",  "tfg_ch4", "current",  "dt"])

files = [i for i in os.listdir(r"raw/mesh_snapdmesh_4065/data_00") if "pilatus" in i]

z = {}


Nos = []
for Pt in positions.Pt_No[:-1]:

    if round(positions.stz[Pt],2) == round(positions.stz[Pt+1],2):
        Nos.append(Pt)
    
    
    else:
        Nos.append(Pt)
        z[round(positions.stz[Pt],2)]=Nos
        Nos=[]

Pt+=1
Nos.append(Pt)
z[round(positions.stz[Pt],2)]=Nos

x_mean = np.mean(positions.sx)
for key in z.keys():
   
    x = []
    
    I = []
    
    for No in z[key]:
        
        sx = positions.loc[No,"sx"]
    
        if sx>-8 and sx<8:
            x.append(sx-x_mean)
            
            file = [i for i in files if int(i[:-4].split("_")[-1])==No][0]
            
            img = fabio.open(f"raw/mesh_snapdmesh_4065/data_00/{file}")
        
            I.append(float(img.header["Photo"])/float(img.header["Monitor"]))
            
    
    plt.plot(x,I,label=str(key))

axis=plt.gca()
ylims = axis.get_ylim()

plt.fill_betweenx(ylims,(-8-x_mean),(-3.9),alpha=0.3,color="gray")
plt.fill_betweenx(ylims,(3.9),(8-x_mean),alpha=0.3,color="gray")
plt.fill_betweenx(ylims,(-3.9),(-1.1),alpha=0.3,color="darkorange")
plt.fill_betweenx(ylims,(1.1),(3.9),alpha=0.3,color="darkorange")
plt.fill_betweenx(ylims,(-1.1),(1.1),alpha=0.3,color="green")
plt.text(0,np.mean(ylims),"Window",ha="center",weight="bold",color="green")
plt.text(np.mean([-3.9,-1.1]),np.mean(ylims),"Channel",ha="center",weight="bold",color="darkorange")
plt.text(np.mean([3.9,1.1]),np.mean(ylims),"Channel",ha="center",weight="bold",color="darkorange")
plt.text(np.mean([3.9,8-x_mean]),np.mean(ylims),"Chip",ha="center",weight="bold",color="gray")
plt.text(np.mean([-3.9,-8-x_mean]),np.mean(ylims),"Chip",ha="center",weight="bold",color="gray")


#plt.legend()

plt.xlabel("x position (mm)")
plt.ylabel("Normalized Photodiode Current (a. u.)")
plt.savefig("figures/Profiles.tif",bbox_inches="tight",dpi=300)
plt.show()
plt.close()


positions = r"raw/mesh_snapdmesh_4065/mesh_snapdmesh_4065.dat"

positions = pd.read_table(positions,skiprows = 6,header= None,sep=" ",skipfooter=1,engine="python", names = ["Pt_No","sx",  "stz",  "uxtimer",  "pilatus_timer",  "rayonix_timer",  "pilatus_roi",  "rayonix_roi",  "tfg_timer",  "tfg_ch1",  "tfg_ch2",  "tfg_ch3",  "tfg_ch4", "current",  "dt"])

files = [i for i in os.listdir(r"raw/mesh_snapdmesh_4065/data_00") if "pilatus" in i]

z = {}


Nos = []
for Pt in positions.Pt_No[:-1]:

    if round(positions.stz[Pt],2) == round(positions.stz[Pt+1],2):
        Nos.append(Pt)
    
    
    else:
        Nos.append(Pt)
        z[round(positions.stz[Pt],2)]=Nos
        Nos=[]

Pt+=1
Nos.append(Pt)
z[round(positions.stz[Pt],2)]=Nos

x_mean = np.mean(positions.sx)

colors = ["b","darkorange","g"]
k = 0
for key in z.keys():
   
    x = []
    
    I = []
    
    for No in z[key]:
        
        sx = positions.loc[No,"sx"]
    
        if sx>-8 and sx<8:
            x.append(sx-x_mean)
            
            file = [i for i in files if int(i[:-4].split("_")[-1])==No][0]
            
            img = fabio.open(f"raw/mesh_snapdmesh_4065/data_00/{file}")
        
            I.append(float(img.header["Photo"])/float(img.header["Monitor"]))
            
    
    plt.plot(x,I,colors[k],label=str(key))
    
    axis=plt.gca()
    
    
    
    ylims = axis.get_ylim()
    
    if key==-28.12:
        plt.fill_betweenx(ylims,(-3.9),(-1.1),alpha=0.3,color="darkorange")
        plt.fill_betweenx(ylims,(1.1),(3.9),alpha=0.3,color="darkorange")
        plt.fill_betweenx(ylims,(-1.1),(1.1),alpha=0.3,color="green")
        plt.fill_betweenx(ylims,(-8-x_mean),(-3.9),alpha=0.3,color="gray")
        plt.fill_betweenx(ylims,(3.9),(8-x_mean),alpha=0.3,color="gray")
        
        plt.text(0,np.mean(ylims),"Window",ha="center",weight="bold",color="green")
        plt.text(np.mean([-3.9,-1.1]),np.mean(ylims),"Channel",ha="center",weight="bold",color="darkorange")
        plt.text(np.mean([3.9,1.1]),np.mean(ylims),"Channel",ha="center",weight="bold",color="darkorange")
        plt.text(np.mean([3.9,8-x_mean]),np.mean(ylims),"Chip",ha="center",weight="bold",color="gray")
        plt.text(np.mean([-3.9,-8-x_mean]),np.mean(ylims),"Chip",ha="center",weight="bold",color="gray")
    
    elif key==-27.12:
        
        plt.fill_betweenx(ylims,(-3.5),(-1.1),alpha=0.3,color="darkorange")
        plt.fill_betweenx(ylims,(1.1),(3.5),alpha=0.3,color="darkorange")
        plt.fill_betweenx(ylims,(-1.1),(1.1),alpha=0.3,color="green")
        plt.fill_betweenx(ylims,(-8-x_mean),(-3.5),alpha=0.3,color="gray")
        plt.fill_betweenx(ylims,(3.5),(8-x_mean),alpha=0.3,color="gray")
        
        plt.text(0,np.mean(ylims),"Window",ha="center",weight="bold",color="green")
        plt.text(np.mean([-3.5,-1.1]),np.mean(ylims),"Channel",ha="center",weight="bold",color="darkorange")
        plt.text(np.mean([3.5,1.1]),np.mean(ylims),"Channel",ha="center",weight="bold",color="darkorange")
        plt.text(np.mean([3.5,8-x_mean]),np.mean(ylims),"Chip",ha="center",weight="bold",color="gray")
        plt.text(np.mean([-3.5,-8-x_mean]),np.mean(ylims),"Chip",ha="center",weight="bold",color="gray")
    
    elif key==-27.62:
        
        plt.fill_betweenx(ylims,-0.2,0.2,alpha=0.3,color="cyan")
        
        plt.fill_betweenx(ylims,(-1.1),(-0.2),alpha=0.3,color="green")
        plt.fill_betweenx(ylims,(0.2),(1.1),alpha=0.3,color="green")
        plt.fill_betweenx(ylims,(-8-x_mean),(-1.1),alpha=0.3,color="gray")
        plt.fill_betweenx(ylims,(1.1),(8-x_mean),alpha=0.3,color="gray")
        
        plt.text(0,np.mean(ylims),"Window",ha="center",weight="bold",color="green")
        #plt.text(np.mean([-3.5,-1.1]),np.mean(ylims),"Channel",ha="center",weight="bold",color="darkorange")
        #plt.text(np.mean([3.5,1.1]),np.mean(ylims),"Channel",ha="center",weight="bold",color="darkorange")
        plt.text(np.mean([1.1,8-x_mean]),np.mean(ylims),"Chip",ha="center",weight="bold",color="gray")
        plt.text(np.mean([-1.1,-8-x_mean]),np.mean(ylims),"Chip",ha="center",weight="bold",color="gray")
    
    
    
    plt.xlabel("x position (mm)")
    plt.ylabel("Normalized Photodiode Current (a. u.)")
    plt.savefig(f"figures/{key}.tif",bbox_inches="tight",dpi=300)
    plt.show()
    plt.close()
    k+=1



mesh(positions)












