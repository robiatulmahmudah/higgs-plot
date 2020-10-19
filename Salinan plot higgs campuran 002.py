# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 16:18:34 2020

@author: TUF GAMING
"""

import matplotlib.pyplot as plt
import numpy as np


"definisi fungsi"

Mp=54000.0
alph= np.sqrt(2.0/3.0)*(1.0/Mp)
xi= 1.0*10000
lam=0.01
M=np.sqrt((4.0*np.pi)/3.0)*(Mp/xi)
ph0=1.2*Mp
vh0=0.0
hh0=Mp/np.sqrt(xi)
gh0=0.0
t0=0.0

def U(ph,hh):
    U1=lam*np.exp(-2.0*alph*ph)*(hh**4)
    U2=3.0*Mp*Mp*M*M/4.0
    U3=1.0-((1.0 + 1.0*xi*hh*hh/(Mp**2) )*np.exp(-alph*ph) )
    U= U1+(U2*(U3**2.0))
    return U


def Uph(ph,hh):
    u1=-2.0*alph*lam*np.exp(-2.0*alph*ph)*(hh**4)
    u2=1.5*M*M*Mp*Mp*alph*(1.0+(1.0*xi*hh*hh)/(Mp*Mp))*np.exp(-alph*ph)
    u3=-1.5*M*M*Mp*Mp*alph*((1.0+(1.0*xi*hh*hh)/(Mp*Mp))**2 )*np.exp(-2.0*alph*ph)
    U=u1+u2+u3
    return U

def Uhh(ph,hh):
    u1=4.0*lam*np.exp(-2.0*alph*ph)*hh*hh*hh
    u2=-6.0*M*M*xi*np.exp(-alph*ph)*hh
    u3=1.5*M*M*(1.0+(1.0*xi*hh*hh/(Mp*Mp) ) )*4.0*xi*hh*np.exp(-2.0*alph*ph)
    U=u1+u2+u3
    return U

def Hub(ph,vh,hh,gh):
    Hk=(0.5*vh*vh)+(0.5*gh*gh*np.exp(-alph*ph))+U(ph,hh)
    Hm=Hk/(3.0*Mp*Mp)
    H=np.sqrt(Hm)
    return H

"langsung definisikan parh dan parp pada dua fungsi di bawah ini"
parp0=Mp/(10**15)
parh0=Mp/(10**15)


def f1(ph,vh,hh,gh,t):
    f1= -(3.0*Hub(ph,vh,hh,gh)*vh)-(0.5*alph*np.exp(-alph*ph)*gh*gh)-Uph(ph,hh)
    return f1

def f2(ph,vh,hh,gh,t):
    f2= vh
    return f2

def g1(ph,vh,hh,gh,t):
    g1= -(3.0*Hub(ph,vh,hh,gh)*gh)+(alph*vh*gh)-(np.exp(alph*ph)*Uhh(ph,hh))
    return g1

"print(g1(ph0,vh0,hh0,gh0,t0))"

def g2(ph,vh,hh,gh,t):
    g2= gh
    return g2

def Tk(ph,hh):
    N3=24.0*np.exp(-1.62*ph/Mp)*M*hh/Mp*(xi*(0.5*np.exp(0.81*ph/Mp)-1-xi*hh/Mp**2)-lam*hh*hh/(3.0*M*M))
    N1=-2.0*np.exp(-1.62*ph/Mp)*M*M*((1+xi*h*h/(M*M))*(0.5*np.exp(0.81*ph/Mp)-1-xi*h*h/(Mp*Mp))-lam*hh**4/(3.0*M*M*Mp*Mp))
    N2=-3.0*np.exp(-1.62*ph/Mp)*M*M*(xi*(np.exp(0.81*ph/Mp)-1-3*xi*h*h/(Mp*Mp))+lam*hh*hh/(3*M*M))
    T= N1*N2-N3**2
    return T

"definisikan langkah iterasi"
h=Mp/(10.0**10)

"bagian runge kutta"
imax=10000000
i=0
sumbux=[t0]
sumbuy=[ph0]
sumbuz=[hh0]
while i<imax:
    f11=h*f1(ph0,vh0,hh0,gh0,t0)
    f21=h*f2(ph0,vh0,hh0,gh0,t0)
    g11=h*g1(ph0,vh0,hh0,gh0,t0)
    g21=h*g2(ph0,vh0,hh0,gh0,t0)
    
    "print('step 1',f11,f21,g11,g21,t0)"
    
    f12=h*f1(ph0+0.5*f21,vh0+0.5*f11,hh0+0.5*g21,gh0+0.5*g11,t0+0.5*h)
    f22=h*f2(ph0+0.5*f21,vh0+0.5*f11,hh0+0.5*g21,gh0+0.5*g11,t0+0.5*h)
    g12=h*g1(ph0+0.5*f21,vh0+0.5*f11,hh0+0.5*g21,gh0+0.5*g11,t0+0.5*h)
    g22=h*g2(ph0+0.5*f21,vh0+0.5*f11,hh0+0.5*g21,gh0+0.5*g11,t0+0.5*h)
    
    "print('step 2',f12,f22,g12,g22)"
    
    f13=h*f1(ph0+0.5*f22,vh0+0.5*f12,hh0+0.5*g22,gh0+0.5*g12,t0+0.5*h)
    f23=h*f2(ph0+0.5*f22,vh0+0.5*f12,hh0+0.5*g22,gh0+0.5*g12,t0+0.5*h)
    g13=h*g1(ph0+0.5*f22,vh0+0.5*f12,hh0+0.5*g22,gh0+0.5*g12,t0+0.5*h)
    g23=h*g2(ph0+0.5*f22,vh0+0.5*f12,hh0+0.5*g22,gh0+0.5*g12,t0+0.5*h) 
    
    "print('step 3',f13,f23,g13,g23)"
    
    f14=h*f1(ph0+f23,vh0+f13,hh0+g23,gh0+g13,t0+h)
    f24=h*f2(ph0+f23,vh0+f13,hh0+g23,gh0+g13,t0+h)
    g14=h*g1(ph0+f23,vh0+f13,hh0+g23,gh0+g13,t0+h)
    g24=h*g2(ph0+f23,vh0+f13,hh0+g23,gh0+g13,t0+h)
    
    "print('step 4',f14,f24,g14,g24)"
    
    phi=ph0+(f21+2.0*f22+2.0*f23+f24)/6.0
    vhi=vh0+(f11+2.0*f12+2.0*f13+f14)/6.0
    hhi=hh0+(g21+2.0*g22+2.0*g23+g24)/6.0
    ghi=gh0+(g11+2.0*g12+2.0*g13+g14)/6.0
    ti=t0+h
    "print(phi,vhi,hhi,ghi)"
    
    sumbux.append(ti)
    sumbuy.append(phi)
    sumbuz.append(hhi)
    
    
    ph0=phi
    vh0=vhi
    hh0=hhi
    gh0=ghi
    t0=ti
    i=i+1
    
plt.plot(sumbux,sumbuy)
plt.plot(sumbux,sumbuz)
plt.xlim(3,30)
A=np.array(sumbux)
B=np.array(sumbuy)
C=np.array(sumbuz)
sumbu=[]

def Tk(ph,hh):
    sumbux
    N3=24.0*np.exp(-1.62*ph/Mp)*M*hh/Mp*(xi*(0.5*np.exp(0.81*ph/Mp)-1-xi*hh/Mp**2)-lam*hh*hh/(3.0*M*M))
    N1=-2.0*np.exp(-1.62*ph/Mp)*M*M*((1+xi*h*h/(M*M))*(0.5*np.exp(0.81*ph/Mp)-1-xi*h*h/(Mp*Mp))-lam*hh**4/(3.0*M*M*Mp*Mp))
    N2=-3.0*np.exp(-1.62*ph/Mp)*M*M*(xi*(np.exp(0.81*ph/Mp)-1-3*xi*h*h/(Mp*Mp))+lam*hh*hh/(3*M*M))
    T= N1*N2-N3**2
    return T

for X, Y, Z in zip(A,B,C):
    sumbu.append([X,Y,Z])
    
outfile = open('hasil.dat', 'w')
for row in sumbu:
    for column in row:
        outfile.write('%12.9f' % column)
    outfile.write('\n')
outfile.close()