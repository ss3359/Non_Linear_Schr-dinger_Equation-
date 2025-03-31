

import math 
import random 
import numpy as np
import pandas as pd 
from mpmath import mp
from sympy import symbols,diff, Matrix
from sympy import re
from sympy import sqrt,acos,cos,sin,atan,tan
from sympy import Inverse,DotProduct, simplify
from sympy import Subs

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


'''
This is the non-linear Schrodinger Equations in terms of the real and imaginary components
Real: 
    Vt+(1/2)Uxx + (u^2 + v^2)u=0 => 
Imaginary: 
    -Ut+(1/2)Vxx+(u^2 + v^2)v=0 
Initial Conditions: 
    u(x,0)=f(x)
    v(x,0)=g(x)
Boundary Conditions: 
    u(0,t)=u(L,t)=0
    v(0,t)=v(L,t)=0
    
'''
'''
Vt+(1/2)Uxx + (u^2 + v^2)u=0 => Vt=-[(1/2)Uxx+(u^2 + v^2)u] corresponds to the real component
-Ut+(1/2)Vxx+(u^2 + v^2)v=0  => Ut=(1/2)Vxx+(u^2 +v^2)v corresponds to the imaginary component


L=10
0<t<10
0<x<10

n=100; 
dx=10/100
dt=10/100
'''


    

#Functions/Classes to use
L=1
N_x=500
dx=L/N_x
dt=0.1*(dx**2)
# x_max=1
# x_min=0

# t_max=1
# t_min=0
t=0
x=0

# iter_T=int((t_max-t_min)/dt)
# iter_X=int((x_max-x_min)/dx)
T=np.array([t+(i*dt) for i in range(N_x)])
X=np.array([x+(j*dx) for j in range(N_x)])

A= 1.0 # Amplitude 
sigma=5.0 #Wavelength 
k=2.0 #Wavenumber

def f(x): 
    return A*np.exp(-(x-0.5)**2)*np.cos(k*x) #The function correspoding to the real component. 
def g(x): 
    return A*np.exp(-(x-0.5)**2)*np.sin(k*x) #The function corresponding to the imaginary component.

def D2_Central_Diff(Grid,i,j): 
    if(j==0): 
        return (Grid[i][j+1]-Grid[i][j])/dx 
    elif(j==len(Grid[0])-1):
        return (Grid[i][j]-Grid[i][j-1])/dx
    else: 
        return (Grid[i][j+1]-2*Grid[i][j]+Grid[i][j-1])/dx**2

def Ut(Grid1,Grid2,i,j,K): 
    non_linear_term=((Grid1[i][j]+K)**2)+((Grid2[i][j]+K)**2)
    return (1/2)*D2_Central_Diff(Grid2,i,j)+np.round(non_linear_term,7)*(Grid1[i][j]+K)  

def Vt(Grid1,Grid2,i,j,K): 
    non_linear_term=((Grid1[i][j]+K)**2)+((Grid2[i][j]+K)**2)
    return -(1/2)*D2_Central_Diff(Grid1,i,j)-np.round(non_linear_term,7)*(Grid2[i][j]+K) 

def PrintComplexNumber(u,v): 
    if(v<0): 
        print(f"{u} {v}i")
    else: 
        print(f" {u}+{v}i ")


class Schrodinger: 
    def __init__(self,X,T,L):
        self.T=T
        self.X=X
        self.L=L
        
           
    def RK_Method(self): 
       k1=0
       k2=0
       k3=0
       k4=0

       j1=0
       j2=0
       j3=0
       j4=0

       T=self.T
       X=self.X
       
       New_U_Grid=np.zeros((len(T),len(X)))
       New_V_Grid=np.zeros((len(T),len(X)))
       
       for i in range(len(T)-1):
           for j in range(len(X)-1):
               #Boundary Conditions
               if(j==0 or j==len(X)-1): 
                  New_U_Grid[i][j]=0
                  New_V_Grid[i][j]=0
                  continue
               #Initial Conditions
               elif(i == 0): 
                   New_U_Grid[i][j]=f(X[j])
                   New_V_Grid[i][j]=g(X[j])
                   continue
               
       U_Grid=New_U_Grid
       V_Grid=New_V_Grid

       for i in range(len(T)-1):
           for j in range(len(X)-1):
               
               k1=dt*Ut(U_Grid,V_Grid,i,j,0)             
               j1=dt*Vt(U_Grid,V_Grid,i,j,0)

               k2=dt*Ut(U_Grid,V_Grid,i,j,k1/2)
               j2=dt*Vt(U_Grid,V_Grid,i,j,j1/2) 

               k3=dt*Ut(U_Grid,V_Grid,i,j,k2/2)
               j3=dt*Vt(U_Grid,V_Grid,i,j,j2/2)

               k4=dt*Ut(U_Grid,V_Grid,i,j,k3)
               j4=dt*Vt(U_Grid,V_Grid,i,j,j3)

            #    print(f"Iter: {i,j}")
            #    print(f"j1:{j1}, j2:{j2}, j3:{j3}, j4:{j4}")
            #    print(f"k1:{k1}, k2:{k2}, k3:{k3}, k4:{k4}")
            #    print()
            
               U_Grid[i+1][j]=U_Grid[i][j]+np.round((1/6)*(k1+(2*k2)+(2*k3)+k4),7)
               V_Grid[i+1][j]=V_Grid[i][j]+np.round((1/6)*(j1+(2*j2)+(2*j3)+j4),7)

       return U_Grid,V_Grid

def main(): 


    S=Schrodinger(X,T,L)
    U_Grid,V_Grid=S.RK_Method()

    Position=[]
    Lengths=[]
    for i in range(len(T)): 
        for j in range(len(X)):
            print()
            print(f"ψ({T[i]},{X[j]})=", end="\t")
            PrintComplexNumber(U_Grid[i][j], V_Grid[i][j])
            Position.append(j)
            Lengths.append(math.sqrt(U_Grid[i][j]**2 +V_Grid[i][j]**2))
    plt.plot(U_Grid,V_Grid)
    plt.title("Complex Plane ψ(x,t)")
    plt.xlabel("Real")
    plt.ylabel("Imaginary")
    plt.show()

main()

    
  
