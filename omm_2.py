# -*- coding: utf-8 -*-
"""
Created on Sat May 22 20:40:34 2021

@author: andre
"""

import numpy as np
import matplotlib.pyplot as plt

M = 50
X = 10

N = 50
Y = 5

J = 50
T = 1

u = np.zeros((M,N,2*J + 1))
x = np.linspace(0,X,M)
y = np.linspace(0,Y,N)
t = np.linspace(0,T,J)

h1 = X/(M-1)
h2 = Y/(N-1)
tau = T/(J-1)

g1 = tau/(h1**2)
g2 = tau/(h2**2)

def F1(m,n,j):
    return 2*g2*u[m,n+1,j-1] + (1 - 4*g2)*u[m,n,j-1] + 2*g2*u[m,n-1,j-1]

def progonka_1(n,j):
    d = np.zeros(M)
    sig = np.zeros(M)
    d[1] = 0
    sig[1] = 0
    A = 2*g1
    B = 1 + 4*g1
    C = 2*g1
    for i in range(1,M-1):
        Fi = -F1(i,n,j)
        d[i+1] = C/(B - A*d[i])
        sig[i+1] = (Fi - A*sig[i])/(A*d[i] - B)
    u[M-1,n,j]=0
    for i in range(M-1,0,-1):
        u[i-1,n,j] = d[i]*u[i,n,j] + sig[i]
        
def F2(m,n,j):
    return 2*g1*u[m+1,n,j-1] + (1 - 4*g1)*u[m,n,j-1] + 2*g1*u[m-1,n,j-1]

def progonka_2(m,j):
    d = np.zeros(N)
    sig = np.zeros(N)
    d[1] = 1
    sig[1] = 0
    A = 2*g2
    B = 1+4*g2
    C = 2*g2
    for i in range(1,N-1):
        Fi = -F2(m,i,j)
        d[i+1] = C/(B - A*d[i])
        sig[i+1] = (Fi - A*sig[i])/(A*d[i] - B)
    u[m,N-1,j] = sig[N-1]/(1 - d[N-1])
    for i in range(N-1,0,-1):
        u[m,i-1,j] = d[i]*u[m,i,j] + sig[i]

# u[:,:,0] = np.sin(np.pi*x*h1) * np.cos(2*np.pi*y*h2)
# print(u)
for m in range(0,M-1):
    for n in range (0,N-1):
        u[m,n,0] = np.sin(np.pi*m*h1) * np.cos(2*np.pi*n*h2)

fig = plt.figure(figsize=(8,6))
plt.pcolormesh(x,y,u[:,:,0],shading='nearest')
plt.colorbar()


for j in range(1,2*J,2):
    for n in range(1,N-1):
        progonka_1(n,j)
    for m in range(1,M-1):
        progonka_2(m,j+1)


fig = plt.figure(figsize=(8,6))
plt.pcolormesh(x,y,u[:,:,0],shading='auto')
plt.colorbar()
        
fig = plt.figure(figsize=(8,6))
plt.pcolormesh(x,y,u[:,:,6],shading='auto')
plt.colorbar()
        
fig = plt.figure(figsize=(8,6))
plt.pcolormesh(x,y,u[:,:,-1],shading='auto')
plt.colorbar()

def U(x0,y0,t0):
    return np.exp(-80*np.pi*np.pi*t0)*np.sin(np.pi*x0) * np.cos(2*np.pi*y0)

e=0
U0 = np.zeros((M,N,2*J + 1))
for m in range(0,M-1):
    for n in range(0,N-1):
        for j in range(0,2*J,2):
            U0[m,n,j] = U(m*h1,n*h2,j*tau)
            e+=(u[m,n,j]-U(m*h1,n*h2,j*tau))**2

fig = plt.figure(figsize=(8,6))
plt.pcolormesh(x,y,U0[:,:,5],shading='auto')
plt.colorbar()














