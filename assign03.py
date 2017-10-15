# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 13:18:16 2017

@author: zq806676
"""
import numpy as np
import matplotlib.pyplot as plt

def mark(ymin,ymax,N,pa=1e5,pb=200,f=1e-4,dens=1,l=2.4e6):
    p = np.zeros(N + 1)
    u = np.zeros(N + 1)
    e = np.zeros(N + 1)
    y = np.zeros(N + 1)
    grad = np.zeros(N + 1)
    exac = np.zeros(N + 1)
    dy = (ymax - ymin) / N

    for i in range(0,N + 1):
        y[i] = ymin + dy * i
        p[i] = pa + pb * np.cos(y[i] * np.pi / l) #calculate the pressure on each point#

##Calculate the begining point and ending point##
    grad[N] = (p[N] - p[N - 1]) / dy
    grad[0] = (p[1] - p[0]) / dy
##Calculate the body part##
    for i in range(1, N):
        grad[i] = (p[i + 1] - p[i - 1]) / (dy * 2)
        
    for i in range(0,N + 1):
        u[i] = -grad[i] / (dens * f) #calculate the numerical winds#
        exac[i] = pb * np.pi * np.sin(np.pi * y[i] / l) / (dens * f * l) #the exact wind#
        e[i] = abs(exac[i] - u[i]) #error#

    return u,exac,e,y/dy

U1,EXAC1,E1,y = mark(0.0,1e6,10)

print("the numerically evaluated winds are",U1,
      "the exact winds are",EXAC1,
      "the erros are",E1)

plt.figure(1)
plt.plot(y,U1,label='numerical')
plt.plot(y,EXAC1,label='exact')
plt.ylabel("winds")
plt.xlabel("points")
plt.legend()
plt.show()

plt.figure(2)
plt.plot(y,E1,label="error")
plt.ylabel("error")
plt.xlabel("points")
plt.legend()
plt.show()
