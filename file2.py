import numpy as np
import matplotlib.pyplot as plt

##using the numerical methods in question1##
def mark01(ymin,ymax,N,pa=1e5,pb=200,f=1e-4,dens=1,l=2.4e6):
    p = np.zeros(N + 1)
    u = np.zeros(N + 1)
    e = np.zeros(N + 1)
    y = np.zeros(N + 1)
    grad = np.zeros(N + 1)
    exac = np.zeros(N + 1)
    dy = (ymax - ymin) / N

    for i in range(0,N + 1):
        y[i] = ymin + dy * i
        p[i] = pa + pb * np.cos(y[i] * np.pi / l)

    grad[N] = (p[N] - p[N - 1]) / dy
    grad[0] = (p[1] - p[0])/dy
    for i in range(1,N):
        grad[i] = (p[i + 1] - p[i - 1]) / (dy * 2)

    for i in range(0,N + 1):
        u[i] = -grad[i] / (dens * f)
        exac[i] = pb * np.pi * np.sin(np.pi * y[i] / l)/(dens * f * l)
        e[i] = abs(exac[i] - u[i])

    return u,exac,e,y/dy

##modify the methods used on beginning point and ending point, change another method to calculate the numerical solutions of winds on main parts##
def mark02(ymin,ymax,N,pa=1e5,pb=200,f=1e-4,dens=1,l=2.4e6):
    p = np.zeros(N + 1)
    u = np.zeros(N + 1)
    e = np.zeros(N + 1)
    y = np.zeros(N + 1)
    grad = np.zeros(N + 1)
    exac = np.zeros(N + 1)
    dy = (ymax - ymin) / N

    for i in range(0,N + 1):
        y[i] = ymin + dy * i
        p[i] = pa + pb * np.cos(y[i] * np.pi / l)

    grad[N] = (3 * p[N] - 4 * p[N - 1] + p[N - 2]) / (2 * dy)
    grad[N-1] = (p[N] - p[N - 2]) / (2 * dy)
    grad[0] = (-3 * p[0] + 4 * p[1] - p[2]) / (2 * dy)
    grad[1] = (p[2] - p[0]) / (2 * dy)
    for i in range(2,N - 1):
        grad[i] = ((p[i - 2] - 27 * p[i - 1] + 27 * p[i] - p[i + 1]) + (p[i - 1] - 27 * p[i] + 27 * p[i + 1] - p[i + 2])) / (2 * 24 * dy)

    for i in range(0,N + 1):
        u[i] = -grad[i] / (dens * f)
        exac[i] = pb * np.pi * np.sin(np.pi * y[i] / l)/(dens * f * l)
        e[i] = abs(exac[i] - u[i])

    return u,exac,e,y/dy
U1,EXAC1,E1,y = mark01(0.0,1e6,10)
U2,EXAC2,E2,y= mark02(0.0,1e6,10)
print ("the numerically evaluated winds are",U2,
      "the exact winds are",EXAC2,
      "the erros are",E2)


##display the different solutions using one figure##
plt.figure(1)
plt.plot(y,U1,label="numerical01")
plt.plot(y,U2,label='numerical02')
plt.plot(y,EXAC2,label='exact')
plt.ylabel("winds")
plt.xlabel("points")
plt.legend()
plt.show()

##compare the difference of these two methods##
plt.figure(2)
plt.plot(y,E1,label="error1")
plt.plot(y,E2,label="error2")
plt.ylabel("error")
plt.xlabel("points")
plt.legend()
plt.show()

##display the middle part of lines to make difference visible##
plt.figure(3)
plt.plot(y[1:8],E1[1:8],label="error1")
plt.plot(y[1:8],E2[1:8],label="error2")
plt.ylabel("error")
plt.xlabel("points")
plt.legend()
plt.show()

##find the possible value of the order of accuracy when using different N in the scond method##
def mark03(ymin,ymax,N,pa=1e5,pb=200,f=1e-4,dens=1,l=2.4e6):
    p = np.zeros(N + 1)
    u = np.zeros(N + 1)
    e = np.zeros(N + 1)
    y = np.zeros(N + 1)
    grad = np.zeros(N + 1)
    exac = np.zeros(N + 1)
    dy = (ymax - ymin) / N

    for i in range(0,N + 1):
        y[i] = ymin + dy * i
        if y[i] == (ymin + ymax) / 2.0:
            j = i
        p[i] = pa + pb * np.cos(y[i] * np.pi / l)

    grad[N] = (3 * p[N] - 4 * p[N - 1] + p[N - 2]) / (2 * dy)
    grad[N-1] = (p[N] - p[N - 2]) / (2 * dy)
    grad[0] = (-3 * p[0] + 4 * p[1] - p[2])/ (2 * dy)
    grad[1] = (p[2] - p[0]) / (2 * dy)
    for i in range(2,N - 1):
        grad[i] = ((p[i - 2] - 27 * p[i - 1] + 27 * p[i] - p[i + 1]) + (p[i - 1] - 27 * p[i] + 27 * p[i + 1] - p[i + 2])) / (2 * 24 * dy)

    for i in range(0,N + 1):
        u[i] = -grad[i] / (dens * f)
        exac[i] = pb * np.pi * np.sin(np.pi * y[i] / l) / (dens * f * l)
        e[i] = abs(exac[i] - u[i])

    return e[j],dy

E = np.zeros(5)
dy = np.zeros(5)
n = np.zeros(4)
for i in range(0,5):
    E[i], dy[i] = mark03(0.0,1e6,10*(2**i))

for i in range(0,4):
    n[i] = (np.log(E[i]) - np.log(E[i+1])) / (np.log(dy[i]) - np.log(dy[i+1]))

print(n)
plt.figure(4)
plt.plot(dy,E)
plt.xlabel("dy")
plt.ylabel("error")
plt.grid()
plt.show()