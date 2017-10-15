import numpy as np
import matplotlib.pyplot as plt


def mark(ymin,ymax,N,pa=1e5,pb=200,f=1e-4,dens=1,l=2.4e6):
    p = np.zeros(N+1)
    u = np.zeros(N+1)
    e = np.zeros(N+1)
    y = np.zeros(N+1)
    grad = np.zeros(N+1)
    exac = np.zeros(N+1)
    dy = (ymax - ymin)/N

    for i in range(0,N+1):
        y[i] = ymin + dy*i
        if y[i] == (ymin + ymax)/2.0:
            j = i
        p[i] = pa + pb*np.cos(y[i]*np.pi/l)

    grad[N] = (p[N] - p[N-1])/dy
    grad[0] = (p[1] - p[0])/dy
    for i in range(1,N):
        grad[i] = (p[i+1] - p[i-1])/(dy*2)

    for i in range(0,N+1):
        u[i] = -grad[i] / (dens * f)
        exac[i] = pb*np.pi*np.sin(np.pi*y[i]/l)/(dens*f*l)
        e[i] = abs(exac[i] - u[i])

    return e[j],dy

E = np.zeros(5)
dy = np.zeros(5)
n = np.zeros(4)
for i in range(0,5):
    E[i], dy[i] = mark(0.0,1e6,10*(2**i))

for i in range(0,4):
    n[i] = (np.log(E[i]) - np.log(E[i+1]))/(np.log(dy[i]) - np.log(dy[i+1]))

print(n)
plt.figure(1)
plt.plot(dy,E)
plt.xlabel("dy")
plt.ylabel("error")
plt.grid()
plt.show()