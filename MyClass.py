import math as np
import matplotlib.pyplot as plt
import numpy

Tatulated = []
xarr = []


def inputFunc(x):
    #return 5*x*np.pow((5 * np.pi +2*x),-1/4) # Таня
    return (x / (10 * np.pi * np.sin(x)))

def tabalatedFunc(a, b,h):  # границы [a,b]
    x = 0.00001
    Tatulated.append(inputFunc(x))
    xarr.append(x)
    x = 1
    while x <= b:
        Tatulated.append(inputFunc(x))  # 7 вариант
        xarr.append(x)
        x += h

def plot():
    plt.figure("График какой то функции ")
    plt.plot(xarr, Tatulated)
    plt.grid()
    plt.show()

def  alpha():
    a =[[1,xarr[0],xarr[0]**2,xarr[0]**3],[1,xarr[1],xarr[1]**2,xarr[1]**3],
        [1,xarr[2],xarr[2]**2,xarr[2]**3],[1,xarr[3],xarr[3]**2,xarr[3]**3]]
    b = numpy.zeros((4,1))
    for i in range(4):
        b[i][0] = Tatulated[i]
    return numpy.linalg.solve(a,b)

def interpolPolynom(x):
    a = alpha()
    return a[0]+a[1]*x+a[2]*x**2+a[3]*x**3

if __name__ == "__main__":
    tabalatedFunc(0,3,1)
    print(Tatulated)
    print(xarr)
    ans = interpolPolynom(1)
    print(ans)
