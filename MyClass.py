import math as np
import matplotlib.pyplot as plt
import numpy
import scipy as sp

Tatulated = []
xarr = []

def inputFunc(x): # исходная функция
    #return 5*x*np.pow((5 * np.pi +2*x),-1/4) # Таня
    return (x / (10 * np.pi * np.sin(x)))

def MaxDx4Func(x): # функ-ция считает макс значение 4 производной
    def csc(t,i):
        return (1/np.sin(t))**i

    def cot(t,i):
        return (1 / np.tan(t)) ** i

    return (1/(10*np.pi)) *( 4*(cot(x,3)*(-csc(x,1))-5*cot(x,1)*csc(x,3))+x*(5*csc(x,5)+cot(x,4)*csc(x,1)+
                              18*cot(x,2)*csc(x,3)))

def tabalatedFunc(a, b,h):  # границы [a,b], шаг
    x=a
    if (x==0):
        x += h
        Tatulated.append(inputFunc(x))
        xarr.append(x)
        x += h
    while x <= b:
        Tatulated.append(inputFunc(x))
        xarr.append(x)
        x += h


def plot():
    plt.figure("График какой то функции ")
    plt.plot(xarr, Tatulated)
    plt.grid()
    plt.show()

def  alpha(N):
    a =[[1,xarr[0],xarr[0]**2,xarr[0]**3],[1,xarr[1],xarr[1]**2,xarr[1]**3],
        [1,xarr[2],xarr[2]**2,xarr[2]**3],[1,xarr[3],xarr[3]**2,xarr[3]**3]]
    a = numpy.zeros(N,N)
    for i in range(N):
        a[i][0] = 1
    for i in range(N):
        for j in range(1,N):
            a[i][j] =



    b = numpy.zeros((N,1))
    for i in range(4):
        b[i][0] = Tatulated[i]
    return numpy.linalg.solve(a,b)

def interpolPolynom(x):
    a = alpha()
    return a[0]+a[1]*x+a[2]*x**2+a[3]*x**3

def f1(i):
    return (Tatulated[i+1] - Tatulated[i])/(xarr[i+1]-xarr[i])
def f2(i):
    return  (f1(i+1) - f1(i))/(xarr[i+2]-xarr[i])
def f3(i):
    return ((f2(i+1) - f2(i)) / xarr[i+3] - xarr[i])

def interpolPolynomNewton(x):
    return Tatulated[0] + f1(0)*(x-xarr[0])+f2(0)*(x-xarr[0])*(x-xarr[1])+\
           f3(0)*(x-xarr[0])*(x-xarr[1])*(x-xarr[2])

def tablFuncZad2():
    tablZad2=[]
    for i in range(4):
        tablZad2.append(interpolPolynomNewton(i))
    return tablZad2

def linearSpline():
    y1 = tablFuncZad2()


    return
def OutPutTab(M):
    print("       Табуляция функции")
    print("--------------------------------")
    for i in range(M):
        print('|  x{0} = {1}   |   y{0} ={2}  |'.format(i,numpy.round(xarr[i], 5),
                                      numpy.round(Tatulated[i], 5)))

if __name__ == "__main__":
    tabalatedFunc(0,3,0.157)
    OutPutTab(19)
    MaxD4 = MaxDx4Func(3)



    #
    # tmp = 2.5
    # ans = interpolPolynom(1)
    # print("Интерполяционный многочлен Лагранже в ", tmp, " =", float(interpolPolynom(tmp)))
    #
    #
    # tmp = 2.5
    # print("Интерполяционный многочлен Ньютона в  ", tmp," =", interpolPolynomNewton(tmp))



