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

    return


def lagranz(x,y,t):
    z = 0

    for j in range(len(y)):
        p1 = 1
        p2 = 1
        for i in range(len(x)):
            if i == j:
                p1 = p1 * 1
                p2 = p2 * 1
            else:
                p1 = p1 * (t - x[i])
                p2 = p2 * (x[j] - x[i])
        z = z + y[j] * p1 / p2


    return z

def f1(i):
    return (Tatulated[i+1] - Tatulated[i])/(xarr[i+1]-xarr[i])
def f2(i):
    return  (f1(i+1) - f1(i))/(xarr[i+2]-xarr[i])
def f3(i):
    return ((f2(i+1) - f2(i)) / xarr[i+3] - xarr[i])

def interpolPolynomNewton(X,x,y,h):

    dy0 = y[1]-y[0]
    dy1 = y[2]-y[1]
    dy2 = y[3]-y[2]

    dy_2 = dy1 - dy0

    dy1_2 = dy2 - dy1

    dy_3 = dy1_2-dy_2

    return y[0] + dy0*((X-x[0])/h) + dy_2*((X-x[0])*(X-x[1]))/(2*h**2)+dy_3*((X-x[0])*(X-x[1])*(X-x[2]))/(6*h**3)

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
    x=xarr[10:14]
    y=Tatulated[10:14]

    xnew =numpy.linspace(numpy.min(x), numpy.max(x), 4)
    ynew = [lagranz(x, y, i) for i in xnew]
    print("Интерполяционный многочлен Лагранже ", numpy.round(ynew, 5))
    plt.plot(x, y, x,y,'o', xnew, ynew)
    plt.grid(True)

    ynw = [interpolPolynomNewton(X,x,y,0.157) for X in x]
    print("Интерполяционный многочлен Ньютона ", numpy.round(ynw, 5))

    plt.plot(x, y, x, y, 'o', xnew, ynw)
    plt.grid(True)
    plt.show()

    #


    #
    #
    # tmp = 2.5
    # print("Интерполяционный многочлен Ньютона в  ", tmp," =", interpolPolynomNewton(tmp))



