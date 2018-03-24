import math as np
import matplotlib.pyplot as plt
import numpy

Tatulated = []
xarr = []

def inputFunc(x): # исходная функция
    #return 5*x*np.pow((5 * np.pi +2*x),-1/4) # Таня
    return (x / (10 * np.pi * np.sin(x)))

def csc(t,i):
    return (1/np.sin(t))**i

def cot(t,i):
    return (1 / np.tan(t)) ** i

def Dx (x):
    return (csc(x,1)/(10*np.pi) - (x*cot(x,1)*csc(x,1))/10*(np.pi) )

def MaxDx4Func(x): # функ-ция считает макс значение 4 производной
    return (1/(10*np.pi)) *( 4*(cot(x,3)*(-csc(x,1))-5*cot(x,1)*csc(x,3))+x*(5*csc(x,5)+cot(x,4)*csc(x,1)+
                              18*cot(x,2)*csc(x,3)))

def step ():
    x3 = 2.5
    x0 = 1.5
    h = (x3 - x0) / 3
    R3 = 1
    x_ = (x3 + x0) / 2
    while (abs(R3) > 10 ** -4):
        R3 = 0.919 * x_ * (x_ - x0) * (x_ - x0 - h) * (x_ - x0 - 2 * h) * (x_ - x0 - 3 * h)
        x0 = x0 + 0.01
        h = (x3 - x0) / 3
        x_ = (x3 + x0) / 2
        print('x0 = {0} h = {1} R3={2}  x_ = {3}'.format(x0, h, R3, x_))

def spline3(a,b):
    x=a
    h=b-a
    S=[]
    X1=[]
    while x<b:
        S.append ((((x-b)**2 * (2*(x-a)+h))*inputFunc(a))/h**3 +(((x-a)**2 * (2*(b-x)+h))*inputFunc(b))/h**3 +
        ((x-b)**2 * (x-a) * Dx(a))/h**2 + ((x-a)**2* (x - b)* Dx(b))/h**2 )
        X1.append(x)
        x+=0.001

    return X1,S


def tabalatedFunc(a, b,h):  # границы [a,b], шаг
    x=a
    while x <= b:
        Tatulated.append(inputFunc(x))
        xarr.append(x)
        x += h
        x = numpy.round(x,4)

def plot():
    plt.figure("График какой то функции ")
    plt.plot(xarr, Tatulated)
    plt.grid()
    plt.show()

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

def interpolPolynomNewton(X,x,y,h):

    dy0 = y[1]-y[0]
    dy1 = y[2]-y[1]
    dy2 = y[3]-y[2]

    dy_2 = dy1 - dy0

    dy1_2 = dy2 - dy1

    dy_3 = dy1_2-dy_2

    return y[0] + dy0*((X-x[0])/h) + dy_2*((X-x[0])*(X-x[1]))/(2*h**2)+dy_3*((X-x[0])*(X-x[1])*(X-x[2]))/(6*h**3)

def linearSpline(x,y):
    lineSpline=[]
    for i in range(3):
        k = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
        b = y[i]
        x1 = [x[i], x[i+1]]
        y1 = [b + k * (x1[0] - x[i]), b + k * (x1[1] - x[i])]
        if(i==0):
            lineSpline.append(y1[0])
            lineSpline.append(y1[1])
            plt.plot(x1, y1, 'r')
            plt.grid(True)
        if(i==1):
            plt.plot(x1, y1, 'r')
        if(i==2):
            lineSpline.append(y1[0])
            lineSpline.append(y1[1])
            plt.plot(x1, y1, 'r')


    plt.show()
    return lineSpline

def parabolSpline(x,a,b):
    A = [[1,(x[1]-x[0]),(x[1]-x[0])**2,0,0,0,0,0,0],
         [0,0,0,1,(x[2]-x[1]),(x[2]-x[1])**2,0,0,0],
         [0,0,0,0,0,0,1,(x[3]-x[2]),(x[3]-x[2])**2],
         [0,1,2*(x[1]-x[0]),0,-1,0,0,0,0],
         [0,0,0,0,1,2*(x[2]-x[1]),0,-1,0],
         [0,1,2*(a-x[0]),0,0,0,0,0,0],
         [0,0,0,0,0,0,0,1,2*(b-x[2])],
         [0,0,1,0,0,-1,0,0,0],
         [0,0,0,0,0,1,0,0,-1]]

    B=[inputFunc(x[1]),inputFunc(x[2]),inputFunc(x[3]),0,0,0,0,0,0]

    ans= numpy.linalg.solve(A,B)
    X=x[0]
    P21=[]
    P21x = []
    while X<x[1]:
        P21.append(ans[0]+ans[1]*(X-x[0])+ans[2]*(X-x[0])**2)
        P21x.append(X)
        X=X+0.01
    plt.plot(P21x, P21)
    plt.grid(True)
    X = x[1]
    P22 = []
    P22x = []
    while X<x[2]:
        P22.append(ans[3]+ans[4]*(X-x[1])+ans[5]*(X-x[1])**2)
        P22x.append(X)
        X=X+0.01
    plt.plot(P22x, P22)
    X = x[2]
    P23 = []
    P23x = []
    while X < x[3]:
        P23.append(ans[6] + ans[7] * (X - x[2]) + ans[8] * (X - x[2]) ** 2)
        P23x.append(X)
        X = X + 0.01
    plt.plot(P23x, P23)

    plt.show()





    # i=0
    # a=(y[i]-y[i+1])/(2*x[i]-x[i]**2-x[i+1]**2)
    # b=-2*a*x[i]
    # c= y[i+1]-a*x[i+1]**2+2*a*x[i]
    # y1=a*x**2+b*x+c
    # A = [ [0,2,0,0 ], [x[1]-x[0],(x[1]-x[0])**2,0,0],  [ x[2]-x[0] ,(x[2]-x[0])**2,(x[2]-x[1])**2,0],[ x[3]-x[0] ,(x[3]-x[0])**2,(x[3]-x[1])**2,(x[3]-x[2])**2]]
    # b = numpy.zeros((4, 1))
    # b[0][0]=0
    # for i in range(1,4):
    #     b[i][0] = y[i] - y[0]
    # ans = numpy.linalg.solve(A, b)
    # f = float(y[0] + ans[0]*(X-x[0])+ans[1] * (X - x[0])**2 + ans[2] * (X - x[1])**2 + ans[3] * (X - x[2])**2)

def cubSpline(X,x,y):
    A = [ [x[1]-x[0],(x[1]-x[0])**3,0,0], [x[1]-x[0],(x[1]-x[0])**2,0,0],  [ x[2]-x[0] ,(x[2]-x[0])**2,(x[2]-x[1])**2,0],[ x[3]-x[0] ,(x[3]-x[0])**2,(x[3]-x[1])**2,(x[3]-x[2])**2]]
    b = numpy.zeros((4, 1))
    b[0][0]=0
    for i in range(1,4):
        b[i][0] = y[i] - y[0]
    ans = numpy.linalg.solve(A, b)
    f = float(y[0] + ans[0]*(X-x[0])+ans[1] * (X - x[0])**2 + ans[2] * (X - x[1])**2 + ans[3] * (X - x[2])**2)
    return f



def OutPutTab(M):
    print("       Табуляция функции")
    print("--------------------------------")
    for i in range(M):
        print('|  x{0} = {1}   |   y{0} ={2}  |'.format(i,numpy.round(xarr[i], 5),
                                      numpy.round(Tatulated[i], 16)))

if __name__ == "__main__":
     tabalatedFunc(2.2,2.5,0.10)
     OutPutTab(len(xarr))
     x = xarr
     y = Tatulated

     xnew = numpy.linspace(numpy.min(x), numpy.max(x), 4)
     ynew = [lagranz(x, y, i) for i in xnew]
     print("Интерполяционный многочлен Лагранже ", numpy.round(ynew, 16))
     plt.plot(xnew, ynew)
     plt.grid(True)
     plt.show()
     ynw = [interpolPolynomNewton(X, x, y, 0.1) for X in x]
     print("Интерполяционный многочлен Ньютона ", numpy.round(ynw, 16))
     print("Линейный сплайн ", numpy.round((linearSpline(x, y)), 16))
     ans = spline3(xarr[0],xarr[1])
     ans1 = spline3(xarr[1], xarr[2])
     ans2 = spline3(xarr[2], xarr[3])
     plt.figure("3 spline")
     plt.plot(x,y,ans[0],ans[1],ans1[0],ans1[1],ans2[0],ans2[1])
     plt.grid(True)
     plt.show()
     parabolSpline(xarr,1,2)
#

#     # yn3 = [parabolSpline(X, x, y) for X in x]
#     # print("Параболический сплайн ", numpy.round((yn3), 5))
#     # plt.plot(x, y, x, y, 'o', xnew, yn3,xnew,yn2)
#     # plt.grid(True)
#     #plt.show()




