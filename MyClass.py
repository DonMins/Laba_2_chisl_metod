import math as np
import matplotlib.pyplot as plt
import numpy


def inputFunc(x): # исходная функция
    return 5*x*np.pow((5 * np.pi +2*x),-1/4) # Таня
    #return (x / (10 * np.pi * np.sin(x)))

def csc(t,i):
    return (1/np.sin(t))**i

def cot(t,i):
    return (1 / np.tan(t)) ** i

def Dx (x):
    return (5*(3*x+10*np.pi))/(2*np.pow((5*np.pi+2*x),5/4))
    #return (csc(x,1)/(10*np.pi) - (x*cot(x,1)*csc(x,1))/10*(np.pi) )



def MaxDx4Func(x): # функ-ция считает макс значение 4 производной
    return (2935*x)/(16*np.pow((2*x +5*np.pi),17/4))-(225*x)/(2*np.pow((2*x +5*np.pi),13/4))
    #return ((1/(10*np.pi)) *( 4*(cot(x,3)*(-csc(x,1))-5*cot(x,1)*csc(x,3))+x*(5*csc(x,5)+cot(x,4)*csc(x,1)+
     #                         18*cot(x,2)*csc(x,3))))/24

def step (a,b):
    x3 = b
    x0 = a
    h = (x3 - x0) / 3
    R3 = 1
    x_ = (x3 + x0) / 2
    while (abs(R3) > 10 ** -3):
        R3 = MaxDx4Func(b) * x_ * (x_ - x0) * (x_ - x0 - h) * (x_ - x0 - 2 * h) * (x_ - x0 - 3 * h)
        x0 = x0 + 0.01
        h = (x3 - x0) / 3
        x_ = (x3 + x0) / 2

    return x0,h



def tabalatedFunc(a, b,h,xarr,Tatulated):  # границы [a,b], шаг
    x=a
    while x <= b:
        Tatulated.append(inputFunc(x))
        xarr.append(x)
        x += h
        x = numpy.round(x,4)

def plot(xarr,Tatulated):
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
    ansX = []
    ansY = []
    for i in range(3):
        k = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
        b = y[i]
        x1 = [x[i], x[i+1]]
        y1 = [b + k * (x1[0] - x[i]), b + k * (x1[1] - x[i])]
        if(i==0):
            lineSpline.append(y1[0])
            lineSpline.append(y1[1])
            ansX.append(x1)
            ansY.append(y1)

        if(i==1):
            ansX.append(x1)
            ansY.append(y1)
        if(i==2):
            lineSpline.append(y1[0])
            lineSpline.append(y1[1])
            ansX.append(x1)
            ansY.append(y1)

    return lineSpline,ansX,ansY

def parabolSpline(x,y,h):
    A = [[1, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 1, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 1, 0, 0],
         [1, h, h ** 2, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 1, h, h ** 2, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 1, h, h ** 2],
         [0, 1, 2 * h, 0, -1, 0, 0, 0, 0],
         [0, 0, 0, 0, 1, 2 * h, 0, -1, 0],
         [0, 0, 2, 0, 0, 0, 0, 0, 0]]

    b = [[y[0]], [y[1]],[y[2]],[y[1]], [y[2]], [y[3]], [0], [0], [0]]

    ans = numpy.linalg.solve(A, b)
    X=x[0]
    ansY = []
    ansX = []

    while X < x[1]:
        ansY.append(ans[0]+ans[1]*(X-x[0])+ans[2]*(X-x[0])**2)
        ansX.append(X)
        X=X+0.001
    X = x[1]

    while X < x[2]:
        ansY.append(ans[3]+ans[4]*(X-x[1])+ans[5]*(X-x[1])**2)
        ansX.append(X)
        X=X+0.001

    X = x[2]
    while X < x[3]:
        ansY.append(ans[6] + ans[7] * (X - x[2]) + ans[8] * (X - x[2]) ** 2)
        ansX.append(X)
        X = X + 0.001

    return ansX,ansY


def cubSpline(x,y,h):
    A = [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
         [1, h, h ** 2, h ** 3, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 1, h, h ** 2, h ** 3, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 1, h, h ** 2, h ** 3],
         [0, 1, 2 * h, 3 * h ** 2, 0, -1, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 1, 2 * h, 3 * h ** 2, 0, -1, 0, 0],
         [0, 0, 1, 3 * h, 0, 0, -1, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 1, 3 * h, 0, 0, -1, 0],
         [0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3*h]]

    b = [[y[0]], [y[1]], [y[2]], [y[1]], [y[2]], [y[3]], [0], [0], [0],[0],[0],[0]]

    ans = numpy.linalg.solve(A, b)
    X = x[0]
    ansY = []
    ansX = []
    while X < x[1]:
        ansY.append(ans[0] + ans[1] * (X - x[0]) + ans[2] * (X - x[0]) ** 2 + ans[3]*(X-x[0])**3)
        ansX.append(X)
        X = X + 0.001
    X = x[1]

    while X < x[2]:
        ansY.append(ans[4] + ans[5] * (X - x[1]) + ans[6] * (X - x[1]) ** 2 + ans[7] * (X - x[1])**3)
        ansX.append(X)
        X = X + 0.001

    X = x[2]
    while X < x[3]:
        ansY.append(ans[8] + ans[9] * (X - x[2]) + ans[10] * (X - x[2]) ** 2 + ans[11] * (X - x[2]) ** 3)
        ansX.append(X)
        X = X + 0.001

    return ansX,ansY



def OutPutTab(M,xarr,Tatulated):
    print("       Табуляция функции")
    print("--------------------------------")
    for i in range(M):
        print('|  x{0} = {1}   |   y{0} ={2}  |'.format(i,numpy.round(xarr[i], 5),
                                      numpy.round(Tatulated[i], 16)))

if __name__ == "__main__":
     a=1.5
     b= 2.6
     temp = step(a, b)
     a=temp[0]
     h=temp[1]



     x =[]
     y=[]

     tabalatedFunc(a,b,h,x,y)

     OutPutTab(len(x),x,y)

     xnew = numpy.linspace(numpy.min(x), numpy.max(x), 4)
     ynew = [lagranz(x, y, i) for i in xnew]
     print("Интерполяционный многочлен Лагранже ", numpy.round(ynew, 16))

     ynw = [interpolPolynomNewton(X, x, y,h) for X in x]
     print("Интерполяционный многочлен Ньютона ", numpy.round(ynw, 16))
     LinSp= linearSpline(x, y)
     print("Линейный сплайн ", numpy.round(LinSp[0], 16))

     plt.figure("Спалайны")
     X_Y_parabol_sp = parabolSpline(x,y,h)
     X_Y_cub_sp = cubSpline(x,y,h)

     plt.grid(True)
     x2=[]
     y2=[]
     tabalatedFunc(a, b, 0.001,x2,y2)

     errCubX = []
     errCubY = []
     errPorabY = []

     leg1, leg2, leg3, leg4,leg5 = plt.plot(x2,y2 ,'r', X_Y_cub_sp[0], X_Y_cub_sp[1],X_Y_parabol_sp[0],X_Y_parabol_sp[1], LinSp[1],LinSp[2], 'g-')
     plt.legend((leg1, leg2, leg3, leg4),
                ("Исходный график", "Кубический сплайн", "Параболический сплайн", "Линейный сплайн"))

     plt.show()

     errCubX=[]
     errCubY=[]
     errPorabY=[]

     def xz():
         err=[]
         X=x[0]
         while(X<x[3]):

             X

     for i in range(0,len(x2),5):
         errCubY.append(abs(float(y2[i])-float(X_Y_cub_sp[1][i])))
         errPorabY.append(abs(float(y2[i])-float(X_Y_parabol_sp[1][i])))

         errCubX.append(x2[i])

     plt.figure("Погрешность")
     plt.grid(True)


     plt.plot(errCubX,errCubY,'o',errCubX,errPorabY,'o')

     plt.show()








