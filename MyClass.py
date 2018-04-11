import math as np
import matplotlib.pyplot as plt
import numpy

def inputFunc(x): # исходная функция

    return (x / (10 * np.pi * np.sin(x)))

def csc(t,i): # косеканс
    return (1/np.sin(t))**i

def cot(t,i): # котангенс
    return (1 / np.tan(t)) ** i

def plotDx_4(x,b,xmin,xmax):
    xar= []
    yar = []
    y1 = [0,35]
    x1 = [xmin,xmin]
    x2 = [xmax,xmax]
    while(x<=b):
        yar.append((1/(10*np.pi)) *( 4*(cot(x,3)*(-csc(x,1))-5*cot(x,1)*csc(x,3))+x*(5*csc(x,5)+cot(x,4)*csc(x,1)+
                            18*cot(x,2)*csc(x,3))))
        xar.append(x)
        x+=0.01
    ma = 24*MaxDx4Func(xmax)
    plt.figure("График 4 производной")
    plt.grid(True)
    leg1,leg2,leg3 = plt.plot(xar,yar,x1,y1 , x2 ,y1)
    plt.legend((leg1, leg2, leg3), ("Исходный график", "a", "b"))
    plt.plot(xmax, ma ,'o')

    print("Max 4 Dx = " , ma )



def MaxDx4Func(x): # функ-ция считает макс значение 4 производной

    return ((1/(10*np.pi)) *( 4*(cot(x,3)*(-csc(x,1))-5*cot(x,1)*csc(x,3))+x*(5*csc(x,5)+cot(x,4)*csc(x,1)+
                            18*cot(x,2)*csc(x,3))))/24

def step (a,b): # находим шаг (параметры - какой нибудь начальный отрезок)
    def maxx(x0, h, b):# считаем максимальны x_ для |w_n+1| в лоб
        max = 0
        x_ = x0
        ans=0
        while (x_<= b):
            tmp = abs(x_ * (x_ - x0) * (x_ - x0 - h) * (x_ - x0 - 2 * h) * (x_ - x0 - 3 * h))
            if (max < tmp):
                max = tmp
                ans = x_
            x_ += 0.001
        return ans

    x3 = b # правая граница
    x0 = a  # левая граница
    h = (x3 - x0) / 3  # шаг
    R3 = 1 # погрешность
    x_ = maxx(x0,h,b) # max x_ для |w_n+1|
    max =MaxDx4Func(b) # само max значени |w_n+1|

    while (abs(R3/inputFunc(x_)) > 10 ** -3):
        R3 =  (max* x_ * (x_ - x0) * (x_ - x0 - h) * (x_ - x0 - 2 * h) * (x_ - x0 - 3 * h))
        x0 = x0 + 0.01
        h = (x3 - x0) / 3
        x_ = maxx(x0,h,b)
        print("x0 = ",x0,"R3 = ",R3/inputFunc(x_),"h  =", h)

    return x0,h

def tabalatedFunc(a, b,h,xarr,Tatulated): #Табуляция, параметры (суженые границы [a,b], шаг,массив  x и у для сохранения )
    x=a

    while x <= b:
        Tatulated.append(inputFunc(x))
        xarr.append(x)
        x += h
        x = numpy.round(x,4)

def lagranz(x,y,t): # хз, не помню что тут
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

def interpolPolynomNewton(X,x,y,h):# тоже не помню - ты диктовала

    dy0 = y[1]-y[0]
    dy1 = y[2]-y[1]
    dy2 = y[3]-y[2]

    dy_2 = dy1 - dy0

    dy1_2 = dy2 - dy1

    dy_3 = dy1_2-dy_2

    return y[0] + dy0*((X-x[0])/h) + dy_2*((X-x[0])*(X-x[1]))/(2*h**2)+dy_3*((X-x[0])*(X-x[1])*(X-x[2]))/(6*h**3)

def linearSpline(x,y,h,delt): #параметры - (массивы х,у ,шаг исходный, шаг для оценки погрешности(гораздо меньше чем исходный))
    A = [[1, 0, 0, 0, 0, 0],
         [0, 0, 1, 0, 0, 0],
         [0, 0, 0, 0, 1, 0],
         [1, h, 0, 0, 0, 0],
         [0, 0, 1, h, 0, 0],
         [0, 0, 0, 0, 1, h]]

    b = [[y[0]], [y[1]], [y[2]], [y[1]], [y[2]], [y[3]]]

    ans = numpy.linalg.solve(A, b)
    X = x[0]
    ansY = []
    ansX = []

    while X < x[1]:
        ansY.append(ans[0] + ans[1] * (X - x[0]))
        ansX.append(X)
        X = X + delt

    while X < x[2]:
        ansY.append(ans[2] + ans[3] * (X - x[1]))
        ansX.append(X)
        X = X + delt

    while X <= x[3]:
        ansY.append(ans[4] + ans[5] * (X - x[2]))
        ansX.append(X)
        X = X + delt

    return ansX,ansY

def parabolSpline(x,y,h,delt): # аналогично линейному
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
        X=X+delt


    while X < x[2]:
        ansY.append(ans[3]+ans[4]*(X-x[1])+ans[5]*(X-x[1])**2)
        ansX.append(X)
        X=X+delt


    while X <=x[3]:
        ansY.append(ans[6] + ans[7] * (X - x[2]) + ans[8] * (X - x[2]) ** 2)
        ansX.append(X)
        X = X + delt

    return ansX,ansY


def cubSpline(x,y,h,delt):# аналогично линейному
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
        X = X + delt

    while X < x[2]:
        ansY.append(ans[4] + ans[5] * (X - x[1]) + ans[6] * (X - x[1]) ** 2 + ans[7] * (X - x[1])**3)
        ansX.append(X)
        X = X + delt

    while X <= x[3]:
        ansY.append(ans[8] + ans[9] * (X - x[2]) + ans[10] * (X - x[2]) ** 2 + ans[11] * (X - x[2]) ** 3)
        ansX.append(X)
        X = X + delt

    return ansX,ansY

def errorLagranz(x,y): # погрешность лагранже , параметры  - массив х,у
    err = []
    xerr=[]
    X = x[0]
    while (X < x[3]):
        err.append(abs(inputFunc(X) - lagranz(x, y, X)))
        xerr.append(X)
        X = X + 0.008
    return xerr, err

#т.к лагранже и ньютон совпадают будет плюсовать немного разный шаг - чтобы точки не сливались на графике (0.008 и 0.006)

def errorNyton(x,y,h):
    err = []
    xerr=[]

    X = x[0]
    while (X <= x[3]):
        err.append(abs(inputFunc(X)- interpolPolynomNewton(X, x, y,h)))
        xerr.append(X)
        X = X + 0.006

    return xerr, err

def OutPutTab(M,xarr,Tatulated): # вывод таблицы табуляции (параметры - размер, массив х , массив у)
    print("       Табуляция функции")
    print("-------------------------------------------")
    for i in range(M):
        print('|  x{0} = {1}   |   y{0} ={2}  |'.format(i,numpy.round(xarr[i], 5),
                                      numpy.round(Tatulated[i], 16)))
    print("-------------------------------------------")

if __name__ == "__main__":

 # пусть начальный отрезок такой
     a=0.1
     b= 2.5
# будем сужать пока не будет заданная точность
     temp = step(a, b)

     a=temp[0]
     h=temp[1]

     x = [] # табул. массив х с шагом h
     y = []  # табул.массив y с шагом h


     plotDx_4(0.1,2.55,a,b)

     x2 = [] # табул. массив х с шагом для оценки погрешности
     y2 = [] # табул. массив у с шагом для оценки погрешности

     tabalatedFunc(a,b,h,x,y)

     tabalatedFunc(a, b, 0.001, x2, y2)

     OutPutTab(len(x),x,y)

#Интерполяционный многочлен Лагранже

     xnew = numpy.linspace(numpy.min(x), numpy.max(x), 4)
     ynew = [lagranz(x, y, i) for i in xnew]
     #print("Интерполяционный многочлен Лагранже ", numpy.round(ynew, 16))

# Интерполяционный многочлен Ньютона

     ynw = [interpolPolynomNewton(X, x, y,h) for X in x]
     #print("Интерполяционный многочлен Ньютона ", numpy.round(ynw, 16))

     plt.figure("Многочлены")
     leg1,leg2,leg3 = plt.plot(x2,y2,x,ynew,x,ynw)
     plt.grid(True)
     plt.legend((leg1, leg2,leg3),("Исходный график","Многочлен Лагранже", "Многочлен Ньютона"))

     LinSp= linearSpline(x,y,h,h)
     #print("Линейный сплайн ", numpy.round(LinSp[1], 16))

     plt.figure("Сплайны")
     X_Y_parabol_sp = parabolSpline(x,y,h,0.001)
     X_Y_cub_sp = cubSpline(x,y,h,0.001)
     X_Y_lin_sp = linearSpline(x,y,h,0.001)

     plt.grid(True)



     leg1, leg2, leg3, leg4 = plt.plot(x2,y2 ,'r', X_Y_cub_sp[0], X_Y_cub_sp[1],X_Y_parabol_sp[0],X_Y_parabol_sp[1], LinSp[0],LinSp[1], 'g-')
     plt.legend((leg1, leg2, leg3, leg4),
                ("Исходный график", "Кубический сплайн", "Параболический сплайн", "Линейный сплайн"))

     errCubX = []  # массив погреш х для куб , порабод сплайна
     errCubY = []  # массив погреш у для куб сплайна
     errPorabY = []  # массив погреш у для парабол сплайна
     errLin=[]   # массив погреш у для лин сплайна
     errLinX=[]  # массив погреш х для лин  сплайна

# погрешности лин, парабол и куб сплайна, лень было в функцию выносить  - поэтому тут будет

     for i in range(0, len(x2), 4): # тут 421 точка , чтобы не слипалось возьмем через 4
         errCubY.append(abs(float(y2[i]) - float(X_Y_cub_sp[1][i])))
         errPorabY.append(abs(float(y2[i]) - float(X_Y_parabol_sp[1][i])))
         errCubX.append(x2[i])

     for i in range(0, len(x2), 3): # чтобы линейный сплайн не слипался с парабол - возьмем лин через 3
          errLin.append(abs(float(y2[i]) - float(X_Y_lin_sp[1][i])))
          errLinX.append(x2[i])

# погрешности лагранже и ньютна
     errLagr = errorLagranz(x,y)
     errNyt = errorNyton(x,y,h)

     plt.figure("Погрешность 1 ")
     plt.grid(True)
     leg1,leg2 = plt.plot(errLagr[0],errLagr[1],'o',errNyt[0],errNyt[1],'o')
     plt.legend((leg1,leg2), ("Погрешность мет.Лагранжа","Погрешность мат.Ньютона"))

     plt.figure("Погрешность 2 ")
     plt.grid(True)
     leg1,leg2,leg3 = plt.plot(errCubX,errCubY,'o',errCubX,errPorabY,'o',errLinX,errLin,'o')
     plt.legend((leg1, leg2,leg3),("Погрешность куб. сплайн", "Погрешность парабол. сплайн","Погрешность лин.сплайн"))
     plt.show()








