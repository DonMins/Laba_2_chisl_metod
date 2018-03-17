import math as np
import matplotlib.pyplot as plt

Tatulated  = []
xarr =[]

def inputFunc(x):
    return (x / (10 * np.pi * np.sin(x)))

def tabalatedFunc(a,b): # границы [a,b]
 x = 0.1
 Tatulated.append(x/(inputFunc(x)))
 xarr.append(x)
 x=1
 while x<=b:
    Tatulated.append(x/(10*np.pi *np.sin(x)))  # 7 вариант
    xarr.append(x)
    x+=1
    #Tatulated.append(5*x*np.pow((5 * np.pi +2*x),-1/4))  # 12 вариант

def plot ():
    plt.figure("График какой то функции ")
    plt.plot(xarr,Tatulated)
    plt.grid()
    plt.show()

if __name__ == "__main__":
    atom = ;


