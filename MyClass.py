import math as np
import matplotlib.pyplot as plt
Tatulated  = []
a=0

b=3
#=a+0.01 # 7 вариант
x=a      # 12 вариант
n=0
while x<=b:
    #Tatulated.append(x/(10*np.pi *np.sin(x)))  # 7 вариант
    Tatulated.append(5*x*np.pow((5 * np.pi +2*x),-1/4))  # 12 вариант
    x+=0.01
    n+=1
    
#xarr = [i for i in range(300)] # 7 вариант
xarr = [i for i in range(301)] # 12 вариант
plt.figure("График какой то функции")
plt.plot(xarr,Tatulated)
plt.grid()
plt.show()







