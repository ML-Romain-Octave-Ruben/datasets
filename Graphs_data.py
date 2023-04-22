import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import math
import os
os.environ["Path"]=";".join(os.environ["Path"].split(";")[12:])
import h5py
import scipy.signal
import scipy.io
import scipy.integrate
import sklearn
import sklearn.linear_model
from sklearn.preprocessing import PolynomialFeatures
import tensorflow_datasets as tfds
import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import tensorflow_probability as tfp
import random
print("import done")

BATCH = 2

FN=f"batch{BATCH}_mit.mat"
arrays = {}
mat = h5py.File(FN)

epsilon = 0.05


# Features: I,Qc,Qd,Qdlin,T,Tdlin,V,discharge_dQdV,t
def create_df(cell, features=["t", "I", "V", "T", "Qc"], K = 10):
    k = mat["batch"]["cycles"][cell][0]
    obj = mat[k]
    df = {i:[] for i in features}
    df["cycle"] = []
    for i in range(1, len(obj["I"])):
        for j in features:
            df[j].extend(list(obj[obj[j][i,0]][0])) 
        N = len(obj[obj["I"][i,0]][0])
        df["cycle"].extend([i]*N)
    df =  pd.DataFrame(df)
    t = min(np.array(df.loc[(df.I<-epsilon)& (df.t>20)]["t"]) )
    df = df.loc[df.t<t]
    
    df = df.loc[df.cycle>=K]
    try:
        M = max(max(df.cycle)-K, K+2)
    except:
        print(df)
        print("t:", t)
        print("cell:", cell)
        M = 0
    df = df.loc[df.cycle<=M]
    return df


df = create_df(1)
cycle_step = 1





def clear_outliers2(M, K = 20, start = 20, B=0.1, C = 2, d = 0.1, Nn = 3):
    start = max(start, K)
    n = len(L)
    for i in range(start, len(M)-2):
        e = M[i+1]
        m = max([abs(M[i-j]-e) for j in range(K)])
        if m>B and abs(M[i+2]-e)>max(C*m, d):
            M[i+2] = M[i+1]
       
            
def clear_outliers(L, K = 3, start = 5):
    n = len(L)
    for i in range(start, n):
        V = sum(L[i-K:i])/K
        if abs(V)>0 and abs(L[i]-V)/abs(V)>0.3:
            L[i] = V
    
            
            
def dx_(L, K=1):
    N = len(L)
    return np.array([L[K+i]-L[i] for i in range(len(L)-K)])
            
           
            
           
def least_squares(X: np.array, y: np.array, alpha: float=0.0) -> np.array:
    return np.matmul(np.linalg.inv(np.matmul(X.T, X) + alpha) , np.matmul(X.T, y)) 




def smothing(X, Y, alpha = 0.01, k = 8):
    return scipy.signal.savgol_filter(Y, 2, 1)

    

"""
for i in range(0, int(max(df.cycle)), cycle_step):
    plt.plot(df.loc[df.cycle==i]["t"],df.loc[df.cycle==i]["I"])
plt.xlim((0, 50))
plt.title("Cumulative voltage (V) as a function of time (min)")
plt.xlabel("Time (min)")
plt.ylabel("i = f(t)")

plt.show()





for j in range(min(df.cycle), min(df.cycle)+1000, 50):
    dx = dx_(np.array(df.loc[df.cycle==j]["I"]))
    dy = dx_(np.array(df.loc[df.cycle==j]["V"]))
    plt.plot(np.array(df.loc[df.cycle==j]["I"])[1:], dy/dx)
    plt.ylim((-1000, 1000))
    plt.xlabel("I")
    plt.ylabel("dV/dI")
plt.title("dV/dI as a function of I")
plt.ylim((-100, 100))
plt.xlim((-0.2, 0.2))
plt.show()


"""


"""
a = 0
b = len(mat["batch"]["cycles"])
b=2
step=1

for I in range(a, b, step):
    Index = []
    df = create_df(I)
    L = []
    for i in range(int(min(df.cycle)), int(max(df.cycle)), cycle_step):
        Int = scipy.integrate.trapezoid(df.loc[df.cycle==i]["I"], df.loc[df.cycle==i]["t"])
        if Int>30:
            L.append(Int)
            Index.append(i)
    
    clear_outliers(L)
    L = scipy.signal.savgol_filter(L, 5, 3)
    plt.plot(Index, L, label=str(I))
#♣plt.ylim((40, 65))
#plt.xlim((222, 227))
plt.legend()
plt.title(f"plot of the integral of Current (A) for all cycles from {a} to {b-1} as a function of time (min)")
plt.xlabel("cycle")
plt.ylabel("Integral of I (A*min)")
plt.show()
"""

"""""""""

a = 0
b = len(mat["batch"]["cycles"])

step=1

for I in range(a, b, step):
    Index = []
    df = create_df(I)
    L = []
    for i in range(int(min(df.cycle)), int(max(df.cycle)), cycle_step):
        Int = scipy.integrate.trapezoid(df.loc[df.cycle==i]["I"], df.loc[df.cycle==i]["t"])
        if Int>30:
            L.append(Int)
            Index.append(i)
    
    clear_outliers2(L)
    L = scipy.signal.savgol_filter(L, 5, 3)
    plt.plot(Index, L, label=str(I))
    
#♣plt.ylim((40, 65))
#plt.xlim((222, 227))
plt.legend()
plt.title(f"plot of the integral of Current (A) for all cycles from {a} to {b-1} as a function of time (min)")
plt.xlabel("cycle")
plt.ylabel("Integral of I (A*min)")
plt.show()









Step = 10
a=0
b=len(mat["batch"]["cycles"])
b = Step-1
while a<b:
    for I in range(a, min(a+Step, b)):
        Index = []
        df = create_df(I)
        L = []
        for i in range(int(min(df.cycle)), int(max(df.cycle)), cycle_step):
            Int = scipy.integrate.trapezoid(df.loc[df.cycle==i]["I"], df.loc[df.cycle==i]["t"])
            if Int>30:
                L.append(Int)
                Index.append(i)
        clear_outliers(L)
        L = scipy.signal.savgol_filter(L, 5, 3)
        Index = np.array(Index)
        L = np.array(L)
        plt.plot(Index, L, label=str(I))
    plt.ylim((40, 65))


    
    plt.legend()
    plt.title(f"S_plot of the integral of current (A) for all cycles from {a} to {min(a+Step, b)} as a function of time (min)")
    plt.xlabel("time (min)")
    plt.ylabel("Integral of I (A*min)")
    plt.savefig(f"Plot of smoothed Integrals current {a} to {min(a+Step, b)}", dpi='figure', format=None)
    a+=Step
    plt.show()
    
    
    
    
    

Step = 10
a=0
b=len(mat["batch"]["cycles"])
b = Step-1
while a<b:
    for I in range(a, min(a+Step, b)):
        Index = []
        df = create_df(I)
        L = []
        for i in range(int(min(df.cycle)), int(max(df.cycle)), cycle_step):
            Int = scipy.integrate.trapezoid(df.loc[df.cycle==i]["I"], df.loc[df.cycle==i]["t"])
            if Int>30:
                L.append(Int)
                Index.append(i)
        Index = np.array(Index)
        L = np.array(L)
        plt.plot(Index, L, label=str(I))
    plt.ylim((40, 65))
    plt.legend()
    plt.title(f"plot of the integral of current (A) for all cycles from {a} to {min(a+Step, b)} as a function of time (min)")
    plt.xlabel("time (min)")
    plt.ylabel("Integral of I (A*min)")
    plt.savefig(f"Plot Integrals current {a} to {min(a+Step, b)}", dpi='figure', format=None)
    a+=Step
    plt.show()
















for i in range(int(min(df.cycle)), int(max(df.cycle)), cycle_step):
    plt.plot(df.loc[df.cycle==i]["t"],df.loc[df.cycle==i]["T"])
plt.xlim((0, 50))
plt.title("Temperature (°C) as a function of time (min)")
plt.xlabel("Time (min)")
plt.ylabel("Temperature (°C)")
plt.show()

"""""""""
"""
for i in range(0, int(max(df.cycle)), cycle_step):
    plt.plot(df.loc[df.cycle==i]["t"],np.cumsum(df.loc[df.cycle==i]["T"]))
plt.xlim((0, 50))
plt.title("Cumulative Temperature (°C) as a function of time (min)")
plt.xlabel("Time (min)")
plt.ylabel("Cumulative Temperature (°C)")

plt.show()
"""














"""""""""

Step = 10
a=0
b=len(mat["batch"]["cycles"])
b = Step-1
while a<b:
    for I in range(a, min(a+Step, b)):
        Index = []
        df = create_df(I)
        L = []
        for i in range(int(min(df.cycle)), int(max(df.cycle)), cycle_step):
            Int = scipy.integrate.trapezoid(df.loc[df.cycle==i]["T"], df.loc[df.cycle==i]["t"])
            L.append(Int)
            Index.append(i)
        clear_outliers(L)
        plt.plot(Index, L, label=str(I))
    plt.legend()
    plt.title(f"plot of the integral of Temperature (°C*min) for all cycles from {a} to {min(a+Step, b)} as a function of time (min)")
    plt.xlabel("cycle ")
    plt.ylabel("Integral of Temperature (°C*min)")
    plt.savefig(f"Plot Integrals Temperature {a} to {min(a+Step, b)}", dpi='figure', format=None)
    a+=Step
    plt.show()


Step = 10
a=0
b=len(mat["batch"]["cycles"])
b = Step-1
while a<b:
    for I in range(a, min(a+Step, b)):
        Index = []
        df = create_df(I)
        L = []
        for i in range(int(min(df.cycle)), int(max(df.cycle)), cycle_step):
            Int = scipy.integrate.trapezoid(df.loc[df.cycle==i]["T"], df.loc[df.cycle==i]["t"])
            L.append(Int)
            Index.append(i)
        plt.plot(Index, L, label=str(I))
    plt.legend()
    plt.title(f"plot of the integral of Temperature (°C*min) for all cycles from {a} to {min(a+Step, b)} as a function of time (min)")
    plt.xlabel("cycle ")
    plt.ylabel("Integral of Temperature (°C*min)")
    plt.savefig(f"Plot Integrals Temperature {a} to {min(a+Step, b)}", dpi='figure', format=None)
    a+=Step
    plt.show()





for i in range(0, int(max(df.cycle)), cycle_step):
    plt.plot(df.loc[df.cycle==i]["t"],df.loc[df.cycle==i]["V"])
plt.xlim((0, 50))
plt.title("Voltage (V) as a function of time (min)")
plt.xlabel("Time (min)")
plt.ylabel("Voltage (V)")
plt.show()

for i in range(0, int(max(df.cycle)), cycle_step):
    plt.plot(df.loc[df.cycle==i]["t"],np.cumsum(df.loc[df.cycle==i]["V"]))
plt.xlim((0, 50))
plt.title("Cumulative voltage (V) as a function of time (min)")
plt.xlabel("Time (min)")
plt.ylabel("Cumulative voltage (V)")

plt.show()







Step = 10
a=0
b=len(mat["batch"]["cycles"])
b = Step-1
while a<b:
    for I in range(a, min(a+Step, b)):
        Index = []
        df = create_df(I)
        L = []
        for i in range(int(min(df.cycle)), int(max(df.cycle)), cycle_step):
            Int = scipy.integrate.trapezoid(df.loc[df.cycle==i]["V"], df.loc[df.cycle==i]["t"])
            L.append(Int)
            Index.append(i)
        clear_outliers(L)
        L = scipy.signal.savgol_filter(L, 5, 3)
        plt.plot(Index, L, label=str(I))
    plt.legend()
    plt.ylim((90, 220))
    plt.title(f"plot of the integral of Voltage (V*min) for all cycles from {a} to {min(a+Step, b)} as a function of time (min)")
    plt.xlabel("cycle")
    plt.ylabel("Integral of Voltage (V*min)")
    plt.savefig(f"Plot Integrals Voltage {a} to {min(a+Step, b)}", dpi='figure', format=None)
    a+=Step
    plt.show()


Step = 10
a=0
b=len(mat["batch"]["cycles"])
b = Step-1
while a<b:
    for I in range(a, min(a+Step, b)):
        Index = []
        df = create_df(I)
        L = []
        for i in range(10, int(max(df.cycle)), cycle_step):
            Int = scipy.integrate.trapezoid(df.loc[df.cycle==i]["V"], df.loc[df.cycle==i]["t"])
            L.append(Int)
            Index.append(i)
        plt.plot(Index, L, label=str(I))
    plt.legend()
    plt.ylim((90, 220))
    plt.title(f"plot of the integral of Voltage (V*min) for all cycles from {a} to {min(a+Step, b)} as a function of time (min)")
    plt.xlabel("cycle")
    plt.ylabel("Integral of Voltage (V*min)")
    plt.savefig(f"Plot Integrals Voltage {a} to {min(a+Step, b)}", dpi='figure', format=None)
    a+=Step
    plt.show()



"""""""""






"""
for i in range(0, int(max(df.cycle)), cycle_step):
    plt.plot(df.loc[df.cycle==i]["I"],df.loc[df.cycle==i]["V"])
plt.title("Voltage (V) as a function of Current (A)")
plt.xlabel("Current (A)")
plt.ylabel("Voltage (V)")
plt.show()
"""


"""""""""

i=25
for j in range(min(df.cycle), min(df.cycle)+500, 10):
    dx = dx_(np.array(df.loc[df.cycle==j]["T"]))
    dy = dx_(np.array(df.loc[df.cycle==j]["V"]))
    plt.plot(np.array(df.loc[df.cycle==j]["T"])[1:], dy/dx)
    plt.ylim((-2000, 2000))
    plt.xlabel("T")
    plt.ylabel("dV/dT")
plt.title("dV/dT as a function of T")
plt.show()

i=25
for j in range(min(df.cycle), min(df.cycle)+500, 10):
    dx = dx_(np.array(df.loc[df.cycle==j]["T"]))
    dy = dx_(np.array(df.loc[df.cycle==j]["V"])*np.array(df.loc[df.cycle==j]["I"]))
    plt.plot(np.array(df.loc[df.cycle==j]["T"])[1:], dy/dx)
    plt.xlabel("T")
    plt.ylabel("d(VI)/dT")
plt.title("d(VI)/dT as a function of T")
plt.show()



for j in range(min(df.cycle), max(df.cycle), 200):
    dx = dx_(np.cumsum(np.array(df.loc[df.cycle==j]["T"])))
    dy = dx_(np.cumsum(np.array(df.loc[df.cycle==j]["V"])))
    Y = dy/dx
    X = np.array(df.loc[df.cycle==j]["I"])[1:]
    X = np.cumsum(X)
    clear_outliers(Y, 3, 10, 1)
    plt.plot(X, Y, label = str(j))
    plt.xlabel("I")
    plt.ylabel("dCumSum(V)/dCumSum(I)")
    
    
    
i=25
for j in range(min(df.cycle), min(df.cycle)+500, 10):
    dx = dx_(np.array(df.loc[df.cycle==j]["T"]))
    dy = dx_(np.array(df.loc[df.cycle==j]["V"]))
    plt.plot(np.array(df.loc[df.cycle==j]["t"])[1:], dy/dx)
    plt.ylim((-2000, 2000))
    plt.xlabel("time")
    plt.ylabel("dV/dT")
plt.title("dV/dT as a function of t")
plt.show()


for j in range(min(df.cycle), min(df.cycle)+500, 10):
    dx = dx_(np.array(df.loc[df.cycle==j]["T"])*2)
    dy = dx_(np.array(df.loc[df.cycle==j]["V"]))
    dy2 = dx_(np.array(df.loc[df.cycle==j]["I"]))
    plt.plot(np.array(df.loc[df.cycle==j]["T"])[1:], dy*dy2/(dx*2))
    plt.xlabel("T")
    plt.ylabel("dVdI/dT²")
    plt.ylim((-2000, 2000))
    plt.xlim((32, 34))
plt.title("dVdI/dT² as a function of T")
plt.show()


for j in range(min(df.cycle), max(df.cycle), 200):
    dx = dx_(np.cumsum(np.array(df.loc[df.cycle==j]["T"])))
    dy = dx_(np.cumsum(np.array(df.loc[df.cycle==j]["V"])))
    Y = dy/dx
    X = np.array(df.loc[df.cycle==j]["t"])[1:]
    #X = np.cumsum(X)
    clear_outliers(Y, 3, 10, 1)
    plt.plot(X, Y, label = str(j))
    plt.xlabel("Temperature")
    plt.ylabel("dCumSum(V)/dCumSum(T)")
    


plt.legend()
plt.title("dCumSum(V)/dCumSum(T) as a function of time")
plt.show()
print(max(df.cycle))


Step = 10
a=0
b=len(mat["batch"]["cycles"])
while a<b:
    for I in range(a, min(a+Step, b)):
        Index = []
        df = create_df(I)
        L = []
        for i in range(10, int(max(df.cycle)), cycle_step):
            Int = scipy.integrate.trapezoid(df.loc[df.cycle==i]["V"], df.loc[df.cycle==i]["I"])
            L.append(Int)
            Index.append(i)
        Index = np.array(Index)
        L = np.array(L)
        L = smothing(Index, L)
        plt.plot(Index, L, label=str(I))
    plt.legend()
    plt.title(f"plot of the integral of Voltage (V) for all cycles from {a} to {min(a+Step, b)} as a function of Current (A)")
    plt.xlabel("Cycle")
    plt.ylabel("Integral of V (V*A)")
    plt.savefig(f"Plot of smoothed Integrals Voltage with respect to current {a} to {min(a+Step, b)}", dpi='figure', format=None)
    a+=Step
    plt.show()


#do not use, too random and seems constant 
Step = 10
a=0
b=len(mat["batch"]["cycles"])
while a<b:
    for I in range(a, min(a+Step, b)):
        Index = []
        df = create_df(I)
        L = []
        for i in range(10, int(max(df.cycle)), cycle_step):
            Int = scipy.integrate.trapezoid(df.loc[df.cycle==i]["V"], df.loc[df.cycle==i]["T"])
            L.append(Int)
            Index.append(i)
        Index = np.array(Index)
        L = np.array(L)
        plt.plot(Index, L, label=str(I))
    plt.legend()
    plt.title(f"plot of the integral of Voltage for all cycles from {a} to {min(a+Step, b)} as a function of Temperature")
    plt.xlabel("Cycle")
    plt.ylabel("Integral of V")
    plt.savefig(f"Plot of Integrals Voltage with respect to Temperature {a} to {min(a+Step, b)}", dpi='figure', format=None)
    a+=Step
    plt.show()


Step = 10
a=0
b=len(mat["batch"]["cycles"])
while a<b:
    for I in range(a, min(a+Step, b)):
        Index = []
        df = create_df(I)
        L = []
        for i in range(10, int(max(df.cycle)), cycle_step):
            Int = scipy.integrate.trapezoid(df.loc[df.cycle==i]["V"], df.loc[df.cycle==i]["T"])
            L.append(Int)
            Index.append(i)
        Index = np.array(Index)
        L = np.array(L)
        L = smothing(Index, L)
        plt.plot(Index, L, label=str(I))
    plt.legend()
    plt.title(f"plot of the integral of Voltage for all cycles from {a} to {min(a+Step, b)} as a function of Temperature")
    plt.xlabel("Cycle")
    plt.ylabel("Integral of V")
    plt.savefig(f"Plot of smoothed Integrals Voltage with respect to Temperature {a} to {min(a+Step, b)}", dpi='figure', format=None)
    a+=Step
    plt.show()
    
    
"""""""""
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    
"""
    ------------------------------------------------
            CREATION OF THE DATA SET:
    ------------------------------------------------
"""




"""
for i in range(50, 51, cycle_step):
    plt.plot(df.loc[df.cycle==i]["t"],df.loc[df.cycle==i]["V"])
plt.title("V(t)")
plt.xlabel("Time (min)")
plt.ylabel("V")

plt.show()

#?# Qc gets bigger than 1???
for i in range(50, 51, cycle_step):
    plt.plot(df.loc[df.cycle==i]["t"],df.loc[df.cycle==i]["Qc"])
plt.title("Qc(t)")
plt.xlabel("Time (min)")
plt.ylabel("Qc")

plt.show()

"""





N_cells = 10#len(mat["batch"]["cycles"])


#Get Integrales of Intensity:
"""
Integrale_Intensity = [0 for _ in range(N_cells)]
    
for I in range(0, N_cells, 1):
    Index = []
    df = create_df(I, ["I", "t"])
    L = []
    for i in range(1, int(max(df.cycle)), cycle_step):
        Int = scipy.integrate.trapezoid(df.loc[df.cycle==i]["I"], df.loc[df.cycle==i]["t"])
        if Int>30:
            L.append(Int)
            Index.append(i)
    
    clear_outliers(L)
    L = scipy.signal.savgol_filter(L, 5, 3)
    Integrale_Intensity[I]=L
"""



'''
INFOS:
    C is constant: C = 1.1A
    
    The batterie chages at:
        C1 until Q1
        Then constant Voltage
        Then C2 untile Qc=0.8
        Then 1C
    /!\ C1 C2 and Q1 are specific to the batterie
        
HOW WE GET C1, C2 and Q1:
    --> We only work on the 50-th cycle for I and Qc (this is because we discard the first ones, and we
        still have some outliers at first, so it might be better working with those ones at first)
    C1:
        -We get the first value of intensity => wont work, isn't instantly there :/
        ->We get the biggest value on the sub part of I with is monotonous (increasing)
    Q1:
        ->Use what we did to get C1
    C2:
        ->We get I when Qc is close to 0.79 (not 0.8 since rught after it change, so it's safer like so)
'''

'''
TO do:
    Remove cycle 41, 36, 35, 23, 22, 38, 37 (weird Current behaviour)
    Keep C1=C2?
'''

"""
#Get C1, Q1 and C2:

Cycle_worked_on = 50
C1s = [0 for _ in range(N_cells)]
C2s = [0 for _ in range(N_cells)]
Q1s = [0 for _ in range(N_cells)]
    
for Index in range(0, N_cells, 1):
    df = create_df(Index, ["I", "Qc", "t"])
    I = list(df.loc[df.cycle==Cycle_worked_on]["I"])
    Qc = list(df.loc[df.cycle==Cycle_worked_on]["Qc"])
    t = list(df.loc[df.cycle==Cycle_worked_on]["t"])
    L = [(x, y, z) for x, y, z in zip(t, Qc, I)]
    L = sorted(L, key = lambda X: X[0])
    for i in range(len(L)):
        t[i], Qc[i], I[i] = L[i]
    
    #Now we have t, Qc and I being sorted with respect to t
    """""""
    plt.plot(list(range(len(I))), I)
    plt.title(str(Index))
    plt.show()
    
    plt.plot(list(range(len(Qc))), Qc)
    plt.title(str(Index))
    plt.show()
    """""""
    #Get C1 and Q1
    C1 = 0
    for i in range(len(I)):
        #-0.05 so that is stabe under small changes of I
        if I[i]>=C1-0.05:
            C1=I[i]
        else:
            Q1 = Qc[i]
            break
    
    #Get C2:
    for i in range(len(Qc)):
        if Qc[i]>0.79:
            C2 = I[i]
            break
    C1s[Index], Q1s[Index], C2s[Index] = C1, Q1, C2
"""
#We get the y value (what we want to predict)
"""
Y = []
    
for I in range(0, N_cells, 1):
    for i in range(int(min(df.cycle)), int(max(df.cycle)), cycle_step):
        Y.append(int(max(df.cycle))-i)
Y = np.array(Y)
print(Y.shape)

"""


#Now lets do everything:

banned_cycles = [41, 36, 35, 23, 22, 38, 37]




"""

S =  "EL150800460514 Prim. Test 2017-05-12 1852 3.6C(80%)-3.6C EL150800460486 Train 2017-05-12 2160 3.6C(80%)-3.6C EL150800460623 Prim. Test 2017-05-12 2237 3.6C(80%)-3.6C EL150800464977 Train 2017-05-12 1434 4C(80%)-4C EL150800464865 Prim. Test 2017-05-12 1709 4C(80%)-4C EL150800464883 Train 2017-05-12 1074 4.4C(80%)-4.4C EL150800463886 Prim. Test 2017-05-12 636 4.8C(80%)-4.8C EL150800465027 Train 2017-05-12 870 4.8C(80%)-4.8C EL150800460468 Prim. Test 2017-05-12 1054 5.4C(40%)-3.6C EL150800463882 Train 2017-05-12 788 5.4C(40%)-3.6C EL150800463838 Prim. Test 2017-05-12 880 5.4C(50%)-3C EL150800453113 Train 2017-05-12 719 5.4C(50%)-3C EL150800460653 Prim. Test 2017-05-12 862 5.4C(50%)-3.6C EL150800460522 Train 2017-05-12 857 5.4C(50%)-3.6C EL150800453240 Prim. Test 2017-05-12 691 5.4C(60%)-3C EL150800464881 Train 2017-05-12 788 5.4C(60%)-3C EL150800464002 Prim. Test 2017-05-12 534 5.4C(60%)-3.6C EL150800463871 Train 2017-05-12 559 5.4C(60%)-3.6C EL150800460477 Prim. Test 2017-05-12 1014 5.4C(70%)-3C EL150800460630 Train 2017-05-12 1017 5.4C(70%)-3C EL150800460617 Prim. Test 2017-05-12 854 5.4C(80%)-5.4C EL150800460659 Train 2017-05-12 870 5.4C(80%)-5.4C EL150800463198 Prim. Test 2017-05-12 842 6C(30%)-3.6C EL150800460507 Train 2017-05-12 860 6C(30%)-3.6C EL150800460644 Prim. Test 2017-05-12 917 6C(40%)-3C EL150800460615 Train 2017-05-12 709 6C(40%)-3C EL150800460481 Prim. Test 2017-05-12 876 6C(40%)-3.6C EL150800460640 Train 2017-05-12 731 6C(40%)-3.6C EL150800460436 Prim. Test 2017-05-12 757 6C(50%)-3C EL150800460525 Train 2017-05-12 742 6C(50%)-3C EL150800460622 Prim. Test 2017-05-12 703 6C(50%)-3.6C EL150800460506 Train 2017-05-12 704 6C(50%)-3.6C EL150800460601 Prim. Test 2017-05-12 648 6C(60%)-3C EL150800453773 Train 2017-05-12 617 6C(60%)-3C EL150800460642 Prim. Test 2017 -05 -12 625 7C(30%) -3.6C EL150800463229 Train 2017 -05 -12 966 7C(30%) -3.6C EL150800460647 Prim. Test 2017 -05 -12 1051 7C(40%) -3C EL150800460618 Train 2017 -05 -12 702 7C(40%) -3C EL150800460636 Prim. Test 2017 -05 -12 651 7C(40%) -3.6C EL150800460485 Train 2017 -05 -12 616 7C(40%) -3.6C EL150800460656 Prim. Test 2017 -05 -12 599 8C(15%) -3.6C EL150800460518 Train 2017 -06 -30 300 1C(4%) -6C EL150800460605 Prim. Test 2017 -06 -30 148 2C(10%) -6C EL150800460602 Train 2017 -06 -30 438 2C(2%) -5C EL150800460673 Prim. Test 2017 -06 -30 335 2C(7%) -5.5C EL150800460655 Train 2017 -06 -30 444 3.6C(22%) -5.5C EL150800460635 Prim. Test 2017 -06 -30 480 3.6C(2%) -4.85C EL150800460634 Train 2017 -06 -30 511 3.6C(30%) -6C EL150800460451 Prim. Test 2017 -06 -30 561 3.6C(9%) -5C EL150800460466 Train 2017 -06 -30 477 4C(13%) -5C EL150800460510 Prim. Test 2017 -06 -30 458 4C(31%) -5C EL150800463208 Train 2017 -06 -30 483 4C(40%) -6C EL150800460449 Prim. Test 2017 -06 -30 485 4C(4%) -4.85C EL150800460478 Train 2017 -06 -30 494 4.4C(24%) -5C EL150800460480 Prim. Test 2017 -06 -30 487 4.4C(47%) -5.5C EL150800440551 Train 2017 -06 -30 461 4.4C(55%) -6C EL150800460652 Prim. Test 2017 -06 -30 502 4.4C(8%) -4.85C EL150800460603 Train 2017 -06 -30 489 4.65C(19%) -4.85C EL150800463245 Prim. Test 2017 -06 -30 513 4.65C(44%) -5C EL150800460501 Train 2017 -06 -30 527 4.65C(69%) -6C EL150800460597 Prim. Test 2017 -06 -30 495 4.8C(80%) -4.8C EL150800460611 Train 2017 -06 -30 461 4.8C(80%) -4.8C EL150800460596 Prim. Test 2017 -06 -30 471 4.8C(80%) -4.8C EL150800460614 Train 2017 -06 -30 468 4.9C(27%) -4.75C EL150800460610 Prim. Test 2017 -06 -30 509 4.9C(61%) -4.5C EL150800460604 Train 2017 -06 -30 498 4.9C(69%) -4.25C EL150800460527 Prim. Test 2017 -06 -30 481 5.2C(10%) -4.75C EL150800460608 Train 2017 -06 -30 492 5.2C(37%) -4.5C EL150800460631 Prim. Test 2017 -06 -30 519 5.2C(50%) -4.25C EL150800460641 Train 2017 -06 -30 520 5.2C(58%) -4C EL150800460492 Prim. Test 2017 -06 -30 499 5.2C(66%) -3.5C EL150800460628 Train 2017 -06 -30 463 5.2C(71%) -3C EL150800460528 Prim. Test 2017 -06 -30 535 5.6C(25%) -4.5C EL150800460511 Train 2017 -06 -30 478 5.6C(38%) -4.25C EL150800460649 Prim. Test 2017 -06 -30 465 5.6C(47%) -4C EL150800460627 Train 2017 -06 -30 459 5.6C(58%) -3.5C EL150800460526 Prim. Test 2017 -06 -30 499 5.6C(5%) -4.75C EL150800460513 Train 2017 -06 -30 429 5.6C(65%) -3C EL150800460473 Prim. Test 2017 -06 -30 466 6C(20%) -4.5C EL150800460498 Train 2017 -06 -30 462 6C(31%) -4.25C EL150800460474 Prim. Test 2017 -06 -30 457 6C(40%) -4C EL150800460613 Train 2017 -06 -30 487 6C(4%) -4.75C EL150800460678 Prim. Test 2017 -06 -30 429 6C(52%) -3.5C EL150800460599 Prim. Test 2017 -06 -30 713 6C(60%) -3C EL150800737329 Sec. test 2018 -04 -12 1009 5C(67%) -4C EL150800737313 Sec. test 2018 -04 -12 1063 5.3C(54%) -4C EL150800737280 Sec. test 2018 -04 -12 1115 5.6C(36%) -4.3C EL150800739476 Sec. test 2018 -04 -12 1048 5.6C(19%) -4.6C EL150800737213 Sec. test 2018 -04 -12 828 5.6C(36%) -4.3C EL150800737229 Sec. test 2018 -04 -12 667 3.7C(31%) -5.9C EL150800737307 Sec. test 2018 -04 -12 1836 4.36C(80%) -4.36 C EL150800737233 Sec. test 2018 -04 -12 828 5C(67%) -4C EL150800737270 Sec. test 2018 -04 -12 1039 5.3C(54%) -4C EL150800737277 Sec. test 2018 -04 -12 1078 4.36C(80%) -4.36 C EL150800737276 Sec. test 2018 -04 -12 817 5.6C(19%) -4.6C EL150800737314 Sec. test 2018 -04 -12 932 5.6C(36%) -4.3C EL150800737319 Sec. test 2018 -04 -12 816 5.6C(19%) -4.6C EL150800737387 Sec. test 2018 -04 -12 858 5.6C(36%) -4.3C EL150800737386 Sec. test 2018 -04 -12 876 5.9C(15%) -4.6C EL150800737345 Sec. test 2018 -04 -12 1638 4.36C(80%) -4.36 C EL150800737378 Sec. test 2018 -04 -12 1315 5.3C(54%) -4C EL150800737274 Sec. test 2018 -04 -12 1146 5.6C(19%) -4.6C EL150800737275 Sec. test 2018 -04 -12 1155 5.6C(36%) -4.3C EL150800737315 Sec. test 2018 -04 -12 813 5C(67%) -4C EL150800737366 Sec. test 2018 -04 -12 772 3.7C(31%) -5.9C EL150800737285 Sec. test 2018 -04 -12 1002 5.9C(60%) -3.1C EL150800737368 Sec. test 2018 -04 -12 825 5C(67%) -4C EL150800737259 Sec. test 2018 -04 -12 989 5.3C(54%) -4C EL150800737287 Sec. test 2018 -04 -12 1028 5.6C(19%) -4.6C EL150800737251 Sec. test 2018 -04 -12 850 5.6C(36%) -4.3C EL150800737234 Sec. test 2018 -04 -12 541 3.7C(31%) -5.9C EL150800737320 Sec. test 2018 -04 -12 858 5.9C(15%) -4.6C EL150800737380 Sec. test 2018 -04 -12 935 5.3C(54%) -4C EL150800737279 Sec. test 2018 -04 -12 731 5.9C(60%) -3.1C EL150800737304 Sec. test 2018 -04 -12 1284 5C(67%) -4C EL150800737350 Sec. test 2018 -04 -12 1158 5.3C(54%) -4C EL150800739477 Sec. test 2018 -04 -12 1093 5.6C(19%) -4.6C EL150800737365 Sec. test 2018 -04 -12 923 5.6C(36%) -4.3C EL150800737334 Sec. test 2018 -04 -12 1935 5C(67%) -4C EL150800737361 Sec. test 2018 -04 -12 1156 5.3C(54%) -4C EL150800737390 Sec. test 2018 -04 -12 796 5.6C(19%) -4.6C EL150800739495 Sec. test 2018 -04 -12 786 5.6C(36%) -4.3C EL150800737369 Sec. test 2018 -04 -12 940 5.3C(54%) -4C EL150800739484 Sec. test 2018 -04 -12 1801 4.36C(80%) -4.36 C"
PARAM = S.split()
d = 0
for i in range(len(PARAM)):
	i = i-d
	if "EL1" in PARAM[i]:
		PARAM.pop(i)
		d+=1
	elif "2017" in PARAM[i]:
		PARAM.pop(i)
		d+=1
	elif "Prim" in PARAM[i]:
		PARAM.pop(i)
		d+=1
	elif "C" not in PARAM[i]:
		PARAM.pop(i)
		d+=1
	elif "C"!=PARAM[i][-1]:
		PARAM[i] = PARAM[i] + PARAM[i+1]
		PARAM.pop(i+1)
		d+=1

PARAM = [x.replace(")-", "(") for x in PARAM]
PARAM = [x.split("(") for x in PARAM]

PARAM = [[float(x[0][:-1]), float(x[1][:-1]), float(x[2][:-1])] for x in PARAM]
print(len(PARAM))
PARAM = PARAM[N_cells:N_cells*2]
##WRONG!!!!! THERE FOR NOW BUT TO CHANGE!!!!!






Integrale_Intensity = []
Integrale_Voltage = []
Integrale_Temperature = []
C1s = []
C2s = []
Q1s = []
Cycle_worked_on = 50
Y = []
for I in range(0, N_cells, 1):
    if I in banned_cycles:
        continue
    else:
        IN = I
        Index = []
        df = create_df(I, ["I", "t", "Qc", "T", "V"])
        C1, Q1, C2 = PARAM[I]
    
    
        L = []
        L_T = []
        L_V = []
        for i in range(int(min(df.cycle)), int(max(df.cycle)), cycle_step):
            Int = scipy.integrate.trapezoid(df.loc[df.cycle==i]["I"], df.loc[df.cycle==i]["t"])
            L.append(Int)
            Int_T = scipy.integrate.trapezoid(df.loc[df.cycle==i]["T"], df.loc[df.cycle==i]["t"])
            L_T.append(Int_T)
            Int_V = scipy.integrate.trapezoid(df.loc[df.cycle==i]["V"], df.loc[df.cycle==i]["t"])
            L_V.append(Int_V)
            
            
            Y.append(mat[mat["batch"]["cycle_life"][IN][0]][0][0])
            C1s.append(C1)
            C2s.append(C2)
            Q1s.append(Q1)
    
        clear_outliers(L)
        L = scipy.signal.savgol_filter(L, 5, 3)
        [Integrale_Intensity.append(x) for x in L]
        [Integrale_Temperature.append(x) for x in L_T]
        [Integrale_Voltage.append(x) for x in L_V]
        

X = [[C1s[i], C2s[i], Q1s[i], Integrale_Intensity[i], Integrale_Temperature[i], Integrale_Voltage[i], Y[i]] for i in range(len(Y))]

X = pd.DataFrame(X, columns = ["C1", "C2", "Q1", "Int_I_t", "Int_T_t", "Int_V_t", "Y"])
print(X)
X.to_csv("input.csv")
"""



X = pd.read_csv("input.csv", index_col=0)
print(X)




"""
    ------------------------------------------------
            THE NORMAL BNN
    ------------------------------------------------
"""



tf.keras.utils.set_random_seed(42)
LR=0.01
model = tf.keras.Sequential()
for i in [8, 8]:
    model.add(layers.Dense(units=i, activation ="relu"))
model.add(layers.Dense(units=1, activation ="relu"))

model.compile(
    loss='mean_squared_error', 
    optimizer=tf.optimizers.Adam(LR), 
    metrics=['mean_squared_error'])


num_epochs = 10
train_size = int(len(X)*0.8)

n = len(X)
Train_indexes = random.sample(list(range(n)), train_size)
X_train, X_test, Y_train, Y_test = [], [], [], []
X = X.T
for i in range(n):
    l = list(X[i])
    y = l.pop(-1)
    if i in Train_indexes:
        X_train.append(l)
        Y_train.append(y)
    else:
        X_test.append(l)
        Y_test.append(y)

model.fit(X_train, Y_train, epochs=num_epochs, batch_size=32)


for i in range(0, len(X_test), len(X_test)//10):
    print(model.predict([X_test[i]])[0][0], Y_test[i])



plt.plot(list(range(len(X_test))), model.predict(X_test), c = "blue")
plt.plot(list(range(len(X_test))), Y_test, c = "red")
plt.show()








tf.keras.utils.set_random_seed(42)
LR=0.01
model = tf.keras.Sequential()
for i in [8, 8]:
    model.add(layers.Dense(units=i, activation ="relu"))
model.add(layers.Dense(units=1, activation ="relu"))

model.compile(
    loss='mean_squared_error', 
    optimizer=tf.optimizers.Adam(LR), 
    metrics=['mean_squared_error'])


num_epochs = 1000
train_size = int(len(X)*0.8)

n = len(X)

model.fit(X_train, Y_train, epochs=num_epochs, batch_size=32)


for i in range(0, len(X_test), len(X_test)//10):
    print(model.predict([X_test[i]])[0][0], Y_test[i])



plt.plot(list(range(len(X_test))), model.predict(X_test), c = "blue")
plt.plot(list(range(len(X_test))), Y_test, c = "red")
plt.show()







tf.keras.utils.set_random_seed(42)
LR=0.01
model = tf.keras.Sequential()
for i in [8, 8, 8, 8]:
    model.add(layers.Dense(units=i, activation ="relu"))
model.add(layers.Dense(units=1, activation ="relu"))

model.compile(
    loss='mean_squared_error', 
    optimizer=tf.optimizers.Adam(LR), 
    metrics=['mean_squared_error'])


num_epochs = 1000
train_size = int(len(X)*0.8)

n = len(X)


model.fit(X_train, Y_train, epochs=num_epochs, batch_size=32)


for i in range(0, len(X_test), len(X_test)//10):
    print(model.predict([X_test[i]])[0][0], Y_test[i])



plt.plot(list(range(len(X_test))), model.predict(X_test), c = "blue")
plt.plot(list(range(len(X_test))), Y_test, c = "red")
plt.show()


"""
    ------------------------------------------------
            THE PROBA BNN
    ------------------------------------------------
"""


hidden_units = [8, 8]
LR = 0.01
tf.keras.utils.set_random_seed(42)

def prior(kernel_size, bias_size, dtype=None):
    n = kernel_size + bias_size
    prior_model = keras.Sequential(
        [
            tfp.layers.DistributionLambda(
                lambda t: tfp.distributions.MultivariateNormalDiag(
                    loc=tf.zeros(n), scale_diag=tf.ones(n)
                )
            )
        ]
    )
    return prior_model




def posterior(kernel_size, bias_size, dtype=None):
    n = kernel_size + bias_size
    posterior_model = keras.Sequential(
        [
            tfp.layers.VariableLayer(
                tfp.layers.MultivariateNormalTriL.params_size(n), dtype=dtype
            ),
            tfp.layers.MultivariateNormalTriL(n),
        ]
    )
    return posterior_model



def create_probablistic_bnn_model(train_size):
    model = tf.keras.Sequential()
    model.add(layers.BatchNormalization())
    for units in hidden_units:
        model.add(tfp.layers.DenseVariational(
            units=units,
            make_prior_fn=prior,
            make_posterior_fn=posterior,
            kl_weight=1 / train_size,
            activation="relu",
            ))


    # Create a probabilisticå output (Normal distribution), and use the `Dense` layer
    # to produce the parameters of the distribution.
    # We set units=2 to learn both the mean and the variance of the Normal distribution.
    model.add(layers.Dense(units=2))
    model.add(tfp.layers.IndependentNormal(1))
    
    return model






def negative_loglikelihood(targets, estimated_distribution):
    return -estimated_distribution.log_prob(targets)





def run_experiment(model, loss, X_train, X_test, Y_train, Y_test):

    model.compile(
        optimizer=keras.optimizers.RMSprop(learning_rate=LR),
        loss=loss,
        metrics=[keras.metrics.RootMeanSquaredError()],
    )

    print("Start training the model...")
    model.fit(X_train, Y_train, epochs=num_epochs, batch_size=1)
    print("Model training finished.")
    _, rmse = model.evaluate(X_train, Y_train, verbose=0)
    print(f"Train RMSE: {round(rmse, 3)}")

    print("Evaluating model performance...")
    _, rmse = model.evaluate(X_test, Y_test, verbose=0)
    print(f"Test RMSE: {round(rmse, 3)}")

num_epochs = 10
train_size = len(X)*0.8
prob_bnn_model = create_probablistic_bnn_model(train_size)



run_experiment(prob_bnn_model, negative_loglikelihood,  X_train, X_test, Y_train, Y_test)


X_test = np.array(X_test)
X_test.reshape((len(X_test), len(X_test[0]), 1))
prediction_distribution = prob_bnn_model(X_test[0:len(X_test):len(X_test)//10])
prediction_mean = prediction_distribution.mean().numpy().tolist()
prediction_stdv = prediction_distribution.stddev().numpy()

# The 95% CI is computed as mean ± (1.96 * stdv)
upper = (prediction_mean + (1.96 * prediction_stdv)).tolist()
lower = (prediction_mean - (1.96 * prediction_stdv)).tolist()
prediction_stdv = prediction_stdv.tolist()
sample = len(prediction_mean)
for idx in range(sample):
    print(
        f"Prediction mean: {round(prediction_mean[idx][0], 2)}, "
        f"stddev: {round(prediction_stdv[idx][0], 2)}, "
        f"95% CI: [{round(upper[idx][0], 2)} - {round(lower[idx][0], 2)}]"
        f" - Actual: {Y_test[idx]}"
    )











