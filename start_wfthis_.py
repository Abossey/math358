import numpy as np
import matplotlib.pyplot as plt
#from scipy.optimize import fsolve
from scipy.optimize import newton

def euler(f,y0,a,h):
    """ Calculates the Euler solution of the IVP
     y'=f(t,y), with y(a)=y0, a[0]<=t<=a[1]. Using the 
     explicit Euler formula. If h<=1, then h is the step size
     else h is the number of nodes to use """
    
    if h <= 1:
	t = np.arange(a[0],a[1]+h,h)
    else:
	t = np.linspace(a[0],a[1],h)
	
    n = len(t)
    h = t[1] - t[0]
    y = np.zeros_like(t)
    y[0] = y0
    i = 1
    while i < n:
	k = f(t[i-1],y[i-1])
	y[i] = y[i-1] + h*k
	i += 1
	
    return t,y

def beuler(f,y0,a,h):
    """ Calculates the Euler solution of the IVP
     y'=f(t,y), with y(a)=y0, a[0]<=t<=[1]. Using the implicit
     Backward Euler formula. If h<=1, then h is the step size
     else h is the number of nodes to use """



    if h <= 1:
	t = np.arange(a[0],a[1]+h,h)
    else:
	t = np.linspace(a[0],a[1],h)
	
    n = len(t)
    h = t[1] - t[0]
    y = np.zeros_like(t)
    y[0] = y0
    i = 1
    while i < n:
	f2 = lambda x: y[i-1]+h*f(t[i],x)-x
	#y[i] = fsolve(f2,y[i-1])
	y[i] = newton(f2,y[i-1])
	i += 1
	
    return t,y