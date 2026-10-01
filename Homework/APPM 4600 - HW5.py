#############
#APPM 4600 Homework 5
#############

import numpy as np
import math
import time
from numpy.linalg import inv 
from numpy.linalg import norm 

def driver():

    x0 = np.array([1,1])
    
    Nmax = 100
    tol = 1e-10
    [xstar,ier,its] = fixedpoint(x0,tol,Nmax)

    print(xstar)
    print('Fixed Point: the error message reads:',ier)
    print('Fixed Point: number of iterations is:',its)


    t = time.time()
    for j in range(50):
      [xstar,ier,its] =  Newton(x0,tol,Nmax)
    elapsed = time.time()-t
    print(xstar)
    print('Newton: the error message reads:',ier) 
    print('Newton: took this many seconds:',elapsed/50)
    print('Netwon: number of iterations is:',its)
     
   
def evalF(x): 
# vector function that you want to find the roots of

    F = np.zeros(2)
    
    F[0] = 3*x[0]**2 - x[1]**2
    F[1] = 3*x[0]*x[1]**2 - x[0]**3 - 1
    
 
    return F

def fixedpoint(x0, tol, Nmax):

    B = np.array([[1/6, 1/18],
                  [0,   1/6]])

    for its in range(Nmax):

        F = evalF(x0)

        x1 = x0 - B.dot(F)

        if norm(x1-x0) < tol:
            xstar = x1
            ier = 0
            return [xstar, ier, its]

        x0 = x1

    xstar = x1
    ier = 1
    return [xstar, ier, its]



def evalJ(x): 
# Jacobian of the vector function you want to find the roots of
    
    J = np.array([
        [6*x[0], -2*x[1]],
        [3*x[1]**2 - 3*x[0]**2, 6*x[0]*x[1]]
    ])
    
    return J


def Newton(x0,tol,Nmax):

    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''

    for its in range(Nmax):
       J = evalJ(x0)
       Jinv = inv(J)
       F = evalF(x0)
       
       x1 = x0 - Jinv.dot(F)
       
       if (norm(x1-x0) < tol):
           xstar = x1
           ier =0
           return[xstar, ier, its]
           
       x0 = x1
    
    xstar = x1
    ier = 1
    return[xstar,ier,its]
           

        
if __name__ == '__main__':
    # run the drivers only if this is called from the command line
    driver()       
