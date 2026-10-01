"""
 This compares Newton's method, Broyden and Lazy Newton for 
 computing the roots of vector valued functions.
 The function and the Jacobian are stored in subroutines and need 
 to be changed for different problems.  
 
"""

############################################# 
"""
Copyright (C) 2025  Adrianna M. Gillman

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
############################################# 


import numpy as np
import math
import time
from numpy.linalg import inv 
from numpy.linalg import norm 

def driver():

    x0 = np.array([1, 0])
    
    Nmax = 100
    tol = 1e-10
    recompute_tol = 1e-5

    
    t = time.time()
    for j in range(50):
      [xstar,ier,its] =  Newton(x0,tol,Nmax)
    elapsed = time.time()-t
    print(xstar)
    print('Newton: the error message reads:',ier) 
    print('Newton: took this many seconds:',elapsed/50)
    print('Netwon: number of iterations is:',its)
     
    t = time.time()
    for j in range(20):
      [xstar,ier,its] =  LazyNewton(x0,tol,Nmax)
    elapsed = time.time()-t
    print(xstar)
    print('Lazy Newton: the error message reads:',ier)
    print('Lazy Newton: took this many seconds:',elapsed/20)
    print('Lazy Newton: number of iterations is:',its)
     
    
    t = time.time()
    for j in range(20):
      [xstar,ier,its] =  SlackerNewton(x0,tol,Nmax,recompute_tol)
    elapsed = time.time()-t
    print(xstar)
    print('Slacker Newton: the error message reads:',ier)
    print('Slacker Newton: took this many seconds:',elapsed/20)
    print('Slacker Newton: number of iterations is:',its)


def evalF(x): 
# vector function that you want to find the roots of

    F = np.zeros(2)
    
    F[0] = 4*x[0]**2 + x[1]**2 - 4
    F[1] = x[0] + x[1] - np.sin(x[0]-x[1])   
 
    return F
    
def evalJ(x): 
# Jacobian of the vector function you want to find the roots of
    
    J = np.array([[8*x[0], 2*x[1]], 
        [1 - np.cos(x[0] - x[1]), 1 + np.cos(x[0] - x[1])]
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
           
def LazyNewton(x0,tol,Nmax):

    ''' Lazy Newton = use only the inverse of the Jacobian for initial guess'''
    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''

    J = evalJ(x0)
    Jinv = inv(J)
    for its in range(Nmax):

       F = evalF(x0)
       x1 = x0 - Jinv.dot(F)
       
       if (norm(x1-x0) < tol):
           xstar = x1
           ier =0
           return[xstar, ier,its]
           
       x0 = x1
    
    xstar = x1
    ier = 1
    return[xstar,ier,its]   
    
def SlackerNewton(x0,tol,Nmax,recompute_tol):

    J = evalJ(x0)
    Jinv = inv(J)
    prev_update = math.inf

    for its in range(Nmax):
    
        F = evalF(x0)
        update = Jinv.dot(F)
        

        if (abs(norm(update) - norm(prev_update)) < recompute_tol):
            J = evalJ(x1)
            Jinv = inv(J)
            update = Jinv.dot(F)


        x1 = x0 - update
        prev_update = update


        if (norm(x1-x0) < tol):
            xstar = x1
            ier =0
            return[xstar, ier,its]
               
        x0 = x1
        
    xstar = x1
    ier = 1
    return[xstar,ier,its]   
     
   
if __name__ == '__main__':
    # run the drivers only if this is called from the command line
    driver()      


######
# Pre lab
######
# For test with x0 = (2,0.5), Newton converges to (1,1) and Lazy Newton diverges.
# For x0 = (3,5), Newton diverges while Lazy newton converges to (-0.47767006  1.33110428)

#Exercise 3.2:
# Question 4: Our Slacker Newton codes performed very similar, we both got them to converge to nearly the same numbers, 
# however it seems like his code converges to 1 from above while mine converges from below. Additionally mine converged
# in 5 iterations while his took 8 iterations.

#Question 5: My code converges in 5 iterations while lazy newton converges in 7, however the time it takes to converges for Slacker
# to converge is slightly slower than lazy newton.


