# APPM 4600 - HW3

#####################################################################
# Problem 2


#Ti is initial soil temp before snap
Ti = 20
#Ts is constant temp during cold period
Ts = -15
#alpha is thermal conductivity
alpha = 0.138e-6
#60 day converted into seconds
t = 518400


import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erf

def f(x):
    return erf(x/(2*np.sqrt(alpha*t))) - 3/7

def fp(x):
   return (1/np.sqrt(np.pi*alpha*t))*np.exp(-x**2/(4*alpha*t))

print(f(1))
print(f(.5))

x = np.linspace(0, 0.5, 1000)

plt.plot(x,f(x))
plt.axhline(0)
plt.xlabel("x (meters)")
plt.ylabel("f(x)")
plt.show()

# def driver():

# # use routines    
    
#     a = 0
#     b = 0.5

# #    f = lambda x: np.sin(x)
# #    a = 0.1
# #    b = np.pi+0.1

#     tol = 1e-13

#     [astar,ier] = bisection(f,a,b,tol)
#     print('the approximate root is',astar)
#     print('the error message reads:',ier)
#     print('f(astar) =', f(astar))


# # define routines
# def bisection(f,a,b,tol):
    
#     fa = f(a)
#     fb = f(b);
#     if (fa*fb>0):
#        ier = 1
#        astar = a
#        return [astar, ier]

# #   verify end points are not a root 
#     if (fa == 0):
#       astar = a
#       ier =0
#       return [astar, ier]

#     if (fb ==0):
#       astar = b
#       ier = 0
#       return [astar, ier]

#     count = 0
#     d = 0.5*(a+b)
#     while (abs(d-a)> tol):
#       fd = f(d)
#       if (fd ==0):
#         astar = d
#         ier = 0
#         return [astar, ier]
#       if (fa*fd<0):
#          b = d
#       else: 
#         a = d
#         fa = fd
#       d = 0.5*(a+b)
#       count = count +1
# #      print('abs(d-a) = ', abs(d-a))
      
#     astar = d
#     ier = 0
#     print('count = ', count)
#     return [astar, ier]
#driver()



        
def driver():
#f = lambda x: (x-2)**3
#fp = lambda x: 3*(x-2)**2
#p0 = 1.2

   
  p0 = 0.5

  Nmax = 100
  tol = 1.e-13

  (p,pstar,info,it) = newton(f,fp,p0,tol, Nmax)
  print('the approximate root is', '%16.16e' % pstar)
  print('the error message reads:', '%d' % info)
  print('Number of iterations:', '%d' % it)


def newton(f,fp,p0,tol,Nmax):
  """
  Newton iteration.
  
  Inputs:
    f,fp - function and derivative
    p0   - initial guess for root
    tol  - iteration stops when p_n,p_{n+1} are within tol
    Nmax - max number of iterations
  Returns:
    p     - an array of the iterates
    pstar - the last iterate
    info  - success message
          - 0 if we met tol
          - 1 if we hit Nmax iterations (fail)
     
  """
  p = np.zeros(Nmax+1);
  p[0] = p0
  for it in range(Nmax):
      p1 = p0-f(p0)/fp(p0)
      p[it+1] = p1
      if (abs(p1-p0) < tol):
          pstar = p1
          info = 0
          return [p,pstar,info,it]
      p0 = p1
  pstar = p1
  info = 1
  return [p,pstar,info,it]
        
driver()


#####################################################################

# Problem 3

import numpy as np
    
def driver1():

# test functions 
    f1 = lambda x1: x1*(1+ (7 - x1**5)/x1**2)**3

    f2 = lambda x1: x1 - (x1**5 - 7)/x1**2

    f3 = lambda x1: x1 - (x1**5 - 7)/(5*x1**4)

    f4 = lambda x1: x1 - (x1**5 - 7)/12



    Nmax = 100
    tol = 1e-10

# test f1 '''
    try:
        x0 = 1
        [xstar,ier] = fixedpt(f1,x0,tol,Nmax)
        print('the approximate fixed point is:',xstar)
        print('f1(xstar):',f1(xstar))
        print('Error message reads:',ier)
    except OverflowError:
        print('f1 diverged: overflow occurred')     
    
#test f2 '''
    try:
        x0 = 1
        [xstar,ier] = fixedpt(f2,x0,tol,Nmax)
        print('the approximate fixed point is:',xstar)
        print('f2(xstar):',f2(xstar))
        print('Error message reads:',ier)
    except OverflowError:
            print('f2 diverged: overflow occurred') 
# test f3 '''
    x0 = 1
    [xstar,ier] = fixedpt(f3,x0,tol,Nmax)
    print('the approximate fixed point is:',xstar)
    print('f3(xstar):',f3(xstar))
    print('Error message reads:',ier)

# test f4 '''
    x0 = 1
    [xstar,ier] = fixedpt(f4,x0,tol,Nmax)
    print('the approximate fixed point is:',xstar)
    print('f4(xstar):',f4(xstar))
    print('Error message reads:',ier)



# define routines
def fixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''

    count = 0
    while (count <Nmax):
       count = count +1
       x1 = f(x0)
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          return [xstar,ier]
       x0 = x1

    xstar = x1
    ier = 1
    return [xstar, ier]
    

driver1()





