
# import libraries
import numpy as np
def driver():
# test functions
    f1 = lambda x: (10/(x+4))**(1/2)
# fixed point is alpha1 = 1.4987....
    #f2 = lambda x: 3+2*np.sin(x)
#fixed point is alpha2 = 3.09...
    Nmax = 100
    tol = 1e-11
# test f1 '''
    x0 = 1.5
    [xstar,ier,guesses] = fixedpt(f1,x0,tol,Nmax)
    print('the approximate fixed point is:',xstar)
    print('f1(xstar):',f1(xstar))
    print('Error message reads:',ier)
    print(order_of_convergence(guesses,1.3652300134140976))
#test f2 '''
    # x0 = 0.0
    # [xstar,ier] = fixedpt(f2,x0,tol,Nmax)
    # print('the approximate fixed point is:',xstar)
    # print('f2(xstar):',f2(xstar))   
    # print('Error message reads:',ier)
# define routines
def fixedpt(f,x0,tol,Nmax):
# ''' x0 = initial guess'''
# ''' Nmax = max number of iterations'''
# ''' tol = stopping tolerance'''
    count = 0

    guesses = []
    while (count <Nmax):
        count = count +1
        x1 = f(x0)
        guesses.append(x1)
        if (abs(x1-x0) <tol):
            xstar = x1
            ier = 0
            return [xstar,ier,np.array(guesses)]
    x0 = x1

    xstar = x1
    ier = 1
    return [xstar,ier,np.array(guesses)]


def order_of_convergence(guesses,p):
    
    error_np1 = np.log(abs(guesses[-1]-p))
    error_n = np.log(abs(guesses[-2]-p))
    error_nm1 = np.log(abs(guesses[-3]-p))
    alpha = (error_np1/error_n)/(error_n/error_nm1)
    return alpha
driver()