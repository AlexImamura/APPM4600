import matplotlib.pyplot as plt


columns = ["Method", "Input", "Iteration", "Idea Behind Method", "Required for Convergence", "Pros", "Cons"]
rows = [["Bisection", "", "", "", "", "", ""], 
        ["Fixed Point", "", "", "", "", "", ""],
        ["Newton", "", "", "", "", "", ""],
        ["Secant", "", "", "", "", "", ""]]

fig, ax = plt.subplots(figsize=(12,6))

ax.set_axis_off()

table = ax.table(
    cellText=rows,
    colLabels=columns,
    cellLoc="center",
    loc="center",
)


# Increase row height
for row in range(5):
    for col in range(7):
        table[(row, col)]

plt.title("Chart for root finding", fontsize=14)

plt.show()

import numpy as np

def driver():
    f = lambda x: np.exp(x**2+7*x-30) - 1
    fp = lambda x: np.exp(x**2 + 7*x - 30)*(2*x+7)
    fpp = lambda x: np.exp(x**2 +7*x -30)*((2*x+7)**2 + 2)
    g = lambda x: x - f(x)/fp(x)
    gp = lambda x: 1 - (fp(x)**2 - f(x)*fpp(x))/(fp(x)**2)
    

    Nmax = 100
    tol = 1e-14
    (p0,ier) = bisection(f,a,b,gp)
    (p,pstar,info, it) = newton(f,fp,p0,tol,Nmax)
    print('the approximate root is', '%16.16e' % pstar)
    print('the error message reads:', '%d' % info)
    print('Number of iterations:', '%d' % it)

###############################################################################
[a,b] = [2,4.5]
 
def bisection(f,a,b,gp):
    
#    Inputs:
#     f,a,b       - function and endpoints of initial interval
#      tol  - bisection stops when interval length < tol

#    Returns:
#      astar - approximation of root
#      ier   - error message
#            - ier = 1 => Failed
#            - ier = 0 == success

#     first verify there is a root we can find in the interval 

    fa = f(a)
    fb = f(b);
    if (fa*fb>0):
       ier = 1
       astar = a
       return [astar, ier]

#   verify end points are not a root 
    if (fa == 0):
      astar = a
      ier =0
      return [astar, ier]

    if (fb ==0):
      astar = b
      ier = 0
      return [astar, ier]

    count = 0
    d = 0.5*(a+b)
    while (abs(gp(d)) >= 1):
      fd = f(d)
      if (fd ==0):
        astar = d
        ier = 0
        return [astar, ier]
      if (fa*fd<0):
         b = d
      else: 
        a = d
        fa = fd
      d = 0.5*(a+b)
      count = count +1
#      print('abs(d-a) = ', abs(d-a))
      
    astar = d
    ier = 0
    print('count = ', count)
    return [astar, ier]
      
              


###############################################################################

def newton(f,fp,p0,tol,Nmax):
    p = np.zeros(Nmax+1);
    p[0] = p0
    for it in range(Nmax):
        p1 = p0-f(p0)/fp(p0)
        p[it+1] = p1
        if (abs(p1-p0)<tol):
            pstar = p1
            info = 0
            return [p,pstar,info,it]
        p0 = p1
    pstar = p1
    info = 1
    return [p,pstar,info,it]
driver()



