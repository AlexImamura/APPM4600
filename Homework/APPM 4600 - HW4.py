import numpy as np
import matplotlib.pyplot as plt


def driver():

    f = lambda x: x**6 - x - 1
    fp = lambda x: 6*x**5 - 1

    p0 = 2
    p1 = 1

    Nmax = 100
    tol = 1.e-14

    # Find the exact root
    roots = np.roots([1, 0, 0, 0, 0, -1, -1])
    real_roots = [r.real for r in roots if abs(r.imag) < 1.e-10]
    alpha = max(real_roots)

    print("alpha =", alpha)

    # ----------------
    # Newton Method
    # ----------------

    (p, pstar, info, it) = newton(f, fp, p0, tol, Nmax)

    print("\nNewton Method")
    print("-----------------------------------------------")
    print(f"{'k':<5} {'x_k':<20} {'|x_k - alpha|':<20}")
    print("-----------------------------------------------")

    newton_errors = np.abs(p - alpha)

    for k in range(it + 2):
        error = abs(p[k] - alpha)
        print(f"{k:<5} {p[k]:<20.12f} {error:<20.8e}")

    # ----------------
    # Secant Method
    # ----------------

    secant_guesses = secant(f, p0, p1, Nmax, tol)

    secant_errors = np.array([abs(x - alpha) for x in secant_guesses])

    print("\nSecant Method")
    print("-----------------------------------------------")
    print(f"{'k':<5} {'x_k':<20} {'|x_k - alpha|':<20}")
    print("-----------------------------------------------")

    for k in range(len(secant_guesses)):
        error = abs(secant_guesses[k] - alpha)
        print(f"{k:<5} {secant_guesses[k]:<20.12f} {error:<20.8e}")


    # ----------------
    # Part (b)
    # ----------------

    # Newton
    newton_x = newton_errors[:-1]
    newton_y = newton_errors[1:]

    # Secant
    secant_x = secant_errors[:-1]
    secant_y = secant_errors[1:]

    # Plot
    plt.figure()

    plt.loglog(newton_x, newton_y, 'o-', label='Newton')
    plt.loglog(secant_x, secant_y, 's-', label='Secant')

    plt.xlabel(r'$|x_k-\alpha|$')
    plt.ylabel(r'$|x_{k+1}-\alpha|$')
    plt.title('Order of Convergence')
    plt.legend()
    plt.grid(True)

    plt.show()


def newton(f, fp, p0, tol, Nmax):

    p = np.zeros(Nmax + 1)
    p[0] = p0

    for it in range(Nmax):

        p1 = p0 - f(p0) / fp(p0)
        p[it + 1] = p1

        if abs(p1 - p0) < tol:
            pstar = p1
            info = 0
            return [p, pstar, info, it]

        p0 = p1

    pstar = p1
    info = 1

    return [p, pstar, info, it]


def secant(f, p0, p1, Nmax, tol):

    guesses = [p0, p1]

    for k in range(Nmax):

        pnew = p1 - f(p1) * (p1 - p0) / (f(p1) - f(p0))

        guesses.append(pnew)

        if abs(pnew - p1) < tol:
            break

        p0 = p1
        p1 = pnew

    return guesses


driver()