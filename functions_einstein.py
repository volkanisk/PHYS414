import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


def func(r,y):
    "This function is the vectorized form of differential equations in the question"
    P = y[0]
    rho = np.sqrt(P / 100)
    M = y[1]
    value_p = -M + (4 * np.pi * r**3 * P) / (r * (r - 2*M)) * (rho + P)
    value_m = 4 * np.pi * r**2 * rho
    return [value_p, value_m]


def shooting(guess_1, guess_2, rho_c):
    "This function is a shooting method imlemented on the func function."
    eps = 0.000001
    guess = (guess_1 + guess_2)/2
    x = np.linspace(0.1, guess, 10)
    sol = solve_ivp(fun=func, t_span=[0.1, guess], y0=(rho_c, 0), t_eval=x)

    if (len(sol.y[0]) != 10) :          # I observed that if R is bigger than R of star, elements size decreases
        guess_2 = (guess_1 + guess_2)/2
    elif (len(sol.y[0]) == 10) and (abs(sol.y[0,-1]) > eps) :   # I want the P function to be near to zero
        guess_1 = (guess_1 + guess_2)/2
    elif (len(sol.y[0]) == 10) and (abs(sol.y[0,-1]) < eps):    # If P is near zero, function returns guess_1,2 as equal
        guess_1 = guess
        guess_2 = guess
    return guess_1,guess_2,sol.y[0,-1], sol.y[1,-1]

def M_R(rho_c):
    # Initial guesses, I got those experimentally
    guess_1 = 0.5;
    guess_2 = 1.5
    i = 0
    while (guess_1 != guess_2) and (i < 20):    # Unless guesses are equal keep on getting integrals
        guess_1, guess_2, func_at_r, m_at_r = shooting(guess_1, guess_2, rho_c)
        i += 1
    Mo = 1.989 * 10 ** 30
    M = m_at_r * Mo
    R = guess_1 * 1477
    return M, R     # function returns M, R as kg and m
