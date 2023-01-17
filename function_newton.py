import numpy as np
from read_csv import read_csv
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.integrate import solve_ivp


def read_csv():
    "This function gets the data of column 1-2 skipping the first line"
    data_set = np.genfromtxt("white_dwarf_data.csv", usecols= (1,2) , delimiter=",", skip_header=1  )
    return data_set

def func(xi,y):
    "This function is vectorized form of differential equation"
    Q = y[0]
    z = y[1]
    value = -(Q**1.5 + 2*z/xi)      # 1.5 is the value of n for integer q
    return [z, value]

def shooting(guess_1, guess_2):
    "This function uses shooting method from the function defined above."
    eps = 0.000001
    guess = (guess_1 + guess_2)/2
    x = np.linspace(0.1, guess, 10)
    sol = solve_ivp(fun=func, t_span=[0.1, guess], y0=(1, 0), t_eval=x)
    if (len(sol.y[0]) != 10) :              # I realized that after the R of the star, length of values decreases
        guess_2 = (guess_1 + guess_2)/2
    elif (len(sol.y[0]) == 10) and (abs(sol.y[0,-1]) > eps) :
        guess_1 = (guess_1 + guess_2)/2
    elif (len(sol.y[0]) == 10) and (abs(sol.y[0,-1]) < eps):
        guess_1 = guess
        guess_2 = guess

    return guess_1,guess_2,sol.y[0,-1], sol.y[1,-1]