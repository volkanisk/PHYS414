from functions_einstein import *

# getting M - R values for a rho_c
rho_c = 10**(-3)
guess_1 = 0.5; guess_2 = 0.8
i =0
while (guess_1 != guess_2) and (i <20):
    guess_1, guess_2, func_at_r, m_at_r = shooting(guess_1, guess_2, rho_c)
    i+=1
Mo = 1.989 * 10**30
M = m_at_r * Mo
R = guess_1 * 1477

# getting M- R graph for list of rho_c
rho_c_list = np.linspace(10**(-3), 10**(-2), 100)
M_list = np.zeros(100)
R_list = np.zeros(100)
i = 0
for element in rho_c_list:
    M_list[i], R_list[i] = M_R(element)
    i += 1
plt.scatter(R_list,M_list)
plt.xlabel("R (m)")
plt.ylabel("M (kg)")
plt.show()