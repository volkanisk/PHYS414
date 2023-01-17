from function_newton import *

data = read_csv()       # reading the data from csv

log_m = np.log(data[:,1])       # getting log of mass
log_r = -(log_m - data[:,0])/2  # transforming surface gravity to log R

plt.scatter(log_r,log_m)        # log_r log_m graph for all r and m
plt.xlabel("log_R")
plt.ylabel("mass")


log_data = np.transpose(np.vstack((log_m,log_r)))
small_log_data = log_data[log_data[:,0] < -0.5,:]   # I only used stars that have log_mass < -0.5, I got 0.5 from first graph

plt.scatter(small_log_data[:,1],small_log_data[:,0])    # Plotting for only small masses
plt.xlabel("log_R")
plt.ylabel("mass")
plt.show()


result = stats.linregress(small_log_data[:,1],small_log_data[:,0])  # linear regression on small mass data

n = (result.slope -3) / (result.slope - 1)      # the slope and n value has this relation
q = round(5*n/(1+n))        # from n value I get q, and round it to get integer value
print("q = {}".format(q))   # I get 3 as q value

### IVP part
guess_1 = 3; guess_2 = 4    # initial Xi guesses, which I got experimentally
i =0
while (guess_1 != guess_2) and (i <20):
    guess_1, guess_2, func_at_r, der_at_r = shooting(guess_1, guess_2)      # You can read function detail in function file
    i+=1

Xi = guess_1    #  Getting the Xi of star


M = data[:,1]
R = np.exp(log_r)

rho_c = - M / (4*np.pi) / np.power(R,3) /der_at_r * Xi      # this equation comes from eqn. (6)

plt.scatter(M,rho_c)
plt.ylabel("rho_c")
plt.xlabel("mass")
plt.show()

## Calculating k*
Nn = (4* np.pi)**(2/3) /(5/2) * ( -Xi**2 * der_at_r)**(-1/3) * Xi
G = 6.6743 * 10**(-11)
K = np.matmul(R, np.power(M, (1/3))) * G * Nn
print("K* = {}".format(K))