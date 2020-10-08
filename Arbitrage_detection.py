
# Single period arbitrage opportunity detection------


import numpy as np 
from numpy import random as rd
from scipy.stats import uniform
from scipy.stats import norm



# generate the price vector randomly --- uniform or normal laws



N=10 #number of securities
M=10 # number of possible outcomes- the cardinal of the state event space
R=0.1 # interest rate


def generate_price(vector_length=N):
	return np.array(uniform.rvs(loc=0, scale=1, size=vector_length)).T 

std=2
asset_return=[R]*N # expected return for different assets
asset_vol=[std]*N  # volatility of different assets


#fucntion to make the payoffs related to the asset i and the possible outcome j----it is a kind of a matrix D= Dij

def payoff_matrix(nb_outcomes=M, nb_securities=N, interest_rate=R, expected_return=asset_return,vol=asset_vol):
	D=np.zeros([nb_securities,nb_outcomes])
	D[0,:]=[interest_rate+1]*nb_outcomes

	for i in range(1,N):
		asset_i=norm.rvs(loc=expected_return[i],scale=vol[i],size=nb_outcomes)

		D[i,:]=asset_i

	return D


# we give to this function a payoff matrix and prices of the asset at the current time ---> return  psi vector of probabilities associated to each outcome


# no arbitrage opportunity <-> ] psi>0 p=D*psi

def arbitrage_opportunity(payoffs, price):
	psi=np.linalg.solve(payoffs,price)
	return psi

#to see if there is an arbitrage opportunity, you check only the probability vector if it contains negative element---> 
#which means that there is a problem with the market and we could find an arbitrage opportunity


def isthere_arbitrage(psi):
	for i in psi:
		if i<0:
			return True
	return False


# when we find a problem with the market -- negative probability-> means that there is an OA, that we develop for the state corresponding to the negative probability-
#we start with a negative portfolio value, and we have a non null probability to get a positive value from the market 
# D*strat=payoff  [0,0....,1,....0]-> strat

def opportunity_strategy(payoffs,psi, nb_securities=N):
    arbitrage_payoff=[0]*nb_securities
    for i in range(N):
        if psi[i]<0:
            arbitrage_payoff[i]=1
            print("state i where there is an arbitrage opportunnity")
            print(i)
        break
    print("payoffs induced by doing the following strategy")
    print(arbitrage_payoff)		
    strategy=np.linalg.solve(D.T,arbitrage_payoff)
    return strategy


# we verify that the strategy induce effectively a gain even with a negative starting protfolio value.

def verify_OA_strategy(strategy ,payoffs,price):
    boolean=np.dot(payoffs.T,strategy)>=0
    if np.dot(price,strategy)<=0 :
        for j in boolean:
            if j==True:
                return True
    return False




#--------Test-----------------------------------
    
p=generate_price()

D=payoff_matrix()


psi=arbitrage_opportunity(D,p)

if isthere_arbitrage(psi):
    strat=opportunity_strategy(D,psi, nb_securities=N)
    print("Strategy----------")
    print(strat)

print("Probaility vector------- notice that sum up of its coefficient is %f \n ",np.sum(psi))

print("if there is a matter about probabilities or a negative probability shows up, then there is mostly an opotunity of arbitrage")


print(psi)



print("Payoffs-------------")
    
print(np.dot(D.T,strat))

print("Initial value of the portfolio------------------")
print(np.dot(p,strat))

print("is is true that strategy induce an opportunity of arbitrage-----")
print(verify_OA_strategy(strat,D,p))




