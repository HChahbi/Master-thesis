# Binomial lattice model


import numpy as np


class Binomial_Lattice:

    #we start by intanciate the binomial lattice and its parameters

    def __init__(self, interest_rate, maturity,N_T, volatility,S0,K):

        
        Delta_T=maturity/N_T

        a=(np.exp(-interest_rate*Delta_T)+np.exp((interest_rate+volatility**2)*Delta_T))/2


        d=a-np.sqrt(a**2-1)

        u=1/d

        p=(np.exp(interest_rate*Delta_T)-d)/(u-d)

        self.sigma=volatility
        self.r=interest_rate
        self.delta=Delta_T
        self.N_T=N_T
        self.d=d
        self.u=u
        self.p=p
        self.S=S0
        self.K=K

#payoffs of call and put

    def call_payoff(S,K):
        if (S-K)>0:
             return S-K
        return 0

    def put_payoff(S,K):
        return np.max((K-S),0)

# the price at the node n (time) m(outcome)

    def node_price(self,n,m):
        return self.S*(self.u**m)*(self.d**(n-m))


# we do the pricing backward, which means that we evaluate the option at the maturity for all the possible outcomes
#then we compute the value of the option for the previous time as the discounted risk free expected value

    def European_call_pricing(self):
        V_f=np.zeros(self.N_T+1)
        for m in range(self.N_T+1):
                p=self.node_price(self.N_T,m)
                V_f[m]=Binomial_Lattice.call_payoff(p,self.K)

        for m in range(self.N_T,0,-1):
            V=np.zeros(m)
            for n in range(m):
                V[n]=np.exp(-self.r*self.delta)*(self.p*V_f[n+1]+(1-self.p)*V_f[n])
                print(V[n])
            V_f=V

        return V_f[0]

    def European_put_pricing(self):
        V_f=np.zeros(self.N_T+1)
        for m in range(self.N_T+1):
                p=self.node_price(self.N_T,m)
                V_f[m]=Binomial_Lattice.put_payoff(p,self.K)

        for m in range(self.N_T,0,-1):
            V=np.zeros(m)
            for n in range(m):
                V[n]=np.exp(-self.r*self.delta)*(self.p*V_f[n+1]+(1-self.p)*V_f[n])
            V_f=V
        return V_f[0]
# the difference between a European and american vanille, is that the american one, u can exercise it whenever u want before the maturity

#the value of the option at the node n,m is the max of the normal payoff anf the discounted value

    def American_call_pricing(self):
        V_f=np.zeros(self.N_T+1)
        for m in range(self.N_T+1):
                V_f[m]=Binomial_Lattice.call_payoff(self.node_price(self.N_T,m),self.K)
                
        for m in range(self.N_T,0,-1):
            V=np.zeros(m)
            for n in range(m):   
                 payoff=Binomial_Lattice.call_payoff(self.node_price(m,n),self.K)
                 v=np.exp(-self.r*self.delta)*(self.p*V_f[n+1]+(1-self.p)*V_f[n])
                 if v>payoff:
                      V[n]=v
                 else:
                      V[n]=payoff
            V_f=V

        return V_f[0]

    def American_put_pricing(self):
        V_f=np.zeros(self.N_T+1)
        for m in range(self.N_T+1):

                p=self.node_price(self.N_T,m)
                V_f[m]=Binomial_Lattice.put_payoff(p,self.K)
                
        for m in range(self.N_T,0,-1):

            V=np.zeros(m)
            for n in range(m):
                 payoff=Binomial_Lattice.put_payoff(self.node_price(m,n),self.K)
                 v=np.exp(-self.r*self.delta)*(self.p*V_f[n+1]+(1-self.p)*V_f[n])
                 if v>payoff:
                      V[n]=v
                 else:
                      V[n]=payoff
            V_f=V

        return V_f[0]





#------------------Test---------


r=0.03
T=4/12
N_T=2
vol=0.30
S0=78
K=80

BL=Binomial_Lattice(0.03,4/12,2,0.30,78,80)

call_price=BL.European_call_pricing()
print("the call price is ")
print(call_price)

#amrican put option

A=Binomial_Lattice(-0.05+0.02,18/12,3,0.12,0.79,0.8)

american_call=A.American_put_pricing()

print("american call price is:")
print(american_call)