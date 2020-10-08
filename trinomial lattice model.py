#trinomial lattice model

import numpy as np


class trinomial_Lattice:
     # for trinomial tree, there are essentially two equations that allows to find the triple probabilities: expected return and variance equation
     # deltax measures the variance of the price for a given time 

        def __init__(self,r,vol,mu,N_T,deltax,T,S0,K):
             
            self.r=r
            self.vol=vol
            self.mu=mu
            self.N_T=N_T
            self.deltax=deltax
            self.delta=T/N_T
            self.S=S0
            self.K=K

            #computation of the lattice probabilities
            if self.mu==0:
                 #driftless lattice
                self.p_u=((np.exp(-2*self.deltax)-1)*(np.exp(self.r*self.delta)-1)-(np.exp(-self.deltax)-1)*(np.exp((2*self.r+self.vol**2)*self.delta)-1))/((np.exp(self.deltax)-1)*(np.exp(-2*self.deltax)-1)-(np.exp(self.deltax)-1)*(np.exp(2*self.deltax)-1))

                self.p_d=((np.exp(self.r*self.delta)-1)-(np.exp(self.deltax)-1))*self.p_u/((np.exp(-self.deltax)-1))

                self.p_m=1-self.p_d-self.p_u
                
            else:
                
                #drifted lattice with symetrical tree p+=p-                
                p=self.vol**2*self.delta/2*self.deltax**2
                # p<1, is a kind of the stability of the pde of the black sholes equation.
                self.p_d=p
                self.p_u=p
                self.p_m=1-self.p_u-self.p_u
                self.mu=self.r-1/self.delta*np.log(2*p*(np.cosh(self.deltax)-1)+1)
        
        

#payoffs of call and put

        def call_payoff(S,K):
            if S>K:
                 return S-K
            return 0

        def put_payoff(S,K):
            if K>S:
                 return K-S
            return 0

# the price at the node m(time) n(outcome)

        def node_price(self,m,n):
            return self.S*np.exp(m*self.mu*self.delta+n*self.deltax)

#for each m delta, they will be 2m+1 delta x---- n goes from -m to m
         
# we follow a backward pricing procedure, in other words, we start from the end and we use the probability transfer matrix to get back and find the value of the option
         #at the time 0
         
        def European_call_pricing(self):

            V_f=np.zeros(2*self.N_T+1)

            for n in range(-self.N_T,self.N_T+1):

                    p=self.node_price(self.N_T,n)

                    V_f[n+self.N_T]=trinomial_Lattice.call_payoff(p,self.K)
          
            
            NT=2*self.N_T+1
            T=np.zeros([NT,NT])

            T[0,0]=self.p_m
            T[0,1]=self.p_u

            for i in range(1,NT-1):
                T[i,i-1]=self.p_d
                T[i,i]=self.p_m
                T[i,i+1]=self.p_u

            T[NT-1,NT-2]=self.p_d
            T[NT-1,NT-1]=self.p_m
            for m in range(self.N_T-1,-1,-1):

                V_f=np.exp(-self.r*self.delta)*np.dot(T,V_f)
                
            print(V_f)

            return V_f[0]

        def European_put_pricing(self):

            V_f=np.zeros(2*self.N_T+1)

            for n in range(-self.N_T,self.N_T+1):

                    p=self.node_price(self.N_T,n)

                    V_f[n+self.N_T]=trinomial_Lattice.put_payoff(p,self.K)
            
            NT=2*self.N_T+1
            T=np.zeros([NT,NT])

            T[0,0]=self.p_m
            T[0,1]=self.p_u

            for i in range(1,NT-1):
                 
                T[i,i-1]=self.p_d
                T[i,i]=self.p_m
                T[i,i+1]=self.p_u
                
            T[NT-1,NT-2]=self.p_d
            T[NT-1,NT-1]=self.p_m

            for m in range(self.N_T-1,-1,-1):

                V_f=np.exp(-self.r*self.delta)*np.dot(T,V_f)
                print(V_f)

            return V_f[0]
       
        def American_put_pricing(self):

            V_f=np.zeros(2*self.N_T+1)

            for n in range(-self.N_T,self.N_T+1):

                    p=self.node_price(self.N_T,n)

                    V_f[n+self.N_T]=trinomial_Lattice.put_payoff(p,self.K)
            
            NT=2*self.N_T+1
            T=np.zeros([NT,NT])

            T[0,0]=self.p_m
            T[0,1]=self.p_u

            for i in range(1,NT-1):
                T[i,i-1]=self.p_d
                T[i,i]=self.p_m
                T[i,i+1]=self.p_u

            T[NT-1,NT-2]=self.p_d
            T[NT-1,NT-1]=self.p_m

            for m in range(self.N_T-1,-1,-1):
                 
                V=np.exp(-self.r*self.delta)*np.dot(T,V_f)
                
                for n in range(2*self.N_T+1):
                     
                    S=self.node_price(m,n)
                    v=trinomial_Lattice.call_payoff(S,self.K)

                    if v>V[n]:
                         
                         V_f[n]=v
                    else:
                         V_f[n]=V[n]
            return V_f[0]
                





        def American_call_pricing(self):
            
            V_f=np.zeros((2*self.N_T+1))
            
            print('dsfsdfs')

            for n in range(-self.N_T,self.N_T+1):

                    p=self.node_price(self.N_T,n)

                    V_f[n+self.N_T]=trinomial_Lattice.call_payoff(p,self.K)
            
            NT=2*self.N_T+1
            T=np.zeros([NT,NT])

            T[0,0]=self.p_m
            T[0,1]=self.p_u

            for i in range(1,NT-1):
                T[i,i-1]=self.p_d
                T[i,i]=self.p_m
                T[i,i+1]=self.p_u

            T[NT-1,NT-2]=self.p_d
            T[NT-1,NT-1]=self.p_m

            for m in range(self.N_T-1,-1,-1):
                 

                V=np.exp(-self.r*self.delta)*np.dot(T,V_f)
                for n in range(2*self.N_T+1):
                    S=self.node_price(m,n)
                    
                    v=trinomial_Lattice.call_payoff(S,self.K)

                    if v>V[n]:
                         
                         V_f[n]=v
                    else:
                         V_f[n]=V[n]
                         
            return V_f[0]


            
#------------------Test----------------------------


r=0.03
T=4/12
N_T=2
vol=0.30
S0=78
K=80
mu=0.05
deltax=0.5

TL=trinomial_Lattice(r,vol,mu,N_T,deltax,T,S0,K)

call_price=TL.European_call_pricing()
print("the call price is ")
print(call_price)

#amrican put option

A=trinomial_Lattice(r,vol,mu,N_T,deltax,T,S0,K)

american_put=A.American_put_pricing()

print("american call price is:")
print(american_put)