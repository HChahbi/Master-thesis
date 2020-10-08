




#Black-Scholes Model ----------------


from scipy.stats import norm
import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd


class Option:

    def __init__(self,intrest_rate, vol,maturity,S0,K):

        self.intrest_rate=intrest_rate
        self.vol=vol
        self.maturity=maturity
        self.S=S0
        self.K=K


    def d_d(self):

        d_plus=(np.log(self.S/self.K)+(self.intrest_rate+0.5* self.vol**2)*self.maturity)/(self.vol*np.sqrt(self.maturity))
        d_minus=d_plus-self.vol*np.sqrt(self.maturity)
        return [d_plus,d_minus]


    def call_price(self):

        d=self.d_d()

        c_price=self.S*norm.pdf(d[0])-self.K*np.exp(-self.intrest_rate*self.maturity)*norm.pdf(d[1])

        return c_price

    def put_price_using_parity(self):

        c=self.call_price()

        p_price=c-self.S+self.K*np.exp(-self.intrest_rate*self.maturity)

        return p_price

    def put_price(self):
        d=self.d_d()
        p_price=K*np.exp(-self.intrest_rate*self.maturity)*(1-norm.pdf(d[1]))-self.S*(1-norm.pdf(d[0]))
        return p_price

        #butterfly spread---- is a combinaison of two short call with slightly diffferent srike prices and two long call with the same strike price


    def butterfly_price(self,epsilon):
        self.K-=epsilon
        c_k_epsilon_=self.call_price()
        self.K+=2*epsilon
        c_k_epsilon=self.call_price()
        self.K-=epsilon
        c_k=self.call_price()
        return c_k_epsilon_+c_k_epsilon-2*c_k
    
    
    def call_greeks(self):
        d=self.d_d()
        
        call=self.call_price()
        
        delta=norm.pdf(d[0])
        
        gamma=np.exp(-d[0]**2/2)/self.vol*self.S*np.sqrt(2*np.pi*self.maturity)
        
        rho=self.K*self.maturity*np.exp(-self.intrest_rate*self.maturity)*norm.pdf(d[1])
        vegga=self.S*np.sqrt(self.maturity/2*np.pi)*np.exp(-d[0]**2/2)
        theta=(self.vol**2*self.S**2/2)*gamma+self.intrest_rate*(self.S*delta-call)
        return [delta,gamma,rho,vegga,theta]    
    
    
    def put_greeks(self):
        d=self.d_d()
        put=self.put_price()
        delta=norm.pdf(d[0])-1
        gamma=np.exp(-d[0]**2/2)/self.vol*self.S*np.sqrt(2*np.pi*self.maturity)
        rho=self.K*self.maturity*np.exp(-self.intrest_rate*self.maturity)*norm.pdf(d[1])-self.K*self.maturity*np.exp(-self.intrest_rate*self.maturity)
        vegga=self.S*np.sqrt(self.maturity/2*np.pi)*np.exp(-d[0]**2/2)
        theta=(self.vol**2*self.S**2/2)*gamma+self.intrest_rate*(self.S*delta-put)
        return [delta,gamma,rho,vegga,theta]
    
    
    def butterfly_greeks(self,epsilon):
            self.K+=epsilon
            greeks_k_epsilon=np.array(self.call_greeks())
            self.K-=2*epsilon
            greeks_k_epsilon_=np.array(self.call_greeks())
            self.K+=epsilon
            greeks_k=np.array(self.call_greeks())
            greeks=greeks_k_epsilon+greeks_k_epsilon_-2*greeks_k
            return greeks
        
        
    def call_greeks_r(self,N_iterations,rmax):
        col=["r","Delta","Gamma","Vegga","Rho","Theta"]
        r=range(0,int(N_iterations*rmax),1)
        k=0
        A=np.zeros([int(N_iterations*rmax),6]) 
        for i in r:
            self.intrest_rate=i/N_iterations 
            g=[self.intrest_rate]+self.call_greeks()
            
            A[k,:]= np.array(g)
            k=k+1
        W=pd.DataFrame(A,columns=col)

        return W  
    
    
    def call_greeks_vol(self,N_iterations,volmax):    
        col=["vol","Delta","Gamma","Vegga","Rho","Theta"]
        r=range(0,int(N_iterations*volmax),1)
        k=0
        A=np.zeros([int(N_iterations*volmax),6]) 
        for i in r:
            self.vol=i/N_iterations 
            g=[self.vol]+self.call_greeks()
            A[k,:]= np.array(g)
            k=k+1
        W=pd.DataFrame(A,columns=col)

        return W  
    
    
    def call_greeks_T(self,N_iterations, Tmax):
        
        col=["T","Delta","Gamma","Vegga","Rho","Theta"]
        r=range(0,int(N_iterations*Tmax),1)
        k=0
        A=np.zeros([int(N_iterations*Tmax),6]) 
        for i in r:
            self.maturity=i/N_iterations
            g=[self.maturity]+self.call_greeks()      
            A[k,:]= np.array(g)
            k=k+1
        W=pd.DataFrame(A,columns=col)

        return W  
    
    def call_greeks_K(self,N_iterations, Kmax):
        col=["K","Delta","Gamma","Vegga","Rho","Theta"]
        r=range(0,int(N_iterations*Kmax),1)
        k=0
        A=np.zeros([int(N_iterations*Kmax),6]) 
        for i in r:
            self.K=i/N_iterations 
            g=[self.K]+self.call_greeks()
            
            A[k,:]= np.array(g)
            k=k+1
        W=pd.DataFrame(A,columns=col)

        return W  



# ---------------------Test------------------------


r=0.05
S=100
K=100
T=1/12
sigma=0.25


option=Option(r,sigma,T,S,K)

c_price=option.call_price()
p_price=option.put_price()

put_g=option.put_greeks()

call_g=option.call_greeks()


print("The call price is %f :\n",c_price)


print("The put price is %f: \n",p_price)

print("We present the greeks in the format: [delta,gamma,rho,vegga,theta]")

print("The greeks for the underlying call \n")
print(call_g)

print("The greeks for the underlying put\n")

print(put_g)

delta_k=0.02

print("Butterfly spread with a spread deltaK=%f \n", delta_k)


B_price=option.butterfly_price(delta_k)


print("The price of a butterfly spread is  %f\n", B_price)


B_greeks=option.butterfly_greeks(delta_k)

print("The greeks for butterfly spread \n")

print(B_greeks)





def greeks_r_plot(N_iterations,rmax):

            W=vanille.call_greeks_r(N_iterations,rmax)
            
            
            
            
            plt.figure()
            plt.plot(np.array(W["r"]), np.array(W['Delta']), label='Delta')
            plt.title("Delta VS r") 
            plt.show()
            
            plt.figure()
            plt.plot(np.array(W["r"]), np.array(W['Gamma']), label='Gamma')
            plt.title("Gamma VS r") 
            
            plt.show()
            
            plt.figure()
            plt.plot(np.array(W["r"]), np.array(W['Rho']), label='Rho')
            plt.title("Rho Vs r") 
            
            plt.show()
            
            plt.figure()
            plt.plot(np.array(W["r"]), np.array(W["Vegga"]), label='Vegga')
            plt.title("Vegga Vs r") 
            
            
            plt.show()
            
            plt.figure()
            plt.plot(np.array(W["r"]),np.array(W["Theta"]), label='Theta')
            plt.title("Theta Vs r") 
            
            plt.show()




def greeks_vol_plot(N_iterations,volmax):

            W=vanille.call_greeks_vol(N_iterations,volmax)
            
            
            
            plt.figure()
            plt.plot(np.array(W["vol"]), np.array(W['Delta']), label='Delta')
            plt.title("Delta VS vol") 

            plt.show()
            
            plt.figure()
            plt.plot(np.array(W["vol"]), np.array(W['Gamma']), label='Gamma')
            plt.title("Gamma VS vol") 

            plt.show()
            
            plt.figure()
            plt.plot(np.array(W["vol"]), np.array(W['Rho']), label='Rho')
            plt.title("Rho Vs vol") 

            plt.show()
            
            plt.figure()
            plt.plot(np.array(W["vol"]), np.array(W["Vegga"]), label='Vegga')
            plt.title("Vegga Vs r") 

            plt.show()
            
            plt.figure()
            plt.plot(np.array(W["vol"]),np.array(W["Theta"]), label='Theta')
            plt.title("Theta Vs r") 

            plt.show()


def greeks_T_plot(N_iterations,Tmax):

            W=vanille.call_greeks_T(N_iterations,Tmax)
            
            
            
            plt.figure()
            plt.plot(np.array(W["T"]), np.array(W['Delta']), label='Delta')
            plt.title("Delta VS maturity") 

            plt.show()
            
            plt.figure()
            plt.plot(np.array(W["T"]), np.array(W['Gamma']), label='Gamma')
            plt.title("Gamma VS maturity") 

            plt.show()
            
            plt.figure()
            plt.plot(np.array(W["T"]), np.array(W['Rho']), label='Rho')
            plt.title("Rho Vs maturity") 

            plt.show()
            
            plt.figure()
            plt.plot(np.array(W["T"]), np.array(W["Vegga"]), label='Vegga')
            plt.title("Vegga Vs maturity") 

            plt.show()
            
            plt.figure()
            plt.plot(np.array(W["T"]),np.array(W["Theta"]), label='Theta')
            plt.title("Theta Vs maturity") 

            plt.show()

def greeks_K_plot(N_iterations,Kmax):

            W=vanille.call_greeks_K(N_iterations,Kmax)
            
            
            
            plt.figure()
            plt.plot(np.array(W["K"]), np.array(W['Delta']), label='Delta')
            plt.title("Delta VS strike price") 

            plt.show()
            
            plt.figure()
            plt.plot(np.array(W["K"]), np.array(W['Gamma']), label='Gamma')
            plt.title("Gamma VS strike price") 

            plt.show()
            
            plt.figure()
            plt.plot(np.array(W["K"]), np.array(W['Rho']), label='Rho')
            plt.title("Rho Vs strike price") 

            plt.show()
            
            plt.figure()
            plt.plot(np.array(W["K"]), np.array(W["Vegga"]), label='Vegga')
            plt.title("Vegga Vs strike price") 

            plt.show()
            
            plt.figure()
            plt.plot(np.array(W["K"]),np.array(W["Theta"]), label='Theta')
            plt.title("Theta Vs strike price") 

            plt.show()






#--------------------------plotting------------------------------------
r=0.10
S=95
K=100
T=1
sigma=0.13
vanille=Option(r,sigma,T,S,K)
N_iterations=200




print("----------------call greeks plotting for different interest -----------")

rmax=0.2



greeks_r_plot(N_iterations,rmax)


print("----------------call greeks plotting for different volatility -----------")






volmax=0.40

greeks_vol_plot(N_iterations,volmax)


print("----------------call greeks plotting for different maturity-----------")


Tmax=10

greeks_T_plot(N_iterations,Tmax)


print("----------------call greeks plotting for different strike prices-----------")


Kmax=3*K

greeks_T_plot(N_iterations,Kmax)





