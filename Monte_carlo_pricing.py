# Monte carlo for basket princing 

    

import numpy as np 
import pandas as pd 



# =============================================================================
# find the implied volatility using the dichotomy method- we seek the volatility in the BS model that matches the market price 
# =============================================================================
def implied_volatility(market_price,S, K, r,T):

    vol_l=0.10
    vol_h=0.45
    vol=(vol_h+vol_l)/2
    c_price=BS.call_price(S,K,vol,r,T)

    while abs(c_price-market_price)<0.1:

        if c_price<market_price:
            vol_l=vol
        else:
            vol_h=vol

        vol=(vol_h+vol_l)/2
        c_price=BS.call_price(S,K,vol,r,T)

    return vol


# =============================================================================
# transform the the price data into a centered, quite normalized return data--- we suppose that we are working under the black scholes framework
# =============================================================================

def asset_returns(data,intrest_rate, asset_vol):

    N_time,N_assets=data.shape()
    Return_data=np.array((N_time-1,N_assets))
    for i in range(N_assets):
        for k in range(N_time-1):
            centring=intrest_rate-1/2*asset_vol[i]**2
            Return_data[k][i]=(np.log(data[k+1]/data[k])-centring)/np.sqrt(1/N_time)

    return pd.DataFrame(Return_data)



# =============================================================================
# using a real data to get the correlation between the different assets of our basket
#this function returns the matrix of correlation between different assets..
# =============================================================================

def corr(data,asset_vol):
    #we suppose that data: each column corresponds to the return of an asset-- 
    #rows correspond to the time series of the assets 

    N_time,N_assets=data.shape()

    Correlation=np.array((N_assets,N_assets))

    for i in range(N_assets):
        for j in range(i,N_assets):
            sigma=0
            for k in range(N_time):

                sigma+=data[k][i]*data[k][j]

            Correlation[i][j]=sigma/N_time *(asset_vol[i]*asset_vol[j])

    return Correlation



# =============================================================================
# type of possible payoffs
# =============================================================================


#choose the best asset price within the basket at the maturity
def Simple_chooser(asset_price):
    return np.max(asset_price)

#choose the best call option payoff on the different assets 
def Chooser_call(asset_price,strike_price):
    selection=[]
    for s in asset_price:
         if s>strike_price:
              selection.append(s-strike_price)
    if len(selection)==0:
         return 0
    else:
         return np.max(selection)

#choose the best put option payoff on the different assets
def Chooser_put(asset_price,strike_price):
    selection=[]
    for s in asset_price:
         if s<strike_price:
              selection.append(strike_price-s)
    if len(selection)==0:
         return 0
    else:
         return np.max(selection)


# =============================================================================
# we use the explicit formula of BS price at the maturity-- we get the multivariate normal variable and we do the computation.
#the function return the price of different assets corresponding to one sample of the multivariate normal variable

# =============================================================================

    
def price_T(r,T,asset_price_0,correlation_matrix):
    N_assets=len(asset_price_0)
    # we could generate a iid normal variable X, and using cholesky decomposition of the correlation matrix U'U, we get the multivariate gaussian vector with the underlying correlation
    #matrix as a combinaison of the iid varible Z=U'X

    x=np.random.multivariate_normal(mean=np.zeros(N_assets),cov=correlation_matrix, size=1).T
    S_T=[]

    for i in range(N_assets):
        S=asset_price_0[i]*np.exp((r-1/2*correlation_matrix[i][i])*T+np.sqrt(T)*x[i][0])
       
        S_T.append(S)
    return S_T

# =============================================================================
# We value the option on the basket based on Montee carlo simulation, using a given payoff
# =============================================================================

def Option_value(r,T,asset_price_0,correlation_matrix,N_simulation,payoff_type=1,strike_price=90):
    #we generate N_simulation  scenario, we do the average of different payoffs.
   option_value=0
   if payoff_type==1:
         for i in range(N_simulation):
                  S_T=price_T(r,T,asset_price_0,correlation_matrix)
                  payoff=Simple_chooser(S_T)
                  option_value+=payoff
   elif payoff_type==2:
          for i in range(N_simulation):
                  S_T=price_T(r,T,asset_price_0,correlation_matrix)
                  payoff=Chooser_call(S_T,strike_price)
                  option_value+=payoff
   else:
          for i in range(N_simulation):
                  S_T=price_T(r,T,asset_price_0,correlation_matrix)
                  payoff=Chooser_put(S_T,strike_price)
                  option_value+=payoff
   #the discounting of the option value
   option_value=option_value/N_simulation*(np.exp(-r*T))

   return option_value


#generate a correlation matrix for the test
    
def generate_corr(N_assets):
     X=np.random.uniform(low=0.1, high=1, size=N_assets*N_assets).reshape(N_assets, N_assets)
     corr= (X+X.T)/2
     return corr

# =============================================================================
# -----------------Test----------------------------------------------------------
# =============================================================================
     
N_assets=10

Correlation=generate_corr(N_assets)


r=0.20
T=1
asset_price_0=np.random.uniform(low=10, high= 100, size=N_assets).reshape(N_assets,1)

N_simulation=100


Value0=Option_value(r,T,asset_price_0,Correlation,N_simulation,2,90)

print("The value of the option on the underlying basket, using the payoff of Simpe_chooser is %f",Value0)



#Call option pricing

S0=[100,90]
vol=np.array([[0.6,0],[0,0.7]])
strike_price=97
call=Option_value(r,T,S0,vol,N_simulation,2,strike_price)
print(call)

