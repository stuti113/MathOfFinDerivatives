#-*- coding: utf-8 -*-



import numpy as np
import matplotlib.pyplot as plt
import math
def bs_finite_exp():
    print()
    print("Black-Scholes PDE, Numerical Estimation, Explicit Finite Differences")
    print()
    #Set the current stock price S curr, strike price K, volatility sigma, risk-free rate r, time to expiration T
    S_curr = 105.99
    K = 105.00
    sigma = 0.3128
    r = 0.0017
    T = 351.0/365.0
    divYield = 0.017
    #Set number of asset price steps J
    J = 200
    #Set number of time steps M and time increment deltaT
    M = 4000
    #Set size of deltaS and deltaT
    S_max = 4*S_curr
    deltaS = S_max/J
    deltaT = T/M
    
    print("Strike Price K = %s" %K)
    print("Volatility sigma = %s" %sigma)
    print("Risk-Free Rate r = %s" %r)
    print("Time to expiration T = %s" %T)
    print("Dividend Yield divYield = %s" %divYield)
    print()
    
    print("Asset Price Steps J = %s" %J)
    print("Time Increment Steps M = %s" %M)
    print("Size of Asset Price Steps = %s" %deltaS)
    print("Size of Time Steps = %s" %deltaT)
    print()
    
    #Create (J-1) X (J-1) Matrix A
    #Create column vector of asset prices S
    #Create column vector of option prices C_hat
    A=np.zeros((J-1,J-1),dtype=float)
    S=np.zeros((J-1,1),dtype=float)
    C_hat_put =np.zeros((J-1,1),dtype=float)
    C_hat_call =np.zeros((J-1,1),dtype=float)
    
    S[0] = deltaS
    j=1 #row counter
    for j in range(1,J-1):
         S[j] = S[j-1] + deltaS

    print("S = %s" %S)
    print()
    
    j=0 #row counter
    
    for j in range(0,J-1):
        C_hat_call[j] = max(S[j]-K,0)
        C_hat_put[j] = max(K-S[j],0)    
    print("C_hat_call = %s" %C_hat_call)
    print("C_hat_put = %s" %C_hat_put)
    print()
    
    j=0 #row counter
    k=0 #column counter
    for j in range (0,J-1):
        for k in range (0,J-1):
            if k == j:
                A[j,j] = 1 - (sigma**2)*((j+1)**2)*deltaT - r*deltaT
            elif k == j-1:
                A[j,k] = 0.5 * ((sigma**2)*((j+1)**2)*deltaT - (r - divYield)*(j+1)*deltaT)
            elif k == j+1:
                A[j,k] = 0.5 * ((sigma**2)*((j+1)**2)*deltaT + (r - divYield)*(j+1)*deltaT)
     
    print("A = %s" %A)
    print()              
                     
    C_0 = 0 #for Call
    C_J = 0 #for Put
    C_start_call = C_hat_call
    C_start_put = C_hat_put
    m=1 #time increment counter
    while m <= M:
        #Boundary values calculated separately
        #Must incorporate lowest and highest option prices, which are not in matrix
        C_min_Call = 0.5*((sigma**2)*((1)**2)*deltaT - (r - divYield)*(1)*deltaT)*C_0 +(1 - (sigma**2)*((1)**2)*deltaT - r*deltaT)*C_hat_call[0] + 0.5*((sigma**2)*((1)**2)*deltaT + (r - divYield)*(1)*deltaT)*C_hat_call[1]
        C_max_Call = 0.5*((sigma**2)*((J-1)**2)*deltaT - (r - divYield)*(J-1)*deltaT)*C_hat_call[J-3] +(1 - (sigma**2)*((J-1)**2)*deltaT - r*deltaT)*C_hat_call[J-2] + 0.5*((sigma**2)*((J-1)**2)*deltaT + (r - divYield)*(J-1)*deltaT)*(S_max - K*math.exp(-r*m*deltaT))
        #Perform matrix multiplication Ac_m = c_m+1
        C_hat_call = A.dot(C_hat_call)
        #Update boundary values
        C_hat_call[0] = C_min_Call
        C_hat_call[J-2] = C_max_Call
        #Capture option prices during backward walk for graphing
        if m == M/4:
            C1_call = C_hat_call
        elif m == 2*M/4:
            C2_call = C_hat_call
        elif m == 3*M/4:
            C3_call = C_hat_call
        elif m == 4*M/4:
            C4_call = C_hat_call
            
        m=m+1
    print("C_hat =%s" %C_hat_call)
    print() 
    print("Asset Price = %s" %S[(J-2)/4])
    print("Call Price = %s"%C_hat_call[(J-2)/4])
    plt.plot(S,C_start_call,'b--',S,C1_call,'r--',S,C2_call,'g--',S,C3_call,'r--',S,C4_call,'b--')
   
    
    
    m=1 
    while m <= M:
        C_min_put = 0.5*((sigma**2)*(1**2)*deltaT - (r - divYield)*(1)*deltaT)*(K*math.exp(-r*m*deltaT))+ A[0,0]*C_hat_put[0] + A[0,1]*C_hat_put[1]
        C_max_put = A[J-2,J-3]*C_hat_put[J-3] + A[J-2,J-2]*C_hat_put[J-2] + 0.5*((sigma**2)*((J-1)**2)*deltaT + (r - divYield)*(J-1)*deltaT)*C_J
        C_hat_put = A.dot(C_hat_put)
            #Update boundary values
        C_hat_put[0] = C_min_put
        C_hat_put[J-2] = C_max_put
            #Capture option prices during backward walk for graphing
        if m == M/4:
            C1_put = C_hat_put
        elif m == 2*M/4:
            C2_put = C_hat_put
        elif m == 3*M/4:
            C3_put = C_hat_put
        elif m == M:
            C4_put = C_hat_put
                
        m=m+1
    print("C_hat_put =%s" %C_hat_put)
    print() 
    print("Asset Price = %s" %S[(J-2)/4])
    print("Put Price = %s" %C_hat_put[(J-2)/4])
    #plt.plot(S,C_start_call,'b--',S,C1,'r--',S,C2,'g--',S,C3,'r--',S,C4,'b--')
    plt.plot(S,C_start_put,'b--',S,C1_put,'r--',S,C2_put,'g--',S,C3_put,'r--',S,C4_put,'b--')
    
 #main program starts here
bs_finite_exp()

                
        
        