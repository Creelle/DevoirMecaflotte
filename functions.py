import numpy as np
gamma = 1.4

def ratio_A_star_over_A_sub(A_star, A, M): # for subsonic speed
    iter = 1
    error = 1
    while (iter<50 and error>0.0001):
        M_new = A_star/A *((1+0.5*(gamma-1)*M**2)/(0.5*(gamma+1)))**(0.5*(gamma+1)/(gamma-1))
        error = abs(M_new - M)
        iter += 1
        M = M_new
    return M, error, iter

def ratio_A_star_over_A_super(A_star, A, M):
    iter = 1
    error = 1
    while(iter<50 and error>0.00001):
        M_new = np.sqrt(2/(gamma-1) * ((0.5*(gamma+1)*(A*M/A_star)**(2*(gamma-1)/(gamma+1))) - 1))
        error = abs(M_new - M)
        iter += 1
        M = M_new
    return M, error, iter
print(ratio_A_star_over_A_super(1,1.84,2))

def f(M):
    return 1/gamma*((1-M**2)/M**2 + 1.2*np.log(1.2*M**2/(1+0.2*M**2)))

def Re(rho,u,d,T):
    mu_ref = 1.716e-5 # [Pa.s]
    S = 111.0 # [K]
    T_ref = 273.15
    mu = mu_ref*(T/T_ref)**1.5*(T_ref + S)/(T + S)
    return rho*u*d/mu

def lambd(Re,lambda_esti):
    iter = 1
    error = 1
    value = 1/np.sqrt(lambda_esti)
    while(iter<50 and error>0.00001):
        new = -3.0*np.log(2.03*value/Re)/np.log(10)
        error = abs(new - value)
        iter += 1
        value = new
    return 1/value**2, error, iter

def shock_M_sh2(M_sh1):
    return np.sqrt((1 + 0.5*(gamma-1)*M_sh1**2)/(gamma*M_sh1**2 - 0.5*(gamma-1)))

def shock_T_sh2(T_sh1, M_sh1):
    return T_sh1 * (gamma*M_sh1**2 - 0.5*(gamma-1))/(0.5*(gamma+1)) * (1 + 0.5*(gamma-1)*M_sh1**2)/(0.5*(gamma+1)*M_sh1**2)
