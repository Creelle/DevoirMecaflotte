import numpy as np

"""
fanno flows
"""
#1) critical length function

def critical_lenghth(M):
    return 1/1.4*((1-M**2)/M**2+1.2*np.log(1.2*M**2/(1+0.2*M**2)))
print(critical_lenghth(0.24))
def inv_critical_length(a,M_2): # lambda*(x2-x1)/Dh = f(M1)-f(M2) so a either  lambda*(x2-x1)/Dh+f(M2) or a= f(M1)- lambda*(x2-x1)/Dh
    """
    returns M_2 is M**2  and not M ==> a changer je pense qu il y a un probleme
    """
    iter=1
    error=1
    ones= np.ones(len(M_2))

    while iter<1000 and error>0.00001 :
        M_2new = (1+1.4*a-1.2*np.log(1.2*M_2**2/(ones+0.2*M_2**2)))**(-1/2)
        error=np.linalg.norm(M_2new-M_2)
        iter +=1
        M_2=M_2new
    if(iter==1000):
        print('max value reached')
    return M_2
def Re(rho,u,d,T):
    mu_ref = 1.716e-5
    S=111.0 #K
    Tref = 273.15
    mu = mu_ref*(T/Tref)**1.5*(Tref+S)/(T+S)
    return rho*u*d/mu

def perte_lambda(Re,lambda_esti):
    iter3=1
    error3=1
    value = 1/np.sqrt(lambda_esti)
    while iter3<100 and error3>0.00001:
        new = -3.0*np.log(2.03*value/Re)/np.log(10)
        error3 = abs(new-value)
        iter3 +=1
        value = new
    return 1/value**2

def inv_air_ratio(air_ratio, M): # for subsonic speed
    iter = 1
    error = 1
    ones = np.ones(len(M))
    while iter<10000 and error>0.00001 :
        M_new = air_ratio*((ones+0.2*M**2)/1.2)**(3)
        error =np.linalg.norm(M_new-M)
        iter +=1
        M = M_new
    return M

def inv_air_ratio_sup(air_ratio,M): #for supersonic air_ratio =  A*/A
    iter = 1
    error = 1
    ones=np.ones(len(M))
    while iter<10000 and error >0.00001 :
        M_new = np.sqrt(2/0.4*(1.2*(M/air_ratio)**(0.8/2.4)-ones))
        error =np.linalg.norm(M_new-M)
        iter +=1
        M=M_new
    return M
