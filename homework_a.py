import numpy as np
from nozzle import nozzle_height
import useful
import matplotlib.pyplot as plt

"""
données de base
"""
h_throat = 0.003 # 3 [mm]
h_duct = 0.0075 #7.5 [mm]
l_nozzle = 0.12 # [m]
l_duct = 1.5 # [m]
de = 2*h_duct
dt = 0.006 # 6 [mm]
b = 0.05 #50 [mm]
Dh_duct = 4*h_duct*b/(b+2*h_duct) #4*A/perimeter
p_a = 1.01325 # 1 [atm] = 1.01325 bar
T0 = 300 # [K]
# Gas = air (assuming air as a perfect gas)
gamma = 1.4
R = 287.1 # [J/kg/K]

p_e = p_a

"""
1) find po in order to have sonic throat and subsonic after and before

"""
# Throat conditions = sonic conditions
T_star = 2*T0/(gamma+1)
c_star = np.sqrt(gamma*R*T_star)


# Beginning of duct (fanno flow)
Mi = np.array([0.5])
air_ratio = dt/de
Mi = useful.inv_air_ratio(air_ratio,Mi)
print('Mi = ',Mi)
T_i = T0/(1 + 0.5*(gamma-1)*Mi**2)
u_i = Mi*np.sqrt(gamma*R*T_i)

# End of the duct
# Assuming M_e = 0.7
Me = np.array([0.7])

#first estimation of lambda
T_e = T0/(1 + 0.5*(gamma-1)*Me**2)
u_e = Me* np.sqrt(gamma*R*T_e)
rho_e = p_e*10**(5)/R/T_e
Re_de = useful.Re(rho_e, u_e, de, T_e)
lambda_e = 1/(2.33*Re_de**0.1)**2

# Iterating to find the right M_e
error = 1
iter = 1
while(iter<100 and error>0.0001):
    a = useful.critical_lenghth(Mi)-0.5*lambda_e*l_duct/Dh_duct
    M_new = useful.inv_critical_length(a,Me)
    iter += 1
    error = abs(M_new - Me)
    Me = M_new
    T_e = T0 / (1 + 0.5 * (gamma - 1) * Me ** 2)
    u_e =  Me* np.sqrt(gamma * R * T_e)
    rho_e = p_e*10**(5) / R / T_e
    Re_de = useful.Re(rho_e, u_e, de, T_e)
    lambda_e = useful.perte_lambda(Re_de, lambda_e)



# a) conditions outlet duct
p_star_e = p_e*Me*np.sqrt((1+0.2*Me**2)/1.2)  # [bar]
p_0e = p_e * (1 + 0.5 * (gamma - 1) * Me ** 2) ** (gamma / (gamma - 1))
p_0_star_e = p_0e * Me * ((0.5*(gamma+1)/(1 + 0.5*(gamma-1)*Me**2)))**(0.5*(gamma+1)/(gamma-1))
print('conditions outlet duct: ',Me,u_e,T_e,p_e,p_0e)

# b) conditions inlet duct
rho_i = rho_e*u_e/u_i #conservation of mass
p_i = rho_i*T_i*R*10**(-5)#[bar]
#p_i2 = p_star_e/Mi * np.sqrt(1.2/(1+0.2*Mi**2)) # [bar]
p_0i = p_i * (1+0.5*(gamma-1)*Mi**2)**(gamma/(gamma-1))
#p_0i2 = p_0_star_e / Mi * ((1 + 0.5*(gamma-1)*Mi**2)/(0.5*(gamma+1)))**(0.5*(gamma+1)/(gamma-1))
print('conditions inlet duct', Mi,u_i,T_i,p_i,p_0i)

# c) chequ because p0_star (and p_star) are conserved trough the duct
p_0_star_i = p_0i*Mi*(0.5*(gamma+1)/(1 + 0.5*(gamma-1)*Mi**2))**(0.5*(gamma+1)/(gamma-1))
print('p_0_star chequ = ', p_0_star_i, p_0_star_e)

# d) conditions at the throat
p_star = p_0i*(T_star/T0)**(1.4/0.4)
rho_star = p_star*10**(5)/(R*T_star)
print('conditions at throat', 1, c_star, T_star, p_star, p_0i)

# e) answers
print('stagnation pressure [bar]: ',p_0i)
Qm = u_e*rho_e*b*2*h_duct#kg/s
print('massflow [kg/s]',Qm)

"""
2) compute M,p,p0 along the nozzle and the duct
"""

po_star = p_0_star_e
#a) nozzle
x=np.linspace(0,120,1000)
ones=np.ones(len(x))
h= np.zeros(len(x))
for i in range(len(h)):
     h[i]= nozzle_height(x[i])/1000
air_ratios = h_throat/h
Mx = useful.inv_air_ratio(air_ratios,0.5*np.ones(len(x)))
px = p_0i/((1+0.2*Mx**2)**(1.4/0.4))

np.savetxt("results\homework_a_nozzle.txt",np.stack((x,Mx,px,p_0i*ones)))
# fig, ax = plt.subplots()
# ax.set_title('a nozzle')
# ax.plot(x,Mx,'-b',label='M')
# ax.plot(x,px,'-r',label='p[bar]')
# ax.plot(x,p_0i*ones,'-c',label='p0[bar]')
# ax.legend()
# ax.grid()


#b) duct
x2=np.linspace(0,l_duct,1000)
ones=np.ones(len(x2))
a = -0.5*lambda_e*x2/Dh_duct+ones*useful.critical_lenghth(Mi)
Mx_d = useful.inv_critical_length(a,0.2*ones)
pox = po_star*(1/Mx_d*((ones+0.2*Mx_d**2)/1.2)**(1.2/0.4))
px = pox/(ones+0.2*Mx_d**2)**(1.4/0.4)
np.savetxt("results\homework_a_duct.txt",np.stack((x2,Mx_d,px,pox)))#attention au stack

# fig2=plt.figure()
# ax2 =plt.subplot(311)
# ax2.set_title('a duct')
# ax2.plot(x2,Mx_d,'-b',label='M')
# ax22=plt.subplot(312)
# ax22.plot(x2,px,'-r',label='p[bar]')
# ax23=plt.subplot(313)
# ax23.plot(x2,pox,'-c',label='p0[bar]')
# ax2.legend()
# ax2.grid()
# ax22.legend()
# ax22.grid()
# ax23.legend()
# ax23.grid()
#
plt.show()
