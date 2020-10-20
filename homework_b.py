import numpy as np
from nozzle import nozzle_height
import useful
import matplotlib.pyplot as plt

"""
1) données de base
"""

h_throat = 0.003 # 3 [mm]
h_duct = 0.0075 #7.5 [mm]
l_nozzle = 0.12 # [m]
l_duct = 1.5 # [m]
xe = l_duct # [m]
de = 2*h_duct
dt = 0.006 # 6 [mm]
b = 0.05 # 50 [mm]
Dh_duct = 4*h_duct*b/(b+2*h_duct) # 4*A/P
p_a = 1.01325 # 1.01325[bar] = 1[atm]
T0 = 300 # [K]
# Gas = air (assuming air as a perfect gas)
gamma = 1.4
R = 287.1 # [J/kg/K]

p_e = p_a
x_sh = 0.07 # [m]
h_sh = nozzle_height(x_sh*1000)/1000 #[m]

"""
2) resolution de l exercice'
"""

# Throat conditions = sonic conditions
T_star = 2*T0/(gamma+1)
c_star = np.sqrt(gamma*R*T_star)

# Shock relations :
#  -1- Isentropic between reservoir and sh1
#  -2- Isentropic between sh2 and end of the nozzle

#1) determiner les conditions juste avant le shock
M_sh1 = np.array([1.5])
air_ratio1 = h_throat/h_sh #ratio between the throat and the first shock
M_sh1 = useful.inv_air_ratio_sup(air_ratio1,np.array([1.5]))
T_sh1 = T0/(1 + 0.5*(gamma-1)*M_sh1**2)
u_1 = M_sh1*np.sqrt(gamma*R*T_sh1)

#2) determiner les conditions apres le shock en utilisant les shocks relations
M_sh2 = np.sqrt((1 + 0.5*(gamma-1)*M_sh1**2)/(gamma*M_sh1**2 - 0.5*(gamma-1)))
T_sh2 = T0/(1 + 0.5*(gamma-1)*M_sh2**2)
u_2 = M_sh2*np.sqrt(gamma*R*T_sh2)
h_star2 = h_sh * M_sh2 * (0.5*(gamma+1) / (1 + 0.5*(gamma-1)*M_sh2**2))**(0.5*(gamma+1)/(gamma-1))
air_ratio2 = h_star2/h_duct

#3)  determiner les conditions the end of nozzle using air_ratio2
Mi =  useful.inv_air_ratio(air_ratio2,np.array([0.5]))
T_i = T0/(1 + 0.2*Mi**2)
c_i = np.sqrt(gamma*R*T_i)
u_i = Mi*c_i


# 4) Fanno flow duct assuming first that Me = 0.7 (Me has to be greater than M_i since the flows accelerates
# Need to calculate the critical length in order to know if we'll be sonic at the exhaust ?
# Maybe we'll have to iterate from here up to...
# Assuming Me = 0.7 and first estimation of lambda_e
Me = np.array([0.7])
T_e = T0/(1 + 0.5*(gamma-1)*Me**2)
u_e = Me*np.sqrt(gamma*R*T_e)
rho_e = p_e*10**(5)/R/T_e
Re_de = useful.Re(rho_e, u_e, de, T_e)
lambda_e = 1/(2.33*Re_de**0.1)**2 #first estimation

# Iterating to find the right Me
error = 1
iter = 1
while(iter<100 and error>0.0001):
    a = useful.critical_lenghth(Mi)-0.5*lambda_e*l_duct/Dh_duct
    M_new = useful.inv_critical_length(a,Me)
    #M_new = 1 / np.sqrt(1 + (fun.f(Mi) - 0.5*lambda_e * xe/Dh_duct)*gamma - 0.5*(gamma+1)*np.log(1.2*Me**2 / (1+0.2*Me**2)))

    iter += 1
    error = abs(M_new - Me)
    Me = M_new
    T_e = T0 / (1 + 0.5 * (gamma - 1) * Me ** 2)
    u_e =  Me* np.sqrt(gamma * R * T_e)
    rho_e = p_e*10**(5) / R / T_e
    Re_de = useful.Re(rho_e, u_e, de, T_e)
    lambda_e = useful.perte_lambda(Re_de, lambda_e)


#5) conditions at the outlet
p_star = p_e*Me*np.sqrt((1+0.2*Me**2)/1.2)  # [bar] constant in the duct
p_0e = p_e * (1 + 0.5 * (gamma - 1) * Me ** 2) ** (gamma / (gamma - 1))
p_0_star = p_0e * Me * ((0.5*(gamma+1)/(1 + 0.5*(gamma-1)*Me**2)))**(0.5*(gamma+1)/(gamma-1)) #bar constant in the duct
print('conditions outlet duct: ',Me,u_e,T_e,p_e,p_0e)

# 6) conditions inlet duct
rho_i = rho_e*u_e/u_i #conservation of mass
p_i = rho_i*T_i*R*10**(-5)#[bar]
#p_i2 = p_star/Mi * np.sqrt(1.2/(1+0.2*Mi**2)) # [bar]
p_0i = p_i * (1+0.5*(gamma-1)*Mi**2)**(gamma/(gamma-1))
#p_0i2 = p_0_star / Mi * ((1 + 0.5*(gamma-1)*Mi**2)/(0.5*(gamma+1)))**(0.5*(gamma+1)/(gamma-1))
print('conditions inlet duct', Mi,u_i,T_i,p_i,p_0i)

# 7) chequ because p0_star (and p_star) are conserved trough the duct
p_0_star_i = p_0i*Mi*(0.5*(gamma+1)/(1 + 0.5*(gamma-1)*Mi**2))**(0.5*(gamma+1)/(gamma-1))
print('p_0_star chequ = ', p_0_star_i, p_0_star)

#8)  conditions before (1) and after(2) the shock
# p_02 = p_0i
# p_01 = p_0

p_sh2 = p_0i/(1+0.2*M_sh2**2)**(1.4/0.4)
#p_sh22 = p_i * ((1 + 0.5*(gamma-1)*Mi**2)/(1 + 0.5*(gamma-1)*M_sh2**2))**(gamma/(gamma-1))
print('conditions after shock', M_sh2,u_2,T_sh2,p_sh2,p_0i)

p_sh1 = p_sh2 * 0.5*(gamma+1)/(gamma*M_sh1**2 - 0.5*(gamma-1))
p_01 = p_sh1 * (1 + 0.5*(gamma-1)*M_sh1**2)**(gamma/(gamma-1))
print('conditions before shock', M_sh1,u_1,T_sh1,p_sh1,p_01)

#9) conditions at the throat
p_star = p_01*(T_star/T0)**(1.4/0.4)
print('conditions at throat',1.,c_star,T_star,p_star,p_01)

# 10) answers
print('stagnation pressure [bar]: ',p_01)
Qm = u_e*rho_e*b*2*h_duct#kg/s
print('massflow [kg/s]',Qm)

"""
3) compute M,p,p0 along the nozzle and the duct
"""

#a) nozzle
# de 0  jusque pres du chock
x1 =np.linspace(0,28,1000)
ones=np.ones(len(x1))
h1= np.zeros(len(x1))
for i in range(len(h1)):
    h1[i]= nozzle_height(x1[i])/1000
air_ratios1 = h_throat/h1

Mx1 = useful.inv_air_ratio(air_ratios1,0.5*np.ones(len(x1)))
px1 = p_01/((ones+0.2*Mx1**2)**(1.4/0.4))

# pres du throat subsonic
x1b =np.linspace(28,30,500)
ones=np.ones(len(x1b))
h1b= np.zeros(len(x1b))
for i in range(len(h1b)):
    h1b[i]= nozzle_height(x1b[i])/1000
air_ratios1b = h_throat/h1b

Mx1_b = useful.inv_air_ratio(air_ratios1b,1.0*np.ones(len(x1b)))
px1_b = p_01/((ones+0.2*Mx1_b**2)**(1.4/0.4))

#pres du throat supersonic
x1c = np.linspace(30,32,500)
h1c= np.zeros(len(x1c))
for i in range(len(h1c)):
    h1c[i]= nozzle_height(x1c[i])/1000
air_ratios1c = h_throat/h1c

Mx_1c = useful.inv_air_ratio_sup(air_ratios1c,1.0*np.ones(len(x1c)))
px_1c = p_01/((ones+0.2*Mx_1c**2)**(1.4/0.4))

#avant le shock supersonic
x1d = np.linspace(32,70,1000)
h1d= np.zeros(len(x1d))
for i in range(len(h1d)):
    h1d[i]= nozzle_height(x1d[i])/1000
air_ratio1d = h_throat/h1d

Mx_1d = useful.inv_air_ratio_sup(air_ratio1d,1.5*np.ones(len(x1d)))
px_1d = p_01/((1+0.2*Mx_1d**2)**(1.4/0.4))

#apres le shook
x2 =np.linspace(70,120,1000)
h2= np.zeros(len(x2))
for i in range(len(h2)):
    h2[i]= nozzle_height(x2[i])/1000
air_ratios2 = h_star2/h2

Mx2 = useful.inv_air_ratio(air_ratios2,0.5*np.ones(len(x2)))
px2 = p_0i/((1+0.2*Mx2**2)**(1.4/0.4))

np.savetxt("results\homework_b_nozzle.txt",np.stack((np.concatenate((x1,x1b,x1c,x1d,x2)),np.concatenate((Mx1,Mx1_b,Mx_1c,Mx_1d,Mx2)),np.concatenate((px1,px1_b,px_1c,px_1d,px2)))))

#b) nozzle plot
# xa = np.linspace(0,70,1000)
# xb = np.linspace(70,120,1000)
# fig, ax = plt.subplots()
# ax.set_title('b nozzle')
# ax.plot(x1,Mx1,'-b',label='M')
# ax.plot(x1b,Mx1_b,'-b')
# ax.plot(x1,px1,'-r',label='p[bar]')
# ax.plot(x1b,px1_b,'-r')
# ax.plot(x1c,Mx_1c,'-b')
# ax.plot(x1c,px_1c,'-r')
# ax.plot(x1d,Mx_1d,'-b')
# ax.plot(x1d,px_1d,'-r')
# ax.plot(x2,Mx2,'-b')
# ax.plot(x2,px2,'-r')
# ax.plot(xb,p_0i*np.ones(len(xb)),'-c',label='p0[bar]')
# ax.plot(xa,p_01*np.ones(len(xa)),'-c')
# ax.legend()
# ax.grid()

#c) duct
x4=np.linspace(0,l_duct,1000)
ones=np.ones(len(x4))
a = -0.5*lambda_e*x4/Dh_duct+ones*useful.critical_lenghth(Mi)
Mx_d = useful.inv_critical_length(a,0.2*ones)
pox = p_0_star*(1/Mx_d*((ones+0.2*Mx_d**2)/1.2)**(1.2/0.4))
px = pox/(ones+0.2*Mx_d**2)**(1.4/0.4)

np.savetxt("results\homework_b_duct.txt",np.stack((x4,Mx_d,px,pox)))#attention au stack

# fig2=plt.figure()
# ax2 =plt.subplot(311)
# ax2.set_title('a duct')
# ax2.plot(x4,Mx_d,'-b',label='M')
# ax22=plt.subplot(312)
# ax22.plot(x4,px,'-r',label='p[bar]')
# ax23=plt.subplot(313)
# ax23.plot(x4,pox,'-c',label='p0[bar]')
# ax2.legend()
# ax2.grid()
# ax22.legend()
# ax22.grid()
# ax23.legend()
# ax23.grid()
#
# plt.show()
