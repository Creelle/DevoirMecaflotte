import numpy as np
from nozzle import nozzle_height
import useful
from matplotlib import pyplot as plt
"""
données de base
"""

h_throat = 0.003 #3mm
h_duct = 0.0075 #7.5 mm
length_nozzle = 0.12
l_duct = 1.5#58.29722
de = 2*h_duct
b = 0.05 #50mm
Dh_duct = 4*h_duct*b/(b+2*h_duct) #4*A/perimeter
x_shock= 0.12
h_shock = nozzle_height(x_shock*1000)/1000 #[m]

pa = 1.01325 # 1 atm = 1.01325 bar
gamma = 1.4
R=287.1 #J/kg/K
T0 = 300 #K

#1) Determine conditions at the throat = sonic conditions A_t = A*
T_star = T0/(gamma+1)*2
c_star = np.sqrt(gamma*R*T_star)

# Shock relations :
#  -1- Isentropic between reservoir and sh1
#  -2- Isentropic between sh2 and end of the nozzle
#2) determiner les conditions juste avant le shock
air_ratio = h_throat/h_shock #ratio between the throat and the first shock
M1 = useful.inv_air_ratio_sup(air_ratio,np.array([2.0]))
T1 =  T0/(1+0.2*M1**2)
u1= M1*np.sqrt(1.4*R*T1)

#3) determiner la vitesse apres le shock en utilisant les shocks relations
# M2=Mi vu que le shock se produit a l inlet du duct
M2 = np.sqrt((1+0.2*M1**2)/(1.4*M1**2-0.2))
T2 = T0/(1+0.2*M2**2)
u2 = M2*np.sqrt(1.4*R*T2)

#4) fanno flow , iterer pour trouver Me qui nous donne M2 a l'inlet
pe = pa #outlet subsonic

Me = np.array([0.7]) #premiere estimation
Te = T0/(1 + 0.5*(gamma-1)*Me**2)
ue = Me*np.sqrt(gamma*R*Te)
rho_e = pe*10**(5)/R/Te
Re_de = useful.Re(rho_e, ue, de, Te)
lambda_e = 1/(2.33*Re_de**0.1)**2 #premiere estimation de lambda_e
error=1
iter = 1
while error>0.0001 and iter<100 :
    a = useful.critical_lenghth(M2)-0.5*lambda_e*l_duct/Dh_duct
    Me_new = useful.inv_critical_length(a,Me)
    error =abs(Me_new-Me)
    iter +=1
    Me=Me_new

    Te = T0/(1+0.2*Me**2)
    ue = Me*np.sqrt(1.4*R*Te)
    rho_e = pe*10**(5)/(R*Te)
    Re = useful.Re(rho_e,ue,de,Te)
    lambda_e = useful.perte_lambda(Re,lambda_e)

#5) conditions at outlet
poe = pe*(1+0.2*Me**2)**(1.4/0.4)
poe_star = poe/(1/Me*((1+0.2*Me**2)/1.2)**(3)) # poe_star constant dans le duct (fanno flow)
p_star = pe*Me/(1.2/(1+0.2*Me**2))**(1.2/0.4)
print('conditions outlet duct: ',Me,ue,Te,pe,poe)

#6) conditions inlet duct - after shock
rho_2 = ue*rho_e/u2
p2 =  rho_2*R*T2*10**(-5)
# p22 = po2/(1+0.2*M2**2)**(1.4/0.4)
po2 = p2*(T0/T2)**(1.4/0.4)
po2_star = po2/(1/M2*((1+0.2*M2**2)/1.2)**(3))
p_star2 = p2*M2/(1.2/(1+0.2*M2**2))**(1.2/0.4)
print('po_star chequ', poe_star, po2_star)
print('pstar chequ',p_star,p_star2)
print('conditions inlet duct - after shock', M2,u2,T2,p2,po2)


#7) conditions inlet duct - before shock
p1 = p2*1.2/(1.4*M1**2-0.2) #p1_2  =p2*(1+1.4*M2**2)/(1+1.4*M1**2)
po1 = p1*(T0/T1)**(1.4/0.4)
rho_1 = p1*10**(5)/(R*T1)
print('massflow chequ at the shock', rho_1*u1, rho_2*u2)
print('conditions before shock', M1,u1,T1,p1,po1)

#8) conditions at the throat
p_star = po1*(T_star/T0)**(1.4/0.4)
rho_star = p_star*10**(5)/(R*T_star)
print('conditions at throat',1.,c_star,T_star,p_star,po1)
print('reponse a la question b po = ',po1)

#9) answers
print('stagnation pressure [bar]: ',po1)
Qm = ue*rho_e*b*2*h_duct#kg/s
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
px1 = po1/((ones+0.2*Mx1**2)**(1.4/0.4))

# pres du throat subsonic
x1b =np.linspace(28,30,500)
ones=np.ones(len(x1b))
h1b= np.zeros(len(x1b))
for i in range(len(h1b)):
    h1b[i]= nozzle_height(x1b[i])/1000
air_ratios1b = h_throat/h1b

Mx1_b = useful.inv_air_ratio(air_ratios1b,1.0*np.ones(len(x1b)))
px1_b = po1/((ones+0.2*Mx1_b**2)**(1.4/0.4))

#pres du throat supersonic
x1c = np.linspace(30,32,500)
h1c= np.zeros(len(x1c))
for i in range(len(h1c)):
    h1c[i]= nozzle_height(x1c[i])/1000
air_ratios1c = h_throat/h1c

Mx_1c = useful.inv_air_ratio_sup(air_ratios1c,1.0*np.ones(len(x1c)))
px_1c = po1/((ones+0.2*Mx_1c**2)**(1.4/0.4))

#avant le shock supersonic
x1d = np.linspace(32,120,1000)
h1d= np.zeros(len(x1d))
for i in range(len(h1d)):
    h1d[i]= nozzle_height(x1d[i])/1000
air_ratio1d = h_throat/h1d

Mx_1d = useful.inv_air_ratio_sup(air_ratio1d,1.5*np.ones(len(x1d)))
px_1d = po1/((1+0.2*Mx_1d**2)**(1.4/0.4))

np.savetxt("results\homework_c_nozzle.txt",np.stack((np.concatenate((x1,x1b,x1c,x1d)),np.concatenate((Mx1,Mx1_b,Mx_1c,Mx_1d)),np.concatenate((px1,px1_b,px_1c,px_1d)))))
#b) nozzle plot
# xa = np.linspace(0,120,1000)
# fig, ax = plt.subplots()
# ax.set_title('c nozzle')
# ax.plot(x1,Mx1,'-b',label='M')
# ax.plot(x1b,Mx1_b,'-b')
# ax.plot(x1,px1,'-r',label='p[bar]')
# ax.plot(x1b,px1_b,'-r')
# ax.plot(x1c,Mx_1c,'-b')
# ax.plot(x1c,px_1c,'-r')
# ax.plot(x1d,Mx_1d,'-b')
# ax.plot(x1d,px_1d,'-r')
# ax.plot(xa,po1*np.ones(len(xa)),'-c',label='p0[bar]')
# ax.legend()
# ax.grid()

#c) duct
x4=np.linspace(0,l_duct,1000)
ones=np.ones(len(x4))
a = -0.5*lambda_e*x4/Dh_duct+ones*useful.critical_lenghth(M2)
Mx_d = useful.inv_critical_length(a,0.2*ones)
pox = poe_star*(1/Mx_d*((ones+0.2*Mx_d**2)/1.2)**(1.2/0.4))
px = pox/(ones+0.2*Mx_d**2)**(1.4/0.4)

np.savetxt("results\homework_c_duct.txt",np.stack((x4,Mx_d,px,pox)))#attention au stack

# fig2=plt.figure()
# ax2 =plt.subplot(311)
# ax2.set_title('c duct')
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
