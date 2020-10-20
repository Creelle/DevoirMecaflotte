import numpy as np
from nozzle import nozzle_height
import useful
from matplotlib import pyplot as plt

#données de base

h_throat = 0.003 #3mm
h_duct = 0.0075 #7.5 mm
length_nozzle = 0.12
length_duct = 1.5#58.29722
d = 2*h_duct
b = 0.05 #50mm
Dh_duct = 4*h_duct*b/(b+2*h_duct) #4*A/perimeter
x_shock= 0.07
h_shock = nozzle_height(x_shock*1000)/1000 #[m]

p_a = 1.01325 #btm
gamma = 1.4
R=287.1 #J/kg/K
T0 = 300 #K

#b.1) Determine the  sonic conditions A_t = A*
T_star = T0/(gamma+1)*2
c_star = np.sqrt(gamma*R*T_star)

#b.2) determiner les conditions juste avant le shock
air_ratio1 = h_throat/h_shock #ratio between the throat and the first shock
M1 = useful.inv_air_ratio_sup(air_ratio1,np.array([1.5]))
T1 =  T0/(1+0.2*M1**2)
u1= M1*np.sqrt(1.4*R*T1)

#b.3) determiner la vitesse apres le shock en utilisant les shocks relations
M2 = np.sqrt((1+0.2*M1**2)/(1.4*M1**2-0.2))
T2 = T0/(1+0.2*M2**2)
u2 = M2*np.sqrt(1.4*R*T2)
h_star2 = h_shock*M2*(1.2/(1+0.2*M2**2))**(1.2/0.4)
print(h_star2,'here')
air_ratio2 = h_star2/h_duct

#3)  determiner les conditions the end of nozzle using air_ratio2
Mi =  useful.inv_air_ratio(air_ratio2,np.array([0.5]))
Ti = T0/(1+0.2*Mi**2)
ui = Mi * np.sqrt(gamma*R*Ti)
print('here ui ',ui,M2)
#b.3) we need to find pi but we have no information about the fanno flow after so lets find the conditions at the exit of the duct and then return to the nozzle
pe = 1.01325 #btm safe assumption pa=pe
# a.3.1) iterer sur Me

Me = 0.7 #premiere estimation
error=1
iter = 1
while error>0.0001 and iter<100 :
    Te = T0/(1+0.2*Me**2)
    ue = Me*np.sqrt(1.4*R*Te)
    rho_e = pe*10**(5)/(R*Te)

    Re = useful.Re(rho_e,ue,d,Te)
    lambda_duct = useful.perte_lambda(Re,0.05)
    a = useful.critical_lenghth(Mi)-0.5*lambda_duct*length_duct/Dh_duct

    Me_new = useful.inv_critical_length(a,np.array([0.5]))
    error =abs(Me_new-Me)
    iter +=1
    Me=Me_new
    print('Me',Me,'iter',iter, 'error',error)

#b.4) info supplemtaire au outlet
poe = pe*(1+0.2*Me**2)**(1.4/0.4)
poe_star = poe/(1/Me*((1+0.2*Me**2)/1.2)**(3)) # poe_star constant dans le duct (fanno flow)
p_star = pe*Me/(1.2/(1+0.2*Me**2))**(1.2/0.4)
print('conditions outlet duct: ',Me,ue,Te,pe,poe)

# a.5) resultats au inlet du duct
poi = poe_star*(1/Mi*((1+0.2*Mi**2)/1.2)**(3))
pi = poi/(1+0.2*Mi**2)**(1.4/0.4)
p_stari = pi*Mi/(1.2/(1+0.2*Mi**2))**(1.2/0.4)

rho_i = pi*10**(5)/(R*Ti)
print('p_star chequ',p_star,p_stari)
print('massflow chequ duct', rho_i*ui, rho_e*ue) #bizarre
print('conditions inlet duct', Mi,ui,Ti,pi,poi)

#b.6) conditions apres le shock
p2 = poi/(1+0.2*M2**2)**(1.4/0.4)
rho_2 = p2*10**(5)/(R*T2)
print('conditions after shock', M2,u2,T2,p2,poi)
#b.7) conditions avant le shock => shock relations
p1 = p2*1.2/(1.4*M1**2-0.2) #p1_2  =p2*(1+1.4*M2**2)/(1+1.4*M1**2)
po1 = p1*(T0/T1)**(1.4/0.4)
rho_1 = p1*10**(5)/(R*T1)
print('massflow chequ shock', rho_1*u1, rho_2*u2)
print('conditions before shock', M1,u1,T1,p1,po1)

#b.8) conditions at the throat
p_star = po1*(T_star/T0)**(1.4/0.4)
rho_star = p_star*10**(5)/(R*T_star)
print('conditions at throat',1.,c_star,T_star,p_star,po1)
print('reponse a la question b po = ',po1)

#b.7) massflow chequ
Qm = ui*rho_i*b*2*h_duct#kg/s
print('massflow', Qm)

#3.b faire un tres beau graphe

#3.b.1 nozzle
# on resepare en 2pour faire jusqu au throat puis apres jusqu au shock
x1 =np.linspace(0,28,1000)
ones=np.ones(len(x1))
h1= np.zeros(len(x1))
for i in range(len(h1)):
    h1[i]= nozzle_height(x1[i])/1000
air_ratios1 = h_throat/h1

Mx1_a = useful.inv_air_ratio(air_ratios1,0.5*np.ones(len(x1)))
px1_a = po1/((ones+0.2*Mx1_a**2)**(1.4/0.4))

x1b =np.linspace(28,30,500)
ones=np.ones(len(x1b))
h1b= np.zeros(len(x1b))
for i in range(len(h1b)):
    h1b[i]= nozzle_height(x1b[i])/1000
air_ratios1b = h_throat/h1b

Mx1_b = useful.inv_air_ratio(air_ratios1b,1.0*np.ones(len(x1b)))
px1_b = po1/((ones+0.2*Mx1_b**2)**(1.4/0.4))

x2 = np.linspace(30,32,500)
h2= np.zeros(len(x2))
for i in range(len(h2)):
    h2[i]= nozzle_height(x2[i])/1000
air_ratios2 = h_throat/h2

Mx1_2 = useful.inv_air_ratio_sup(air_ratios2,1.0*np.ones(len(x2)))
px1_2 = po1/((ones+0.2*Mx1_2**2)**(1.4/0.4))

x3 = np.linspace(32,70,1000)
h3= np.zeros(len(x3))
for i in range(len(h3)):
    h3[i]= nozzle_height(x3[i])/1000
air_ratios3 = h_throat/h3

Mx1_3 = useful.inv_air_ratio_sup(air_ratios3,1.5*np.ones(len(x3)))
px1_3 = po1/((1+0.2*Mx1_3**2)**(1.4/0.4))

x4 =np.linspace(70,120,1000)
h4= np.zeros(len(x4))
for i in range(len(h4)):
    h4[i]= nozzle_height(x4[i])/1000
air_ratios4 = h_star2/h4

Mx2 = useful.inv_air_ratio(air_ratios4,0.5*np.ones(len(x4)))
px2 = poi/((1+0.2*Mx2**2)**(1.4/0.4))

#3.b.1.plot
xa = np.linspace(0,70,1000)
xb = np.linspace(70,120,1000)
fig, ax = plt.subplots()
ax.set_title('b nozzle')
ax.plot(x1,Mx1_a,'-b',label='M')
ax.plot(x1b,Mx1_b,'-b')
ax.plot(x1,px1_a,'-r',label='p[bar]')
ax.plot(x1b,px1_b,'-r')
ax.plot(x2,Mx1_2,'-b')
ax.plot(x2,px1_2,'-r')
ax.plot(x3,Mx1_3,'-b')
ax.plot(x3,px1_3,'-r')
ax.plot(x4,Mx2,'-b')
ax.plot(x4,px2,'-r')
ax.plot(xb,poi*np.ones(len(xb)),'-c',label='p0[bar]')
ax.plot(xa,po1*np.ones(len(xa)),'-c')
ax.legend()
ax.grid()

#3.b.2 duct
x4=np.linspace(0,length_duct,1000)
ones=np.ones(len(x4))
a = -0.5*lambda_duct*x4/Dh_duct+ones*useful.critical_lenghth(Mi)
Mx_d = useful.inv_critical_length(a,0.2*ones)
pox = poe_star*(1/Mx_d*((ones+0.2*Mx_d**2)/1.2)**(1.2/0.4))
px = pox/(ones+0.2*Mx_d**2)**(1.4/0.4)

fig2=plt.figure()
ax2 =plt.subplot(311)
ax2.set_title('a duct')
ax2.plot(x4,Mx_d,'-b',label='M')
ax22=plt.subplot(312)
ax22.plot(x4,px,'-r',label='p[bar]')
ax23=plt.subplot(313)
ax23.plot(x4,pox,'-c',label='p0[bar]')
ax2.legend()
ax2.grid()
ax22.legend()
ax22.grid()
ax23.legend()
ax23.grid()

plt.show()
