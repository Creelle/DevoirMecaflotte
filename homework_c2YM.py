import numpy as np
from nozzle import nozzle_height
import functions as fun
import matplotlib.pyplot as plt

h_throat = 0.003 # 3 [mm]
h_duct = 0.0075 #7.5 [mm]
l_nozzle = 0.12 # [m]
l_duct = 1.5 # [m]
xe = l_duct # [m]
de = 2*h_duct
dt = 0.006 # 6 [mm]
b = 0.05 # 50 [mm]
Dh_duct = 4*h_duct*b/(b+2*h_duct) # 4*A/P
p_a = 101325.0 # [Pa] = 1 [atm]
T0 = 300 # [K]
# Gas = air (assuming air as a perfect gas)
gamma = 1.4
R = 287.1 # [J/kg/K]

p_e = p_a
x_sh = 0.12 # [m]
h_sh = nozzle_height(x_sh*1000)/1000 #[m]

# Throat conditions
T_t = 2*T0/(gamma+1)
T_star = T_t
c_t = np.sqrt(gamma*R*T_t)
u_t = c_t

# Shock relations :
#  -1- Isentropic between reservoir and sh1
#  -2- Isentropic between sh2 and end of the nozzle
M_sh1 = 1.5
M_sh1, error_M_sh1, iter_M_sh1 = fun.ratio_A_star_over_A_super(h_throat, h_sh, M_sh1)
print('M_sh1 = ', M_sh1)
T_sh1 = T0/(1 + 0.5*(gamma-1)*M_sh1**2)
c_1 = np.sqrt(gamma*R*T_sh1)
u_1 = M_sh1*c_1

M_sh2 = np.sqrt((1 + 0.5*(gamma-1)*M_sh1**2)/(gamma*M_sh1**2 - 0.5*(gamma-1)))
print('M_sh2 = ', M_sh2)
# T_sh2 = fun.shock_T_sh2(T_sh1, M_sh1)
T_sh2 = T0/(1 + 0.5*(gamma-1)*M_sh2**2)
print('T_sh2 = ', T_sh2)
c_2 = np.sqrt(gamma*R*T_sh2)
u_2 = M_sh2*c_2
h_star2 = h_sh * M_sh2 * (0.5*(gamma+1) / (1 + 0.5*(gamma-1)*M_sh2**2))**(0.5*(gamma+1)/(gamma-1))

M_i = M_sh2
M_e = 0.7
T_e = T0/(1 + 0.5*(gamma-1)*M_e**2)
c_e = np.sqrt(gamma*R*T_e)
u_e = c_e*M_e
rho_e = p_e/R/T_e
Re_de = fun.Re(rho_e, u_e, de, T_e)
print('Re_de = ', Re_de)
lambda_est = 1/(2.33*Re_de**0.1)**2
print('lambda_est = ', lambda_est)
lambda_e, error_inv_sqrt_lambda, iter_lambda_e = fun.lambd(Re_de, lambda_est)
print('lambda = ', lambda_e, 'error 1/sqrt(lambda) = ', error_inv_sqrt_lambda, 'iter lambda = ', iter_lambda_e)

# Iterating to find the right M_e
error = 1
iter = 1
while(iter<100 and error>0.0001):
    Me_new = 1 / np.sqrt(1 + (fun.f(M_i) - 0.5*lambda_e * xe/Dh_duct)*gamma - 0.5*(gamma+1)*np.log(1.2*M_e**2 / (1+0.2*M_e**2)))
    iter += 1
    print(iter)
    error = abs(Me_new - M_e)
    print(error)
    M_e = Me_new

    T_e = T0 / (1 + 0.5 * (gamma - 1) * M_e ** 2)
    c_e = np.sqrt(gamma * R * T_e)
    u_e = c_e * M_e
    rho_e = p_e / R / T_e
    Re_de = fun.Re(rho_e, u_e, de, T_e)
    lambda_e, error_inv_sqrt_lambda, iter_lambda_e = fun.lambd(Re_de, lambda_e)

print('M_e = ', M_e)
p_star_e = p_e*M_e*np.sqrt((1+0.5*(gamma-1)*M_e**2)/(0.5*(gamma+1)))
print('p_star_e = ', p_star_e*1e-5)
p_0e = p_e*(1+0.5*(gamma-1)*M_e**2)**(gamma/(gamma-1))
print('p_0e = ', p_0e*1e-5)
p_0star = p_0e*M_e*(0.5*(gamma+1)/(1 + 0.5*(gamma-1)*M_e**2))**(0.5*(gamma+1)/(gamma-1))
# p_0star = p_0e/(1/M_e*((1+0.2*M_e**2)/1.2)**(3))
print('p_0star = ', p_0star*1e-5)
p_0i = p_0star/M_i * ((1 + 0.5*(gamma-1)*M_i**2)/(0.5*(gamma+1)))**(0.5*(gamma+1)/(gamma-1))
print('p_0i = ', p_0i*1e-5)
P_0star_i =  p_0i*M_i*(0.5*(gamma+1)/(1 + 0.5*(gamma-1)*M_i**2))**(0.5*(gamma+1)/(gamma-1))
print('p_0star_i = ', P_0star_i*1e-5)
# rho_2 = u_e*rho_e/u_2
# p2 = rho_2*R*T_sh2*10**(-5)
# p_0i_R = p2*(T0/T_sh2)**(1.4/0.4)
# print('p_0i_R = ', p_0i_R)
p_i = p_0i/(1+0.5*(gamma-1)*M_i**2)**(gamma/(gamma-1))
print('p_i = ', p_i*1e-5)
p_star_i = p_i*M_i*np.sqrt((1+0.5*(gamma-1)*M_i**2)/(0.5*(gamma+1)))
print('p_star_i = ', p_star_i*1e-5)
p_sh2 = p_i
p_sh1 = p_sh2*(0.5*(gamma+1)/(gamma*M_sh1**2 - 0.5*(gamma-1)))
print('p_sh1 = ', p_sh1*1e-5)
p_t = p_sh1 * ((1 + 0.5*(gamma-1)*M_sh1**2)/(0.5*(gamma+1)))**(gamma/(gamma-1))
print('p_t = ', p_t*1e-5)
p_0t = p_t * (0.5*(gamma+1))**(gamma/(gamma-1))
print('p_0t = ', p_0t*1e-5)
p_0sh1 = p_sh1 * (1 + 0.5*(gamma-1)*M_sh1**2)**(gamma/(gamma-1))
print('p_0sh1 = ', p_0sh1*1e-5)
p_0sh1V2 = p_sh1*(T0/T_sh1)**(1.4/0.4)
print('p_0sh1_V2 = ', p_0sh1V2*1e-5)
# print(p_sh1*(T0/T_sh1)**(gamma/(gamma-1)) / 100000)
