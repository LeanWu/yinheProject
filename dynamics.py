import numpy as np
import math
from scipy.integrate import odeint

# 动力学方程
def dynamics(rv,t,mu):
    drdv=np.zeros(6)
    r_norm=np.linalg.norm(rv[:3])
    C=-mu/r_norm**3
    for i in range(6):
        if i<3:
            drdv[i]=rv[i+3]
        else:
            drdv[i]=C*rv[i-3]
    return drdv


# 积分求解
def mypropagation(rv0,dt,mu,t_step):
    if t_step==-1:
        num=2
    else:
        num=int(dt/t_step+1)    
    t=np.linspace(0,dt,num)
    new_dynamics=lambda rv,t:dynamics(rv,t,mu)
    rv=odeint(new_dynamics,rv0,t)
    return rv

# test：right rv:[ 0.86960342 -0.16077931 -0.05226604 -0.26656077 -0.69818329  0.70599091]
# rv0=np.array([0.6, 0.5, -0.6, 0.7, -0.5, 0.3])
# rv=mypropagation(rv0,1,1,1)
# print(rv)

# 加推力后的动力学方程
def dynamics_thrust(rv,t,mu,T,m):
    drdv=np.zeros(6)
    r_norm=np.linalg.norm(rv[:3])
    v_norm=np.linalg.norm(rv[3:6])
    C=-mu/r_norm**3
    for i in range(6):
        if i<3:
            drdv[i]=rv[i+3]
        else:
            drdv[i]=C*rv[i-3]+T/m*rv[i]/v_norm
    return drdv

# 推力积分求解
def mypropagation_thrust(rv0,dt,mu,t_step,T,m):
    if t_step==-1:
        num=2
    else:
        num=int(dt/t_step+1)    
    t=np.linspace(0,dt,num)
    new_dynamics_thrust=lambda rv,t:dynamics_thrust(rv,t,mu,T,m)
    rv=odeint(new_dynamics_thrust,rv0,t)
    return rv