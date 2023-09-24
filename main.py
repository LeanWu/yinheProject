import dynamics

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import poliastro.constants.general as constant

global mu
global Re
mu=constant.GM_earth.value
Re=constant.R_earth.value

# 获得速度
def getV(r,a):
    return math.sqrt(2*mu/r-mu/a)

# 获得周期
def getT(r):
    return 2*math.pi*math.sqrt(r**3/mu)

# 绘制轨道图
def drawOrb(rv):
    num=len(rv)
    x=np.zeros(num)
    y=np.zeros(num)
    z=np.zeros(num)
    for i in range(num):
        x[i]=rv[i,0]
        y[i]=rv[i,1]
        z[i]=rv[i,2]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.scatter(rv[0,0], rv[0,1], rv[0,2], c='red', s=1, label='Highlighted Point')
    ax.plot(x, y, z)
    plt.show()

# 变轨后新轨道
def get_newOrb(dv):
    h1=500e3
    v1=getV(Re+h1,Re+h1)
    T1=getT(Re+h1)
    v=v1-dv
    rv0=np.array([Re+h1,0,0,0,v,0])

    num=1000
    T=T1/2
    t=np.linspace(0,T,num)
    rv=dynamics.mypropagation(rv0,T,mu,T/(num-1))
    return t,rv


# 轨道高度
def geth(rv):
    h=np.zeros(len(rv))
    for i in range(len(rv)):
        h[i]=math.sqrt(rv[i,0]**2+rv[i,1]**2+rv[i,2]**2)-Re
    return h

# 得到120km时间
def getTime120(t,h):
    target_h=120e3
    interp_func = interp1d(h, t, kind='cubic')
    t_120=interp_func(target_h)
    return t_120


# 速度增量范围计算
h1=500e3
h2=120e3
v1=getV(Re+h1,Re+h1)
v2=getV(Re+h1,Re+(h1+h2)/2)
v3=getV(Re+h1,Re+(h1+0)/2)
dv_min=v1-v2
dv_max=v1-v3
print("v_min:",dv_min)
print("v_max:",dv_max)

# region
# 测试初始轨道
# T1=getT(Re+h1)
# num=1000
# rv0=np.array([Re+h1,0,0,0,v1,0])
# rv1=dynamics.mypropagation(rv0,T1,mu,T1/num)
# print(T1)
# drawOrb(rv1)
# endregion
# 单个速度增量测试
dv=dv_min
t,rv=get_newOrb(dv)
# drawOrb(rv)
h=geth(rv)
print("h_min",min(h))
t_120=getTime120(t,h)
print("using time:",t_120)
# plt.plot(t,h)
# plt.show()

# 投射物仿真绘图
num=100
dv_series=np.linspace(dv_min,dv_max,num)
t_120_series=np.zeros(num)
fig, ax = plt.subplots()
for i in range(num):
    t,rv=get_newOrb(dv_series[i])
    h=geth(rv)
    if i==0:
        name="dv_min:"+"{:.2f}".format(dv_series[i])+"m/s"
        ax.plot(t,h,color=(i/num,0.4,0.6),label=name)
    elif i==num-1:
        name="dv_max:"+"{:.2f}".format(dv_series[i])+"m/s"
        ax.plot(t,h,color=(i/num,0.4,0.6),label=name)
    else:
        ax.plot(t,h,color=(i/num,0.4,0.6))
    t_120_series[i]=getTime120(t,h)
plt.rcParams['font.sans-serif']=['SimHei']
ax.set_xlabel("time(s)")
ax.set_ylabel("height(m)")
ax.legend()
ax.set_title("不同速度增量下的投射物体高度变换曲线")
plt.show()
plt.plot(dv_series,t_120_series)
plt.xlabel("dv(m/s)")
plt.ylabel("time(s)")
plt.title("速度增量与投射物高度降低到120km所需时间的关系曲线")
plt.show()