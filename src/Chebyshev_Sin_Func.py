import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import pdb


def Z(x, x0, x1):
    return (2*x - x0 - x1)/(x1 - x0)


def cheb_T(n, x):
    # may want to make it return function T_n(x)
    if x == 1:
        return 1
    elif x == 0:
        return (-1)**n
    elif n == -1:
        return 0
    elif n == 0:
        return 1
    elif n == 1:
        return x
    else:
        return 2*x*cheb_T(n-1, x) - cheb_T(n-2, x)


def cheb_U(n, x):
    if x == 1:
        return n + 1
    elif x == 0:
        return (n + 1)*(-1)**n
    elif n == -1:
        return 0
    elif n == 0:
        return 1
    elif n == 1:
        return 2*x
    else:
        return 2*x*cheb_U(n-1, x) - cheb_U(n-2, x)


def V_hat(x, b, xmin, xmax):
    n = len(b)
    sumt = 0
    for i in range(n):
        sumt += b[i]*cheb_T(i, Z(x, xmin, xmax))
    return sumt


def dV_hat(x, b, xmin, xmax):
    n = len(b)
    sumt = 0
    for i in range(n):
        sumt += b[i]*i*cheb_U(i-1, Z(x, xmin, xmax))
    return sumt*2/(xmax - xmin)


def vander(x, N, xmin, xmax):
    V = [[0 for i in range(N)] for j in range(2*len(x))]
    for i in range(N):
        for j in range(len(x)):
            V[j][i] = cheb_T(i, Z(x[j], xmin, xmax))
        for k in range(len(x)):
            V[k + len(x)][i] = i*cheb_U(i-1, Z(x[k], xmin, xmax))*2/(xmax - xmin)
    return V


def vander1(x, N, xmin, xmax):
    V = [[0 for i in range(N)] for j in range(len(x))]
    for i in range(N):
        for j in range(len(x)):
            V[j][i] = cheb_T(i, Z(x[j], xmin, xmax))
    return V


def approx_objective(b, x, v, s):
    obj = 0
    N = len(x)
    xmin = np.min(x)
    xmax = np.max(x)
    for i in range(N):
        obj += (v[i] - V_hat(x[i], b, xmin, xmax))**2 + (s[i] - dV_hat(x[i], b, xmin, xmax))**2
        # (v[i] - V_hat(x[i], b, xmin, xmax))**2 + (s[i] - dV_hat(x[i], b, xmin, xmax))**2
        # (v[i] - V_hat(x[i], b, xmin, xmax))**2
        # (v[i] - V_hat(x[i], b, xmin, xmax) + (.5)*(s[i] - dV_hat(x[i], b, xmin, xmax)))**2
    return obj


def error(b, x, v, s):
    error = 0
    xx = np.linspace(np.min(x), np.max(x), len(x)*20)
    for i in range(len(xx)):
        error += (np.sin(xx[i]) - V_hat(xx[i], b, np.min(x), np.max(x)))**2
    return error

m = 4
approx_x = np.zeros(m)
z = np.zeros(m)
for i in range(m):
    z[i] = -np.cos((2*(i+1)-1)*np.pi/(2*m))
    approx_x[i] = (z[i]+1)*np.pi
approx_y = np.copy(approx_x)
for i in range(len(approx_x)):
    approx_y[i] = np.sin(approx_x[i])
# approx_y = unit(approx_x)
# approx_y = np.sin(approx_x)
approx_s = np.cos(approx_x)

# output = optimize.fmin(approx_objective, x0=[1]*degree,
#                        args=(approx_x, approx_y, approx_s), maxfun=10000)
V = vander(approx_x, 2*m, np.min(approx_x), np.max(approx_x))
V1 = vander1(approx_x, m+1, np.min(approx_x), np.max(approx_x))
b0 = np.linalg.lstsq(V, np.concatenate((approx_y, approx_s), axis=0))
b1 = np.linalg.lstsq(V1, approx_y.reshape(-1,1))

degree = m+1
b = np.zeros(degree)
b[0]=sum(approx_y)/m
for i in range(degree-1):
    tmp = 0
    for j in range(m):
        tmp += approx_y[j]*cheb_T(i+1, z[j])
    b[i+1] = 2*tmp/m

# print(output)
# print(error(output, approx_x, approx_y, approx_s) - error(b0[0], approx_x, approx_y, approx_s))
# dV_hat(approx_x[0], output, np.min(approx_x), np.max(approx_x))

xx = np.linspace(np.min(approx_x), np.max(approx_x), len(approx_x)*10)
yy = np.sin(xx)

approx_y1 = [0]*len(xx)
for j in range(len(xx)):
    approx_y1[j] = V_hat(xx[j], b1[0], np.min(approx_x), np.max(approx_x))

approx_s1 = [0]*len(xx)
for j in range(len(xx)):
    approx_s1[j] = V_hat(xx[j], b0[0], np.min(approx_x), np.max(approx_x))+.025

approx_s2 = [0]*len(xx)
# for j in range(len(xx)):
#     approx_s2[j] = V_hat(xx[j], output, np.min(approx_x), np.max(approx_x))

plt.ylim([-1.5, 1.5])
plt.plot(xx, yy, xx, approx_y1, xx, approx_s1, lw=2)
plt.plot(approx_x, approx_y, 'ro')
plt.legend(['Sin(x)', 'Lagrange w/ degree=5', 'Hermite w/ degree=9'], loc='best')
# plt.plot(xx, yy, xx, approx_yy, lw=2)
plt.savefig("../plot_Cheb_Sin_Func.png")
plt.close()
# plt.plot(xx,approx_y1, lw=2)
# plt.show()

# x0, x1 = 0, 2*np.pi
# X = np.linspace(x0, x1, 4)
# cheb_T(2)(.5)


# fig = plt.figure()
# ax = fig.gca(projection='3d')
# S = np.linspace(0, 2*np.pi, 100)
# X, Y = np.meshgrid(S, S)
# Z = np.sin(X) + np.sin(Y)
# ax.plot_wireframe(X, Y, Z, rstride=4, cstride=4, color='b')
# plt.show()
