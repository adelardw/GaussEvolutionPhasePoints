from celluloid import Camera
import math
import matplotlib.pyplot as plt
import numpy as np
# {start}
t0 = 0
tn = 2000
n = 750000
h = (tn-t0)/n
tk = np.arange(t0, tn + h, h)
# for the f(t)
c = 1
b = 0.1
a = 2
w = 0.008
t1 = 500
# for the Runge - Kutt (4) methods
def fr(t):
  return (1/(1+math.exp(-2*w*(t-t1))))+(1/(1+math.exp(2*w*(t-3*t1))))-1
def f1(p,x,t): # f1(p,x,t) -> [dx(t)/dt=f1(p,x,t)]
  return p
def f2(p,x,t):  # f2(p,x,t) -> [dp(t)/dt=f2(p,x,t)]
  return c*fr(t)-b*x-a*math.sin(x)
x0p0 = np.random.uniform(0, 1, (100, 2)) #2D - list <-> [x0_0,p0_0;x0_1,p0_1;...;x0_99,p0_99]
s_x = 0.5 # dispersion of x
s_p = 0.5 # dispersion of p
mu_x = 0 # average x
mu_p = 0 # average p
def Gauss_func(sx, sp, x_av, p_av, x, p): #2d GaussFunct(x, p)  == GaussFunct(x) * GaussFunct(p)
  return (1/((2*np.pi)*((sx*sp)**0.5)))*math.exp(-0.5*(((x-x_av)/sx)**2))*math.exp(-0.5*(((p-p_av)/sp)**2))
x = np.arange(-10,10,0.001)
y = [Gauss_func(0.5, 0.5, 0, 0, i, 0) for i in x]
# truncation method
x_start = mu_x-3*s_x
x_end = mu_x+3*s_x
p_start = mu_p-3*s_p
p_end = mu_p+3*s_p
z_start = 0
z_end = max(y)
zr = np.random.uniform(0, 1, (len(x0p0), 1))
xs = []
ps = []
for i in range(len(x0p0)):
  xone = x_start+(x_end-x_start)*x0p0[i][0]
  pone = p_start+(p_end-p_start)*x0p0[i][1]
  z1 = z_start+(z_end-z_start)*zr[i]
  if Gauss_func(s_x, s_p, mu_x, mu_p, xone, pone) > z1:
          i += 1
  else:
    xs.append(xone)
    ps.append(pone)
X_stepI = []
P_stepI = []
# we sorted every point from x0p0,and we have 2 lists:xs and ps <-> for i in xs (or ps):abs(xs[i])<=1(or abs(ps[i]))
# now used Runge-Kutt method for xs[i] and ps[i]
for i in range(len(xs)):
  un = ps[i]
  wn = xs[i]
  X_stepJ = [wn]
  P_stepJ = [un]
  for j in range(len(tk)-1):
    k1 = f2(un, wn, tk[j])*h
    m1 = f1(un, wn, tk[j])*h
    k2 = f2(un+m1/2, wn+k1/2, tk[j]+h/2)*h
    m2 = f1(un+m1/2, wn+k1/2, tk[j]+h/2)*h
    k3 = f2(un+m2/2, wn+k2/2, tk[j]+h/2)*h
    m3 = f1(un+m2/2, wn+k2/2, tk[j]+h/2)*h
    k4 = f2(un+m3, wn+k3, tk[j]+h/2)*h
    m4 = f1(un+m3, wn+k3, tk[j]+h/2)*h
    q1 = (k1+2*k2+2*k3+k4)/6
    q2 = (m1+2*m2+2*m3+m4)/6
    un = un+q1
    wn = wn+q2
    X_stepJ.append(wn)
    P_stepJ.append(un)
  X_stepI.append(X_stepJ)
  P_stepI.append(P_stepJ)
# if you want, you can delete # and see results (or you can check my HamiltomSys work)
# #plt.figure(2)
#plt.plot(X_stepI[0],P_stepI[0])
#plt.xlabel('x')
#plt.ylabel('p (x)')
# This is solution verification
#import scipy
#from scipy import integrate
#def f(y,t):
  #x=y[0]
  #p=y[1]
  #y11=p
  #y22=c*fr(t)-b*x-a*math.sin(x)
  #return [y11,y22]
#sol = integrate.odeint(f, [xs[0],ps[0]], tk)
#X_sol = []
#P_sol = []
#for i in range(len(sol)):
  #X_sol.append(sol[i][0])
  #P_sol.append(sol[i][1])
#plt.figure(3)
#plt.plot(X_sol,P_sol)
#plt.xlabel('x')
#plt.ylabel('p (x)')
# Lists_sorted
SortX = []
SortP = []
for k in range(n):
  Matrix_X = []
  Matrix_P = []
  for j in range(len(xs)):
    Matrix_X.append(X_stepI[j][k])
    Matrix_P.append(P_stepI[j][k])
  SortX.append(Matrix_X)
  SortP.append(Matrix_P)
# GaussEvolution
fig = plt.figure(4)
camera = Camera(fig)
plt.xlabel('x')
plt.ylabel('p (x)')
for k in range(0, n, 62500):
  for j in range(len(xs)):
    plt.scatter(SortX[k][j], SortP[k][j], s=50, color='red')
  camera.snap()
animation = camera.animate()
plt.show()
animation.save('image9.gif', writer = 'ffmpg')