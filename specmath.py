import math as mt
# ~ import time
# ~ import sympy as sp
import numpy as np
from scipy import special


### best exported functions ###

def gammamn(m ,n):
    return gamma5(m, n)

def wigsmalld(j, m, m1, teta):
    return wigner(j, m, m1, teta)

def wigsmalld_diff(j, m, m1, teta):
    return diffwigner(j, m, m1, teta)


### auxilary functions ###

#Считает время выполнения фрагментов кода
# moved to bench_specmath.py
# ~ class timerclass:
    # ~ def __init__(self):
        # ~ self.startTime=time.time()
    # ~ def __del__(self):
        # ~ print(str((time.time()-self.startTime)*1E3)+' msec left')


#Считает факториал
def factorial(x):
        if (x==1 or x==0):
            return 1
        return x*factorial(x-1)


#Считаем коэфф гаммы в gamma1-6 разными способами
def gamma1(m,n):
    assert (type(m)==int and type(n)==int),'Type is wrong!'
    assert (np.fabs(m)<=n),'Absolute value of m is greater than n!'
    return np.sqrt(((2*n+1)*mt.factorial(n-m))/(4*np.pi*n*(n+1)*\
                                                mt.factorial(n+m)))

def gamma2(m,n):
    assert (type(m)==int and type(n)==int),'Type is wrong!'
    assert (np.fabs(m)<=n),'Absolute value of m is greater than n!'
    return np.sqrt(((2*n+1)*factorial(n-m))/(4*np.pi*n*(n+1)*factorial(n+m)))

def gamma3(m,n):
    assert (type(m)==int and type(n)==int),'Type is wrong!'
    assert (np.fabs(m)<=n),'Absolute value of m is greater than n!'
    q=(2*n+1)/(4*np.pi*n*(n+1))
    if m>=0:
        for i in range(n-m+1,n+m+1):
            q*=1/i
        return np.sqrt(q)
    else:
        for i in range(n+m+1,n-m+1):
            q*=i
        return np.sqrt(q)

def gamma4(m,n):
    assert (type(m)==int and type(n)==int),'Type is wrong!'
    assert (np.fabs(m)<=n),'Absolute value of m is greater than n!'
    q=(2*n+1)/(4*np.pi*n*(n+1))
    if m>=0:
        i=n-m+1
        while i!=n+m+1:
            q/=i
            i+=1
        return np.sqrt(q)
    else:
        i=n+m+1
        while i!=n-m+1:
            q*=i
            i+=1
        return np.sqrt(q)

def gamma5(m,n):
    assert (type(m)==int and type(n)==int),'Type is wrong!'
    assert (np.fabs(m)<=n),'Absolute value of m is greater than n!'
    gamma0= (2*n+1)/(4*n*np.pi*(n+1))
    if m>=0:
        i=1
        while i!=m+1:
            gamma0/=((n-i+1)*(n+i))
            i+=1
        return np.sqrt(gamma0)
    else:
        i=-1
        while i!=m-1:
            gamma0*=((n+i+1)*(n-i))
            i-=1
        return np.sqrt(gamma0)

def gamma6(m,n):
    assert (type(m)==int and type(n)==int), 'Type is wrong!'
    assert (np.abs(m)<=n), 'Absolute value of m is greater than n!'
    gamma0 = (2 * n + 1) / (4 * n * np.pi * (n + 1))
    if m>=0:
        def gamma_inner(i):
            if i==0:
                return gamma0
            return (1 / (n - i + 1) / (n+i)) * gamma_inner(i-1)
        return np.sqrt(gamma_inner(m))
    else:
        def gamma_inner(i):
            if i==0:
                return gamma0
            return (n - i) * (n + i + 1) * gamma_inner(i+1)
        return np.sqrt(gamma_inner(m))

#Определяю факториал на отрицательную область через гамма функцию Эйлера
# def gammafunc(p):
#         x=sp.symbols('x')
#         f=(x**(p-1))*sp.exp(-x)
#         integral=sp.integrate(f,(x,0,sp.oo))
#         return integral

#Считаю дэ-функцию Вигнера с ошибкой в пределах суммирования
# def wigner1(j,m,m1,teta):
#     assert (j>=m),'j must be greater than m or equal!'
#     kmin=min(0,m-m1)
#     kmax=max(0,j+m,j-m1)
#     cvar=np.sqrt(mt.factorial(j+m)*mt.factorial(j-m)*mt.factorial(j+m1)*\
#                                                       mt.factorial(j-m1))
#     Sum=0
#     k=kmin
#     if m==0:
#         while k!=kmax+1:
#             Sum+=(-1)**(k-m+m1)/(mt.factorial(j+m-k)*mt.factorial(j-k-m1)*\
#                                  mt.factorial(k-m+m1) * mt.factorial(k))*\
#                 (np.cos(teta/2))**(2*j-2*k+m-m1)*(np.sin(teta/2))**(2*k-m+m1)
#             k+=1
#     else:
#         while k!=kmax+1:
#             Sum+=(-1)**(k-m+m1)/(gammafunc(j+m-k+1)*gammafunc(j-k-m1+1)*\
#                                  gammafunc(k-m+m1+1)*gammafunc(k+1)) * \
#                 (np.cos(teta/2))**(2*j-2*k+m-m1)*(np.sin(teta/2))**(2*k-m+m1)
#             k+=1
#     return cvar*Sum

#Считаю дэ-функцию Вигнера с правильными пределами суммирования
def wigner(j, m, m1, teta):
    assert (j>=m),'j must be greater than m or equal!'
    kmin=max(0, m-m1)
    kmax=min(j+m, j-m1)
    cvar=np.sqrt(mt.factorial(j+m)*mt.factorial(j-m)*mt.factorial(j+m1)*\
                 mt.factorial(j-m1))
    Sum=0
    k=kmin
    while k != kmax + 1:
        Sum += (-1)**(k-m+m1) / (mt.factorial(j+m-k) * mt.factorial(j-k-m1)*\
                             mt.factorial(k-m+m1) * mt.factorial(k))*\
            (np.cos(teta/2))**(2*j-2*k+m-m1) * (np.sin(teta/2))**(2*k-m+m1)
        k += 1
    return cvar * Sum
                                                                               
#Считаю производную дэ-функции Вигнера
def diffwigner(j, m, m1, teta):
    assert (j >= m), 'j must be greater or equal than m'
    kmin = max(0, m-m1)
    kmax = min(j+m, j-m1)
    cvar=np.sqrt(mt.factorial(j+m) * mt.factorial(j-m) * mt.factorial(j+m1)*\
                 mt.factorial(j-m1))
    Sum = 0
    k = kmin
    while k != kmax + 1:
        Sum += ((-1)**(k-m+m1))*((2*j-2*k+m-m1)* \
              ((np.cos(teta/2))**(2*j-2*k+m-m1-1))*(-np.sin(teta/2)) * \
              ((np.sin(teta/2))**(2*k-m+m1)) + \
                  ((np.cos(teta/2))**(2*j-2*k+m-m1))*(2*k-m+m1)*\
                      ((np.sin(teta/2))**(2*k-m+m1-1))*np.cos(teta/2))/ \
            (2*mt.factorial(j+m-k)*mt.factorial(j-k-m1)*mt.factorial(k-m+m1)*\
             mt.factorial(k))
        k += 1
    return cvar * Sum

def bessel_spher_1kind(n,z):
    assert(n>=0 and z>0),'Something wrong in bessel_spher_1kind'
    return special.spherical_jn(n,z)

def bessel_spher_1kind_diff(n,z):
    assert(n>=0 and z>0),'Something wrong in bessel_spher_1kind_diff'
    return special.spherical_jn(n,z,derivative=True)

def bessel_spher_2kind(n,z):
    assert(n>=0 and z!=0),'Something wrong in bessel_spher_2kind'
    return special.spherical_yn(n,z)
    
def bessel_spher_2kind_diff(n,z):
    assert(n>=0 and z!=0),'Something wrong in bessel_spher_2kind_diff'
    return special.spherical_yn(n,z,derivative=True)

def hankel_spher(n, z):
    return bessel_spher_1kind(n, z) + 1j * bessel_spher_2kind(n, z)

def hankel_spher_diff(n, z):
    return bessel_spher_1kind_diff(n, z) + 1j * bessel_spher_2kind_diff(n, z)

def pi(m,n,teta):
    assert(not((teta%np.pi)==0)),'Sin in pi equals zero!'
    return m/np.sin(teta)/wigsmalld(n, m, 0, teta)

class vector_spher:
    def __init__(self,r=0,teta=0,phi=0):
        self.r=r
        self.teta=teta
        self.phi=phi
        
def M(m,n,x,teta,phi):
    return vector_spher(0,gammamn(m,n)*hankel_spher(n, x)*1j*pi(m,n,teta), \
                        -gammamn(m,n)*hankel_spher(n, x)* \
                        wigsmalld_diff(n, m, 0, phi))
        
def RgM(m,n,x,teta,phi):
    return vector_spher(0,gammamn(m,n)*bessel_spher_1kind(n, x)*1j* \
                        pi(m,n,teta),-gammamn(m,n)*bessel_spher_1kind(n, x)* \
                        wigsmalld_diff(n, m, 0, phi))
        
def N(m,n,x,teta,phi):
    return vector_spher(gammamn(m,n)*(n*(n+1)/x)*bessel_spher_1kind(n, x)* \
                        wigsmalld(n,m,0,teta),(1/x)* \
                            (hankel_spher(n, x)+x*hankel_spher_diff(n, x))* \
                                wigsmalld_diff(n, m, 0, teta),(1/x)* \
                                    (hankel_spher(n, x)+ \
                                     x*hankel_spher_diff(n, x))* \
                                        1j*pi(m,n,phi))   
        
def RgN(m,n,x,teta,phi):
    return vector_spher(gammamn(m,n)*(n*(n+1)/x)*hankel_spher(n, x)* \
                        wigsmalld(n,m,0,teta),(1/x)* \
                            (bessel_spher_1kind(n, x)+ \
                             x*bessel_spher_1kind_diff(n, x))* \
                                wigsmalld_diff(n, m, 0, teta),(1/x)* \
                                    (bessel_spher_1kind(n, x)+x* \
                                     bessel_spher_1kind_diff(n, x))* \
                                        1j*pi(m,n,phi))
    
# class vector:
#     def __init__(self,x=0,y=0,z=0):
#         self.x=x
#         self.y=y
#         self.z=z        
        
def fM(m,n,x,y,z):
    r=np.sqrt(x**2+y**2+z**2)
    teta=np.arccos(z/r)
    phi=np.arctan2(y,x)
    return gammamn(m,n)*hankel_spher(n, r)* \
                  (1j*pi(m,n,teta)*np.cos(teta)*np.cos(phi)+ \
                    wigsmalld_diff(n, m, 0, phi)*np.sin(phi)), \
                      gammamn(m,n)*hankel_spher(n, r)* \
                          (1j*pi(m,n,teta)*np.cos(teta)*np.sin(phi)- \
                            wigsmalld_diff(n, m, 0, phi)*np.cos(phi)), \
                              -gammamn(m,n)*hankel_spher(n, r)* \
                                  np.sin(teta)*1j* \
                                  pi(m,n,teta)

def fRgM(m,n,x,y,z):
    r=np.sqrt(x**2+y**2+z**2)
    teta=np.arccos(z/r)
    phi=np.arctan2(y,x)
    return gammamn(m,n)*bessel_spher_1kind(n, r)*\
                  (1j*pi(m,n,teta)*np.cos(teta)*np.cos(phi)+\
                    wigsmalld_diff(n, m, 0, phi)*np.sin(phi)),\
                      gammamn(m,n)*bessel_spher_1kind(n, r)*\
                  (1j*pi(m,n,teta)*np.cos(teta)*np.sin(phi)-\
                    wigsmalld_diff(n, m, 0, phi)*np.cos(phi)),\
                      -gammamn(m,n)*bessel_spher_1kind(n, r)*\
                          np.sin(teta)*1j*pi(m,n,teta)

def fN(m,n,x,y,z):
    r=np.sqrt(x**2+y**2+z**2)
    teta=np.arccos(z/r)
    phi=np.arctan2(y,x)
    assert(r!=0), 'In fN r==0'
    return gammamn(m,n)*(n*(n+1)/r)*bessel_spher_1kind(n, r)* \
                        wigsmalld(n,m,0,teta)*np.sin(teta)*np.cos(phi)+ \
                            (1/r)*(hankel_spher(n, r)+ \
                                    r*hankel_spher_diff(n, r))* \
                                wigsmalld_diff(n,m,0,teta)*np.cos(teta)* \
                                    np.cos(phi)-1j*np.sin(phi)*(1/r)* \
                                        (hankel_spher(n, r)+r* \
                                          hankel_spher_diff(n, r))* \
                                    pi(m,n,phi),gammamn(m,n)* \
                                    (n*(n+1)/r)*bessel_spher_1kind(n,r)* \
                        wigsmalld(n,m,0,teta)*np.sin(teta)*np.sin(phi)+ \
                        (1/r)*(hankel_spher(n,r)+\
                                r*hankel_spher_diff(n,r))* \
                                wigsmalld_diff(n,m,0,teta)*np.cos(teta)* \
                                    np.sin(phi)+1j*np.cos(phi)*(1/r)* \
                        (hankel_spher(n, r)+r*hankel_spher_diff(n, r))* \
                          pi(m,n,phi),gammamn(m,n)*(n*(n+1)/r)* \
                              bessel_spher_1kind(n, r)* \
                        wigsmalld(n,m,0,teta)*np.cos(teta)-\
                        (1/r)*(hankel_spher(n,r)+ \
                                r*hankel_spher_diff(n, r))* \
                                wigsmalld_diff(n, m, 0, teta)*np.sin(teta)

def fRgN(m,n,x,y,z):
    r=np.sqrt(x**2+y**2+z**2)
    teta=np.arccos(z/r)
    phi=np.arctan2(y,x)
    assert(r!=0), 'In fN r==0'
    return gammamn(m,n)*(n*(n+1)/r)*hankel_spher(n, r)* \
                        wigsmalld(n,m,0,teta)*np.sin(teta)*np.cos(phi)+ \
                            (1/r)*(bessel_spher_1kind(n, r)+ \
                                    r*bessel_spher_1kind_diff(n, r))* \
                                wigsmalld_diff(n,m,0,teta)*np.cos(teta)* \
                                    np.cos(phi)-1j*np.sin(phi)*(1/r)* \
                                        (bessel_spher_1kind(n, r)+ \
                                          r*bessel_spher_1kind_diff(n, r))* \
                                    pi(m,n,phi),gammamn(m,n)* \
                                    (n*(n+1)/r)*hankel_spher(n, r)* \
                        wigsmalld(n,m,0,teta)*np.sin(teta)*np.sin(phi)+ \
                        (1/r)*(bessel_spher_1kind(n, r)+ \
                                r*bessel_spher_1kind_diff(n, r))* \
                                wigsmalld_diff(n,m,0,teta)*np.cos(teta)* \
                                    np.sin(phi)+1j*np.cos(phi)*(1/r)* \
                        (bessel_spher_1kind(n, r)+ \
                          r*bessel_spher_1kind_diff(n, r))* \
                          pi(m,n,phi),gammamn(m,n)*(n*(n+1)/r)* \
                              hankel_spher(n, r)* \
                        wigsmalld(n,m,0,teta)*np.cos(teta)-\
                        (1/r)*(bessel_spher_1kind(n, r)+ \
                                r*bessel_spher_1kind_diff(n, r))* \
                                wigsmalld_diff(n, m, 0, teta)*np.sin(teta)        


















































































