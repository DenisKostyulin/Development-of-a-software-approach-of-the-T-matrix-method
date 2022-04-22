# ~ from specmath import specmathclass as sm
import specmath
# ~ from specmath import timerclass as tc
import numpy as np
# ~ import math as mt

# ~ f=sm()

###Замер времени + проверка функций для расчета гаммы мин
#
# print('\ngamma1')
# print(np.fabs(f.gamma1(1,2))-105.13052*1E-3<1E-8)
# t=tc()
# for i in range(1000):
#     f.gamma1(1,2)
# t=0

# print('\ngamma2')
# print(np.fabs(f.gamma2(1,2))-105.13052*1E-3<1E-8)
# t=tc()
# for i in range(1000):
#     f.gamma2(1,2)
# t=0

# print('\ngamma3')
# print(np.fabs(f.gamma3(1,2))-105.13052*1E-3<1E-8)
# t=tc()
# for i in range(1000):
#     f.gamma3(1,2)
# t=0

# print('\ngamma4')
# print(np.fabs(f.gamma4(1,2))-105.13052*1E-3<1E-8)
# t=tc()
# for i in range(1000):
#     f.gamma4(1,2)
# t=0

# print('\ngamma5')
# print(np.fabs(f.gamma5(1,2))-105.13052*1E-3<1E-8)
# t=tc()
# for i in range(1000):
#     f.gamma5(1,2)
# t=0

# print('\ngamma6')
# print(np.fabs(f.gamma6(1,2))-105.13052*1E-3<1E-8)
# t=tc()
# for i in range(1000):
#     f.gamma6(1,2)
# t=0

###Замер времени + расчет по дэ-функции Вигнера

# t=tc()
# for i in range(1000):
#     f.wigner(1,0,0,np.pi/4)
# t=0

# t=tc()
# for i in range(1000):
#     f.wigner(1,-1,0,np.pi/4)
# t=0

def test_wigner():
    teta=np.pi/18
    assert np.isclose(specmath.wigner(1, 1, 0, teta),\
                      np.sin(teta)/np.sqrt(2))
    assert np.isclose(specmath.wigner(1, 0, 0, teta),\
                      np.cos(teta))
    assert np.isclose(specmath.wigner(1, -1, 0, teta),\
                      -np.sin(teta)/np.sqrt(2))
    assert np.isclose(specmath.wigner(2, 2, 0, teta),\
                      np.sqrt(3/8)*np.sin(teta)**2)
    assert np.isclose(specmath.wigner(2, 1, 0, teta),\
                      np.sqrt(3/8)*np.sin(2*teta))
    assert np.isclose(specmath.wigner(2, 0, 0, teta),\
                      (1/2)*((3*np.cos(teta)**2)-1))
    assert np.isclose(specmath.wigner(2, -1, 0, teta),\
                      -np.sqrt(3/8)*np.sin(2*teta))
    assert np.isclose(specmath.wigner(2, -2, 0, teta),\
                      np.sqrt(3/8)*np.sin(teta)**2)   

def test_diffwigner():
    teta=np.pi/18
    assert np.isclose(specmath.diffwigner(1, 1, 0, teta), (np.sqrt(2)/2)*\
                      np.cos(teta))
    assert np.isclose(specmath.diffwigner(1, 0, 0, teta), -np.sin(teta))
    assert np.isclose(specmath.diffwigner(1, -1, 0, teta), -(np.sqrt(2)/2)*\
                      np.cos(teta))
    assert np.isclose(specmath.diffwigner(2, 2, 0, teta), (np.sqrt(6))*\
                      np.sin(teta/2)*np.cos(teta/2)*np.cos(teta))
    assert np.isclose(specmath.diffwigner(2, 1, 0, teta), (-np.sqrt(6)/8)*\
                      ((6*np.sin(teta)**2)-np.cos(2*teta) - 3))
    assert np.isclose(specmath.diffwigner(2, 0, 0, teta), -6*np.sin(teta/2)*\
                      np.cos(teta/2)*np.cos(teta))
    assert np.isclose(specmath.diffwigner(2, -1, 0, teta), (np.sqrt(6)/8)*\
                      ((6*np.sin(teta)**2)-np.cos(2*teta) - 3))
    assert np.isclose(specmath.diffwigner(2, -2, 0, teta), (np.sqrt(6))*\
                      np.sin(teta/2)*np.cos(teta/2)*np.cos(teta))

def test_bessel_spher_1kind():
    z=10
    assert np.isclose(specmath.bessel_spher_1kind(0, z), np.sin(z)/z)
    assert np.isclose(specmath.bessel_spher_1kind(1, z), (np.sin(z)/(z**2))\
                      -np.cos(z)/z)
    assert np.isclose(specmath.bessel_spher_1kind(2, z), ((3/(z**3))-1/z)*\
                      np.sin(z) - 3*np.cos(z)/(z**2))
        
def test_bessel_spher_1kind_diff():
    z=10
    assert np.isclose(specmath.bessel_spher_1kind_diff(0, z),\
                      (z*np.cos(z)-np.sin(z))/(z**2))
    assert np.isclose(specmath.bessel_spher_1kind_diff(1, z),\
                      (2*z*np.cos(z)+np.sin(z)*(-2+z**2))/(z**3))
    assert np.isclose(specmath.bessel_spher_1kind_diff(2, z),\
                     ((9/z**3)-(1/z))*np.cos(z)+((4/z**2)-(9/z**4))*np.sin(z))
        
def test_bessel_spher_2kind():
    z=10
    assert np.isclose(specmath.bessel_spher_2kind(0, z), -np.cos(z)/z)
    assert np.isclose(specmath.bessel_spher_2kind(1, z), (-np.cos(z)/(z**2))\
                      -np.sin(z)/z)
    assert np.isclose(specmath.bessel_spher_2kind(2, z), -((3/(z**3))-1/z)*\
                      np.cos(z) - 3*np.sin(z)/(z**2))
        
def test_bessel_spher_2kind_diff():
    z=10
    assert np.isclose(specmath.bessel_spher_2kind_diff(0, z),\
                      (z*np.sin(z)+np.cos(z))/(z**2))
    assert np.isclose(specmath.bessel_spher_2kind_diff(1, z),\
                      (2*z*np.sin(z)-np.cos(z)*(-2+z**2))/(z**3))
    assert np.isclose(specmath.bessel_spher_2kind_diff(2, z),\
                     ((9/z**3)-(1/z))*np.sin(z)-((4/z**2)-(9/z**4))*np.cos(z))

##Замер времени + расчет производной дэ-функции Вигнера

# ~t=tc()
# ~for i in range(1000):
#      ~diffwigner(1,0,0,np.pi/4)
# ~t=0

# ~t=tc()
# ~for i in range(1000):
#     ~diffwigner(1,-1,0,np.pi/4)
# ~t=0

# ~print('\n')
# ~print(specmath.diffwigner(1,1,0,np.pi/4))

# ~print('\n')
# ~print(specmath.diffwigner(1,0,0,np.pi/4))

# ~print('\n')
# ~print(specmath.diffwigner(1,-1,0,np.pi/4))

# ~print('\n')
# ~print(specmath.diffwigner(2,2,0,np.pi/4))

# ~print('\n')
# ~print(specmath.diffwigner(2,1,0,np.pi/4))

# ~print('\n')
# ~print(specmath.diffwigner(2,0,0,np.pi/4))

# ~print('\n')
# ~print(specmath.diffwigner(2,-1,0,np.pi/4))

# ~print('\n')
# ~print(specmath.diffwigner(2,-2,0,np.pi/4))

# ~print('\n')

if __name__ == '__main__':
    test_wigner()
    test_diffwigner()
    test_bessel_spher_1kind()
    test_bessel_spher_1kind_diff()
    test_bessel_spher_2kind()
    test_bessel_spher_2kind_diff()
    # print(specmath.diffwigner(2, -2, 0, np.pi/4))

m=1
n=2
x=2
teta=np.pi/8
phi=np.pi/3    
print(specmath.fM(m,n,x,teta,phi)[0])
print(specmath.fM(m,n,x,teta,phi)[1])
print(specmath.fM(m,n,x,teta,phi)[2])
print(specmath.fRgM(m,n,x,teta,phi)[0])
print(specmath.fRgM(m,n,x,teta,phi)[1])
print(specmath.fRgM(m,n,x,teta,phi)[2])
print(specmath.fN(m,n,x,teta,phi)[0])
print(specmath.fN(m,n,x,teta,phi)[1])
print(specmath.fN(m,n,x,teta,phi)[2])
print(specmath.fRgN(m,n,x,teta,phi)[0])
print(specmath.fRgN(m,n,x,teta,phi)[1])
print(specmath.fRgN(m,n,x,teta,phi)[2])




