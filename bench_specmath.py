# ~ from specmath import specmathclass as sm
import specmath
# ~ from specmath import timerclass as tc
import numpy as np
# ~ import math as mt
import time
# ~ f=sm()


class Timer:
    def __enter__(self):
        self.start = time.perf_counter()
        return self

    def __exit__(self, *args):
        self.end = time.perf_counter()
        self.interval = self.end - self.start
        print('Time left %.03f sec' % self.interval+'\n')


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


def bench_wigner(repeats=1000):
    print('wigner performance:')

    with Timer() as tc:
        for i in range(repeats):
            specmath.wigner(1, 0, 0, np.pi/4)

    with Timer() as tc:
        for i in range(repeats):
            specmath.wigner(1, -1, 0, np.pi/4)

def bench_diffwigner(repeats=1000):
    print('diffwigner performance:')

    with Timer() as tc:
        for i in range(repeats):
            specmath.diffwigner(1, 0, 0, np.pi/4)

    with Timer() as tc:
        for i in range(repeats):
            specmath.diffwigner(1, -1, 0, np.pi/4)

###Замер времени + расчет производной дэ-функции Вигнера

# ~ t=tc()
# ~ for i in range(1000):
    # ~ f.diffwigner(1,0,0,np.pi/4)
# ~ t=0

# ~ t=tc()
# ~ for i in range(1000):
    # ~ f.diffwigner(1,-1,0,np.pi/4)
# ~ t=0

# ~ print('\n')
# ~ print(f.diffwigner(1,1,0,np.pi/4))

# ~ print('\n')
# ~ print(f.diffwigner(1,0,0,np.pi/4))

# ~ print('\n')
# ~ print(f.diffwigner(1,-1,0,np.pi/4))

# ~ print('\n')
# ~ print(f.diffwigner(2,2,0,np.pi/4))

# ~ print('\n')
# ~ print(f.diffwigner(2,1,0,np.pi/4))

# ~ print('\n')
# ~ print(f.diffwigner(2,0,0,np.pi/4))

# ~ print('\n')
# ~ print(f.diffwigner(2,-1,0,np.pi/4))

# ~ print('\n')

if __name__ == '__main__':

    bench_wigner()
    bench_diffwigner()




