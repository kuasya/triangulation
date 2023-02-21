from math import *
import matplotlib.pyplot as plt
from numpy import *

INTtime=[]
INTrate=[]
KFtime=[]
KFrate=[]

f=open('int_new_1.6.txt', 'r')
for line in f:
    words = line.split() # разбиение строки на слова
    INTtime.append(float(words[0]))
    INTrate.append(float(words[1]))
f.close()
print(INTrate[10],INTtime[10])
N=len(INTrate)
print(N,'INT')

#Построение
plt.step(INTtime,INTrate)
plt.title('int(1.6c)')
plt.xlabel(u'time')
plt.ylabel(u'count')
plt.grid()
plt.show()

N=len(INTrate)
print(N,'INT')
f=open('krf_1.6.txt', 'r')
for line in f:
    words = line.split() # разбиение строки на слова
    KFtime.append(float(words[0]))
    KFrate.append(float(words[1]))
f.close()

L=len(KFrate)
print(L, 'KF')

plt.step(KFtime,KFrate)
plt.title('krf(1.6c)')
plt.xlabel(u'time')
plt.ylabel(u'count')
plt.grid()
plt.show()

#Find background
BackINT=0
i=N//2
while i<N:
    BackINT=INTrate[i]+BackINT
    i=i+1
BackINT=BackINT/(N-N//2)

BackKF=0
i=L//2
while i<L:
    BackKF=KFrate[i]+BackKF
    i=i+1
BackKF=BackKF/(L-L//2)
print('Background:INT and KF')
print(BackINT,BackKF)

s=(27091-BackKF)/(411970.000-BackINT)

i=0
a=0
while i<L//2:
    a=((KFrate[i]-BackKF-s*(INTrate[i]-BackINT))**2)/(KFrate[i]+s*s*INTrate[i])+a
    i=i+1
j=0
k=0

while j<30:
    b=0
    i=0
    while i<L//2:
        b=((KFrate[i]-BackKF-s*(INTrate[i+j]-BackINT))**2)/(KFrate[i]+s*s*INTrate[i+j])+b
        i=i+1
    if b<a:
        a=b
        k=j
    j=j+1
print(s,k,KFtime[0]+71203-71203.024-INTtime[0+k],KFtime[0],INTtime[0+k])