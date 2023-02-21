from math import *
import matplotlib.pyplot as plt
from numpy import *

KFtime=[]
KFrate=[]
KWtime=[]
KWrate=[]

f=open('krf.txt', 'r')
for line in f: 
    words = line.split() # разбиение строки на слова  
    KFtime.append(float(words[0]))
    KFrate.append(float(words[1]))
f.close()
print(KFrate[10],KFtime[10])


print(len(KFrate),'KRF')
i=0
while i<len(KFrate):
    KFtime[i]=KFtime[i]+71203
    i=i+1

f=open('kw.txt', 'r')
for line in f: 
    words = line.split() # разбиение строки на слова  
    KWtime.append(float(words[0]))
    KWrate.append(float(words[1]))
f.close()

L=len(KWrate)
print(len(KWrate), 'KW')
i=0
while i<len(KWrate):
    KWtime[i]=KWtime[i]+71198.539
    i=i+1
plt.step(KFtime,KFrate,color='b',label='Konus-RF')
plt.step(KWtime,KWrate,color='g',label='Konus-Wind')
plt.title('Временные истории для Konus-RF и Konus-Wind (0.064c)')
plt.xlabel(u't,с')
plt.ylabel(u'count')
plt.grid()
plt.legend()
plt.show()

#Find background
BackKF=0
i=len(KFrate)//2
while i<len(KFrate):
    BackKF=KFrate[i]+BackKF
    i=i+1
BackKF=BackKF/(len(KFrate)-len(KFrate)//2)

BackKW=0
i=len(KWrate)//2
while i<len(KWrate):
    BackKW=KWrate[i]+BackKW
    i=i+1
BackKW=BackKW/(len(KWrate)-len(KWrate)//2)
print('Background:KRF and KW')
print(BackKF,BackKW)

s=(7101.85-BackKW)/(2832-BackKF)

i=0   #выбирать самостоятельно исходя из картинки
a=0
while i<len(KWrate)//2:#выбирать самостоятельно исходя из картинки
    a=((KWrate[i]-BackKW-s*(KFrate[i]-BackKF))**2)/(KWrate[i]+s*s*KFrate[i])+a
    i=i+1
j=0
k=0
N=150      #выбирать самостоятельно исходя из картинки
err=zeros(N)
print(a/(len(KWrate)//2))
while j<N:
    b=0
    i=0
    while i<len(KWrate)//2:
        b=((KWrate[i]-BackKW-s*(KFrate[i+j]-BackKF))**2)/(KWrate[i]+s*s*KFrate[i+j])+b
        i=i+1
    #print(b/(L//2))
    err[j]=b/(len(KWrate)//2)
    if b<a:
        a=b
        k=j
    j=j+1
          

print(s,a/(len(KWrate)//2),k,KWtime[0]-KFtime[0+k],KFtime[0+k],KWtime[0])
#print(KWtime[0]+71198.539-71203-KFtime[0+k+int((2*(L//2))**0.5)]-(KWtime[0]+71198.539-71203-KFtime[0+k]))
#print(KWtime[0]+71198.539-71203-KFtime[0+k-int((2*(L//2))**0.5)]-(KWtime[0]+71198.539-71203-KFtime[0+k]))
#print(L//2,int((2*(L//2))**0.5))
j=zeros(N)
i=0
while i<N:
    j[i]=i
    i=i+1
plt.step(j,err)
plt.grid()
plt.show()

