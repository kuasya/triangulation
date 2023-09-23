from math import *
import matplotlib.pyplot as plt
from numpy import *

#Initial data
rw=1476375.925*(10**3)
tetaw=14.954
alfaw=40.731
alfaw=alfaw*pi/180
tetaw=tetaw*pi/180

rf=6923.206*(10**3)
tetaf=5.994
alfaf=293.391
alfaf=alfaf*pi/180
tetaf=tetaf*pi/180

ri=123687.7*(10**3)
tetai=64.697
alfai=10.623
alfai=alfai*pi/180
tetai=tetai*pi/180

rh=322916207.0*(10**3) #metres
tetah=-5.366
alfah=350.500
alfah=alfah*pi/180
tetah=tetah*pi/180

xf=rf*cos(alfaf)*cos(tetaf)  
yf=rf*sin(alfaf)*cos(tetaf)
zf=rf*sin(tetaf)

xw=rw*cos(alfaw)*cos(tetaw)
yw=rw*sin(alfaw)*cos(tetaw)
zw=rw*sin(tetaw)

xi=ri*cos(alfai)*cos(tetai)
yi=ri*sin(alfai)*cos(tetai)
zi=ri*sin(tetai)

xh=rh*cos(alfah)*cos(tetah)
yh=rh*sin(alfah)*cos(tetah)
zh=rh*sin(tetah)

#error for HEND
delta_rh=50.000
delta_alfah=0.001*pi/180
delta_tetah=0.001*pi/180

delta_xh=((delta_rh*cos(alfah)*cos(tetah))**2+(delta_alfah*sin(alfah)*rh*cos(tetah))**2+(delta_tetah*sin(tetah)*rh*cos(alfah))**2)**0.5
delta_yh=((delta_rh*sin(alfah)*cos(tetah))**2+(delta_alfah*cos(alfah)*rh*cos(tetah))**2+(delta_tetah*sin(tetah)*rh*sin(alfah))**2)**0.5
delta_zh=((delta_rh*sin(tetah))**2+(delta_tetah*cos(tetah)*rh)**2)**0.5

deltaTeta=zeros(6)
deltaTeta[0]=0.064 # in second
deltaTeta[1]=0.25
deltaTeta[2]=0.4
deltaTeta[3]=2
#deltaTeta[4]=1.6
deltaTeta[5]=8

c=299792458   #m/s
deltaT=zeros(6)
deltaT[0]=4.823 #s KW-KF
deltaT[1]=524.283 #s HEND-INTEGRAL
deltaT[2]=4.383   #KW-INTEGRAL
deltaT[3]=519.091 #HEND-KW
deltaT[4]=0.218 #KF-INTEGRAL
deltaT[5]=521.256  #HEND-KF

ModuleVec=zeros(6)
ModuleVec[0]=((xf-xw)**2+(yf-yw)**2+(zf-zw)**2)**0.5  # KW-KF
ModuleVec[1]=((xi-xh)**2+(yi-yh)**2+(zi-zh)**2)**0.5  # HEND-INTEGRAL
ModuleVec[2]=((xi-xw)**2+(yi-yw)**2+(zi-zw)**2)**0.5  # KW-INTEGRAL
ModuleVec[3]=((xw-xh)**2+(yw-yh)**2+(zw-zh)**2)**0.5  # HEND-KW
ModuleVec[4]=((xi-xf)**2+(yi-yf)**2+(zi-zf)**2)**0.5  # KF-INTEGRAL
ModuleVec[5]=((xf-xh)**2+(yf-yh)**2+(zf-zh)**2)**0.5  # HEND-KF
print(ModuleVec[0])
print(ModuleVec[1])

delta_MV=zeros(6) #delta_ModuleVec

delta_MV[1]=(((xi-xh)*delta_xh)**2+((yi-yh)*delta_yh)**2+((zi-zh)*delta_zh)**2)**0.5/(ModuleVec[1])
delta_MV[3]=(((xw-xh)*delta_xh)**2+((yw-yh)*delta_yh)**2+((zw-zh)*delta_zh)**2)**0.5/(ModuleVec[3])
delta_MV[5]=(((xf-xh)*delta_xh)**2+((yf-yh)*delta_yh)**2+((zf-zh)*delta_zh)**2)**0.5/(ModuleVec[5])
#print("delta_ModuleVec",delta_MV[1],delta_MV[3],delta_MV[5])

teta1=zeros(6)
alfa0=zeros(6)
teta0=zeros(6)

vector=zeros((6,3)) #запись разницы координат: на первом месте, куда раньше пришел сигнал
vector[0,0]=xw-xf
vector[0,1]=yw-yf
vector[0,2]=zw-zf

vector[1,0]=xh-xi
vector[1,1]=yh-yi
vector[1,2]=zh-zi

vector[2,0]=xw-xi
vector[2,1]=yw-yi
vector[2,2]=zw-zi

vector[3,0]=xh-xw
vector[3,1]=yh-yw
vector[3,2]=zh-zw

vector[4,0]=xf-xi
vector[4,1]=yf-yi
vector[4,2]=zf-zi

vector[5,0]=xh-xf
vector[5,1]=yh-yf
vector[5,2]=zh-zf


N=720
alfa=zeros(N)
i=1
while i<N:
    alfa[i]=alfa[i-1]+(2*pi)/(N-1)
    i=i+1

tetaP=zeros((6,N))
fiP=zeros((6,N))
xP=zeros(N)
yP=zeros(N)
tetaN=zeros((6,N))
fiN=zeros((6,N))
xN=zeros(N)
yN=zeros(N)
delta_teta1=zeros(6)
j=0
while j<6:
    teta1[j]=math.acos(c*deltaT[j]/ModuleVec[j])
    teta0[j]=math.asin(vector[j,2]/ModuleVec[j])
    if vector[j,1]>0:
        alfa0[j]=math.atan2(vector[j,1],vector[j,0])
    else: alfa0[j]=math.atan2(vector[j,1],vector[j,0])+2*pi

    delta_teta1[j]=((c/ModuleVec[j]*deltaTeta[j])**2+(c*deltaT[j]/(ModuleVec[j])**2*delta_MV[j])**2)**0.5/(1-(c*deltaT[j]/ModuleVec[j])**2)**0.5
    teta1[j]=teta1[j]*180/pi
    teta0[j]=teta0[j]*180/pi
    alfa0[j]=alfa0[j]*180/pi
    delta_teta1[j]=delta_teta1[j]*180/pi
    print(j)
    print('%s %s %s %s' %('teta0', 'alfa0','teta1','delta'))
    print('%f %f %f %f\n' %(teta0[j], alfa0[j],teta1[j], delta_teta1[j]))
    #Translate to radian
    teta0[j]=teta0[j]*pi/180
    alfa0[j]=alfa0[j]*pi/180    
    teta1[j]=teta1[j]*pi/180  
    delta_teta1[j]=delta_teta1[j]*pi/180

    i=0
    while i<N:
        xP[i]=-sin(alfa0[j])*sin(teta1[j]+delta_teta1[j])*sin(alfa[i])+cos(alfa0[j])*(sin(teta0[j])*cos(alfa[i])*sin(teta1[j]+delta_teta1[j])+cos(teta0[j])*cos(teta1[j]+delta_teta1[j]))
        yP[i]=cos(alfa0[j])*sin(teta1[j]+delta_teta1[j])*sin(alfa[i])+sin(alfa0[j])*(sin(teta0[j])*cos(alfa[i])*sin(teta1[j]+delta_teta1[j])+cos(teta0[j])*cos(teta1[j]+delta_teta1[j]))
        tetaP[j,i]=math.asin(sin(teta0[j])*cos(teta1[j]+delta_teta1[j])-cos(teta0[j])*sin(teta1[j]+delta_teta1[j])*cos(alfa[i]))
        if yP[i]>=0:
            fiP[j,i]=math.atan2(yP[i],xP[i])
        else: fiP[j,i]=math.atan2(yP[i],xP[i])+2*pi
        fiP[j,i]=fiP[j,i]*180/pi
        tetaP[j,i]=tetaP[j,i]*180/pi
        i=i+1
    i=0
    while i<N:
        xN[i]=-sin(alfa0[j])*sin(teta1[j]-delta_teta1[j])*sin(alfa[i])+cos(alfa0[j])*(sin(teta0[j])*cos(alfa[i])*sin(teta1[j]-delta_teta1[j])+cos(teta0[j])*cos(teta1[j]-delta_teta1[j]))
        yN[i]=cos(alfa0[j])*sin(teta1[j]-delta_teta1[j])*sin(alfa[i])+sin(alfa0[j])*(sin(teta0[j])*cos(alfa[i])*sin(teta1[j]-delta_teta1[j])+cos(teta0[j])*cos(teta1[j]-delta_teta1[j]))
        tetaN[j,i]=math.asin(sin(teta0[j])*cos(teta1[j]-delta_teta1[j])-cos(teta0[j])*sin(teta1[j]-delta_teta1[j])*cos(alfa[i]))
        if yN[i]>=0:
            fiN[j,i]=math.atan2(yN[i],xN[i])
        else: fiN[j,i]=math.atan2(yN[i],xN[i])+2*pi
        fiN[j,i]=fiN[j,i]*180/pi
        tetaN[j,i]=tetaN[j,i]*180/pi
        i=i+1
    
    j=j+1
#print(delta_teta1[0],delta_teta1[1],delta_teta1[2])
plt.scatter(fiP[0,:],tetaP[0,:], color='b', s=0.5)
plt.scatter(fiN[0,:],tetaN[0,:], color='b',s=0.5, label='KW-KRF')
plt.scatter(fiP[1,:],tetaP[1,:], color='r', s=0.5)
plt.scatter(fiN[1,:],tetaN[1,:], color='r',s=0.5, label='HEND-INTEGRAL')
plt.scatter(fiP[2,:],tetaP[2,:], color='g', s=1)
plt.scatter(fiN[2,:],tetaN[2,:], color='g',s=1, label='KW-INTEGRAL')
plt.scatter(fiP[3,:],tetaP[3,:], color='k', s=1)
plt.scatter(fiN[3,:],tetaN[3,:], color='k',s=1, label='HEND-KW')
#plt.scatter(fiP[4,:],tetaP[4,:], color='y', s=1)
#plt.scatter(fiN[4,:],tetaN[4,:], color='y',s=1, label='KF-INTEGRAL')
plt.scatter(fiP[5,:],tetaP[5,:], color='m', s=1)
plt.scatter(fiN[5,:],tetaN[5,:], color='m',s=1, label='HEND-KF')

plt.xlabel(u'alfa') 
plt.ylabel(u'teta') 
plt.grid()
plt.legend()
plt.show()
 
