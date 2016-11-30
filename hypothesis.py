# Hypothesis for vibrations project: shear thinning material
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import spline
from numpy import genfromtxt
import os
from math import sin, cos, pi, atan


plt.close('all')
# change the folder paths for the data
folderPath = 'C:\Users\owner\Desktop\Vibrations\Project\data\T114_50.csv'
T114_50 = genfromtxt(folderPath,skip_header=0, delimiter=',')
T114_50=np.array(T114_50)

folderPath = 'C:\Users\owner\Desktop\Vibrations\Project\data\T114_60.csv'
T114_60 = genfromtxt(folderPath,skip_header=0, delimiter=',')
T114_60=np.array(T114_60)

folderPath = 'C:\Users\owner\Desktop\Vibrations\Project\data\T114_72.csv'
T114_72 = genfromtxt(folderPath,skip_header=0, delimiter=',')
T114_72=np.array(T114_72)

folderPath = 'C:\Users\owner\Desktop\Vibrations\Project\data\T118_50.csv'
T118_50 = genfromtxt(folderPath,skip_header=0, delimiter=',')
T118_50=np.array(T118_50)

# Plot the G' and G" data

# fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()
# ax1.plot(T114_50[:,0],T114_50[:,5],'g-')
# ax2.plot(T114_50[:,0],T114_50[:,6],'b-')
# ax1.set_xlabel('Frequency (rad/s) constant strain')
# ax1.set_ylabel("G' (Pa)", color='g')
# ax2.set_ylabel("G'' (Pa)", color='b')
# plt.title("1:1:4 50%")
#
# fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()
# ax1.plot(T114_60[:,0],T114_60[:,5],'g-')
# ax2.plot(T114_60[:,0],T114_60[:,6],'b-')
# ax1.set_xlabel('Frequency (rad/s) constant strain')
# ax1.set_ylabel("G' (Pa)", color='g')
# ax2.set_ylabel("G'' (Pa)", color='b')
# plt.title("1:1:4 60%")
#
# fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()
# ax1.plot(T114_72[:,0],T114_72[:,5],'g-')
# ax2.plot(T114_72[:,0],T114_72[:,6],'b-')
# ax1.set_xlabel('Frequency (rad/s) constant strain')
# ax1.set_ylabel("G' (Pa)", color='g')
# ax2.set_ylabel("G'' (Pa)", color='b')
# plt.title("1:1:4 72%")
#
# fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()
# ax1.plot(T118_50[:,0],T118_50[:,5],'g-')
# ax2.plot(T118_50[:,0],T118_50[:,6],'b-')
# ax1.set_xlabel('Frequency (rad/s) constant strain')
# ax1.set_ylabel("G' (Pa)", color='g')
# ax2.set_ylabel("G'' (Pa)", color='b')
# plt.title("1:1:8 50%")


#data format:
# ang. frequency	temperature	time	osc. stress	  strain delta	G'	G''	ang. frequency
# rad/s	                C	      s	        Pa		 degrees         Pa	Pa	rad/s
# 0                     1          2        3           4            5   6    7

#assume rotary mass of hypothetical system
mHyp=.0000001 #kg
# not used v = .4   #poisson's ratio
# G' is the elastic shear modulus
# G" is the viscous shear modulus

# Wikipedia on stiffness gives k=GJ/L with J being the torsion constant
# J=pi*D^4/32 for a cylinder

L=.125*.0254*.01 #thickness of sample in meters
D=1.*.0254*1.3 # diameter of sample in meters
J=pi*D**4/32


VelnoFreq=.16*2.*pi/360.*D/2.#*freq

zetaCorrection = 1. # we need to figure out how to actually relate G' and G"
# Calculate the natural frequency and damping ratio

T4_50=[]
T4_60=[]
T4_72=[]
T8_50=[]
noDamp2=[]

wnT4_50=(T114_50[:,5]*J/L/mHyp)**.5
wnT4_60=(T114_60[:,5]*J/L/mHyp)**.5
wnT4_72=(T114_72[:,5]*J/L/mHyp)**.5
wnT8_50=(T118_50[:,5]*J/L/mHyp)**.5

# zT4_50=T114_50[:,6]*J/L/(2.*mHyp*wnT4_50)/zetaCorrection
# zT4_60=T114_60[:,6]*J/L/(2.*mHyp*wnT4_60)/zetaCorrection
# zT4_72=T114_72[:,6]*J/L/(2.*mHyp*wnT4_72)/zetaCorrection
# zT8_50=T118_50[:,6]*J/L/(2.*mHyp*wnT8_50)/zetaCorrection

zT4_50=T114_50[:,6]*J/L*(VelnoFreq*T114_50[:,0])/(2.*mHyp*wnT4_50)/zetaCorrection
zT4_60=T114_60[:,6]*J/L*(VelnoFreq*T114_60[:,0])/(2.*mHyp*wnT4_60)/zetaCorrection
zT4_72=T114_72[:,6]*J/L*(VelnoFreq*T114_72[:,0])/(2.*mHyp*wnT4_72)/zetaCorrection
zT8_50=T118_50[:,6]*J/L*(VelnoFreq*T118_50[:,0])/(2.*mHyp*wnT8_50)/zetaCorrection


#plot the wn and zeta

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(T114_50[:,0],wnT4_50,'g-')
ax2.plot(T114_50[:,0],zT4_50,'b-')
ax1.set_xlabel('Frequency (rad/s) constant strain')
ax1.set_ylabel("Wn (rad/s)", color='g')
ax2.set_ylabel("Zeta", color='b')
plt.title("1:1:4 50%")

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(T114_60[:,0],wnT4_60,'g-')
ax2.plot(T114_60[:,0],zT4_60,'b-')
ax1.set_xlabel('Frequency (rad/s) constant strain')
ax1.set_ylabel("Wn (rad/s)", color='g')
ax2.set_ylabel("Zeta", color='b')
plt.title("1:1:4 60%")

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(T114_72[:,0],wnT4_72,'g-')
ax2.plot(T114_72[:,0],zT4_72,'b-')
ax1.set_xlabel('Frequency (rad/s) constant strain')
ax1.set_ylabel("Wn (rad/s)", color='g')
ax2.set_ylabel("Zeta", color='b')
plt.title("1:1:4 72%")

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(T118_50[:,0],wnT8_50,'g-')
ax2.plot(T118_50[:,0],zT8_50,'b-')
ax1.set_xlabel('Frequency (rad/s) constant strain')
ax1.set_ylabel("Wn (rad/s)", color='g')
ax2.set_ylabel("Zeta", color='b')
plt.title("1:1:8 50%")

for i in range(len(T114_50[:,0])):
    a1=abs(1/(((1-(T114_50[i,0]**2/wnT4_50[i]**2)))**2+(2*zT4_50[i]*T114_50[i,0]/wnT4_50[i])**2)**.5)
    T4_50=np.append(T4_50,a1)

    a2=abs(1/(((1-(T114_60[i,0]**2/wnT4_60[i]**2)))**2+(2*zT4_60[i]*T114_60[i,0]/wnT4_60[i])**2)**.5)
    T4_60=np.append(T4_60,a2)

    a3=abs(1/(((1-(T114_72[i,0]**2/wnT4_72[i]**2)))**2+(2*zT4_72[i]*T114_72[i,0]/wnT4_72[i])**2)**.5)
    T4_72=np.append(T4_72,a3)

    a4=abs(1/(((1-(T118_50[i,0]**2/wnT8_50[i]**2)))**2+(2*zT8_50[i]*T118_50[i,0]/wnT8_50[i])**2)**.5)
    T8_50=np.append(T8_50,a4)

    # a5=abs(1/(1-(T118_50[i,0]**2/wnT4_50[1]**2)))
    # noDamp2=np.append(noDamp2,a5)
    a5=abs(1/(1-(T118_50[i,0]**2/600.**2)))
    noDamp2=np.append(noDamp2,a5)

plt.figure("one")
plt.plot(T114_50[:,0],T4_50,"r",label='T114_50')
plt.plot(T114_60[:,0],T4_50,"b", label = 'T114_60')
plt.plot(T114_72[:,0],T4_72,"k", label = 'T114_72')
plt.plot(T118_50[:,0],T8_50,"g", label = 'T118_50')
plt.plot(T118_50[:,0],noDamp2,"k+", label = 'Normal System without Damping')
plt.legend(loc = 'center left')
plt.ylim((0,4))
plt.title('Harmonically Forced Spring Damper')
plt.show()


#Original Hypothesis

# n=500
# endF=50
# m=.5
#
# w=np.linspace(0,endF,n)
# k=[]
#
# for i in range(n):
#
#     if w[i]<=20:
#
#         k1=1000
#     elif w[i]>20. and w[i]<40.:
#         k1=-48*w[i]+2122
#
#     else:
#         k1=213
#     k=np.append(k,k1)
#
# wn=np.divide(k,m)
# wn=wn**.5
# undamped = []
# for i in range(n):
#     a1=abs(1/(1-(w[i]**2/wn[i]**2)))
#     undamped=np.append(undamped,a1)
#
# speed=np.multiply(undamped,w)
# zetaThin=[]
# for i in range(n):
#     if speed[i] < 15.:
#
#         zetaThinning =1
#     else:
#         zetaThinning=14.904*speed[i]**-1.009
#
#     if zetaThinning<.05:
#         zetaThinning=.05
#
#
#     zetaThin=np.append(zetaThin,zetaThinning)
#
# damped = []
# dampedC = []
# dampedConstK=[]
# noDamp = []
# zetaC=.06
# wnC=(550.0/m)**.5
#
# for i in range(n):
#     a1=abs(1/(((1-(w[i]**2/wn[i]**2)))**2+(2*zetaThin[i]*w[i]/wn[i])**2)**.5)
#     damped=np.append(damped,a1)
#
#     a2=abs(1/(((1-(w[i]**2/wnC**2)))**2+(2*zetaC*w[i]/wnC)**2)**.5)
#     dampedC=np.append(dampedC,a2)
#
#     a3=abs(1/(((1-(w[i]**2/wnC**2)))**2+(2*zetaThin[i]*w[i]/wnC)**2)**.5)
#     dampedConstK=np.append(dampedConstK,a3)
#
#     a4=abs(1/(1-(w[i]**2/wnC**2)))
#     noDamp=np.append(noDamp,a4)
#
# xnew=np.linspace(0,endF,n/2-10)
# smoothDamped=spline(w,damped,xnew)
#
# plt.figure("two")
# plt.plot(xnew,smoothDamped,"r",label='Shear Thinning System (K and Zeta)')
# plt.plot(w,dampedConstK,"b", label = 'Shear Thinning System (Zeta Only)')
# plt.plot(w,dampedC,"k", label = 'Normal System with Damping')
# plt.plot(w,noDamp,"k+", label = 'Normal System without Damping')
# plt.legend(loc = 'center left')
# plt.ylim((0,12))
# plt.title('Harmonically Forced Spring Damper')
# #plt.show()
#
#
#

