#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 23:32:02 2017

@author: ashley
"""



{'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'}
"""

from astropy.io import fits
import math
import matplotlib.pyplot as plt
import numpy as np

"""
file=fits.open('/home/ashley/Link_xray/result/match_all')
#dec=file[1].data.field('dec_1')
redshift=file[1].data.field('z')
flux_b=file[1].data.field('flux_aper_b')
#print('flux_b:',flux_b)
#flux_w=file[1].data.field('flux_aper_w')
name=np.ndarray.tolist(file[1].data.field('name'))
#print('len_name_before',len(name))
name_1=list(set(name))
#print('len_name_after:',len(name_1))



D1=[i*3E8*308568E16/67.73999999999998 for i in redshift]#cm
D=[i*3.24077929e-25 for i in D1]#MPC

L=[]
for i in range(0,len(flux_b),1):
    l=4*math.pi*D1[i]*D1[i]*flux_b[i]
    L.append(l)
#log_L=[math.log10(i) for i in L]
#print('lem_flux_b',len(flux_b))
plt.scatter(D,L,facecolors='none',edgecolor='g')
plt.xlabel('Distance'+'  '+'(Mpc)')
plt.ylabel('L'+'  '+'(erg/s)')
plt.xscale('log')
plt.yscale('log')
plt.show()
plt.savefig('/home/ashley/Link_xray/result/match_allpic',format='png',dpi=200)
file.close()
"""
"""
file2=fits.open('/home/ashley/Link_xray/result/match_all')
#3911
redshift2=np.ndarray.tolist(file2[1].data.field('z'))
flux_b2=np.ndarray.tolist(file2[1].data.field('flux_aper_b'))
index=[i for i,v in enumerate(flux_b2) if v==0]
#print('index',index)
#print(file2[1].data[869][4])
#flux_w2=file2[1].data.field('flux_aper_w')
name2=np.ndarray.tolist(file2[1].data.field('name'))
acis_num=np.ndarray.tolist(file2[1].data.field('acis_num'))
#print('len_name2',len(name2))
#name_2=list(set(name2))
#print('len_name_22',len(name_2))


flux_all=[]
redshift_all=[]
name_all=[]
acis_num_all=[]
for i in range(0,len(flux_b2),1):
    if flux_b2[i]==0:
        i=i+1
    elif flux_b2[i]=='         ':
        
        i=i+1               
    else:
        flux_all.append(flux_b2[i])
        redshift_all.append(redshift2[i])
        name_all.append(name2[i]) 
        acis_num_all.append(acis_num[i])
#print('len(name_all)',len(name_all))
name_all1=list(set(name_all))
#print('len_name_all1',len(name_all1))

D22=[i*3E8*308568E16/67.73999999999998 for i in redshift_all]#cm
D2=[i*3.24077929e-25 for i in D22]

L2=[]
for i in range(0,len(flux_all),1):
    l2=4*math.pi*D22[i]*D22[i]*flux_all[i]
    L2.append(l2)
#print('len:',len(flux_all))
plt.figure()
plt.title('all')
plt.scatter(D2,L2,s=4,facecolors='none',edgecolor='y')
plt.xlabel('Distance'+' '+'Mpc')
plt.ylabel('L'+' '+'erg s-1')
plt.xscale('log')
plt.yscale('log')

plt.savefig('/home/ashley/Link_xray/result/pic_match',format='png',dpi=200)


#=================acis_time>=2 or <2
flux_1=[]
red_1=[]
name_1=[]
flux_2=[]
red_2=[]
name_2=[]
num_2=[]
for i in range(0,len(flux_b2),1):
    if flux_b2[i]==0:
        i=i+1
    elif flux_b2[i]=='         ':
        
        i=i+1               
    else:
        if acis_num[i]==1:
            flux_1.append(flux_b2[i])
            red_1.append(redshift2[i])
            name_1.append(name2[i])
        else:
            flux_2.append(flux_b2[i])
            red_2.append(redshift2[i])
            name_2.append(name2[i])
            num_2.append(acis_num[i])

#print('len(name_1)',len(name_1))
#print('len(name_2)',len(name_2))
D11=[i*3E8*308568E16/67.73999999999998 for i in red_1]
D111=[i*3.24077929e-25 for i in D11]
D21=[i*3E8*308568E16/67.73999999999998 for i in red_2]
D211=[i*3.24077929e-25 for i in D21]

L11=[]
for i in range(0,len(flux_1),1):
    l11=4*math.pi*D11[i]*D11[i]*flux_1[i]
    L11.append(l11)
plt.figure()
#plt.title('acis_num=1')
star=plt.scatter(D111,L11,s=4,facecolors='none',edgecolor='r',label='acis_num=1')
plt.legend(loc=0)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Distance'+' '+'Mpc')
plt.ylabel('L'+' '+'erg s-1')
plt.savefig('/home/ashley/Link_xray/result/pic_match1',format='png',dpi=200)
#===========================acis_num>=2
L22=[]
for i in range(0,len(flux_2),1):
    l22=4*math.pi*D21[i]*D21[i]*flux_2[i]
    L22.append(l22)
#plt.figure()
#plt.title('acis_num>1')
blue=plt.scatter(D211,L22,s=4,facecolors='none',edgecolor='b',label='acis_num>1')
plt.legend(loc=0)
plt.xscale('log')
plt.yscale('log')
#plt.xlabel('Distance'+' '+'Mpc')
#plt.ylabel('L'+' '+'erg s-1')
plt.savefig('/home/ashley/Link_xray/result/pic_match_all',format='png',dpi=200)
"""
    
