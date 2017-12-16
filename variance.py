#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 10:37:01 2017

@author: ashley
"""

#==================this code is used to calculate the variance of one galaxy and then draw the diatance_variance pic.
#==================the asic_num>=2 galaxy can calculate the variance
#======nan
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import pandas as pa
from matplotlib.patches import Rectangle
import matplotlib.patches as pt
import math
"""
file read
"""
#file
file=fits.open('/media/ashley/KINGSTON/resultge2.fits')#match all the data from SDSS and CXO
name=np.ndarray.tolist(file[1].data.field('name_1'))
redshift=np.ndarray.tolist(file[1].data.field('z'))
acis_num=np.ndarray.tolist(file[1].data.field('acis_num_1'))
#file1
file1=fits.open('/media/ashley/KINGSTON/dataforyang2.fits')#3770=acis_num>=1   830 acis_num>=2
flux=np.ndarray.tolist(file1[1].data.field('flux_aper_b'))
flux_hilim=np.ndarray.tolist(file1[1].data.field('flux_aper_hilim_b'))
flux_lolim=np.ndarray.tolist(file1[1].data.field('flux_aper_lolim_b'))
name1=np.ndarray.tolist(file1[1].data.field('name'))
acis_num1=np.ndarray.tolist(file1[1].data.field('acis_num'))


#=============pick up the name of the acis_num >=2
#=============svar >3 xianzhu guangbian 
name_2=[]
redshift_2=[]
acis_2=[]
for i in range(0,len(name),1):
    if acis_num[i] >= 2:#830
        #print(i)
        name_2.append(name[i])
        redshift_2.append(redshift[i])
        acis_2.append(acis_num[i])
    else:
        i+=1

D=[]
fluxratio=[]  
variance=[]
mean_L=[]
max_L=[]
min_L=[]
excess_var=[]
id_fluxeq0=[]
Sd2=[]
svar=[]
#print('len(name_3)',len(name_3))

for i in range(0,len(name_2),1):
    indexname=name1.index(name_2[i])
    acisnum=acis_2[i]
    redshift_tmp=redshift[i]
    
    dist=redshift_tmp*3E8*308568E16/70#cm
    dist_mpc=dist*3.24077929e-25
    D.append(dist_mpc)
    tmp_flux=[]
    #mp_flux_h=[]
    uncertainty_flux=[]
    tmpstd_flux=[]
    std=[]
    if math.isnan(flux[indexname])==1 and math.isnan(flux_lolim[indexname])==1 and math.isnan(flux_hilim[indexname])==1:
        while math.isnan(flux[indexname])==1:
            
            indexname+=1
    
    for j in range(indexname,indexname + acisnum,1):
        
            
        if flux[j]==0.0:
            id_fluxeq0.append(j)
            tmp_flux.append(flux_hilim[j])
            tmpstd_flux.append(flux_hilim[j])
            tmpstd_flux.append(0.0)
            std_tmp=np.std(tmpstd_flux)
            std.append(std_tmp)
        #    uncertainty=flux_hilim[j]/2
       #     uncertainty_flux.append(uncertainty)
         #   j=j+1
        elif flux[j]!=0.0 and math.isnan(flux_lolim[j])==1:
            tmp_flux.append(flux[j])
            uncertainty=flux_hilim[j]-flux[j]
            tmpstd_flux.append(flux_hilim[j])
            tmpstd_flux.append(flux[j])
            std_tmp=np.std(tmpstd_flux)
            std.append(std_tmp)
            #uncertainty_flux.append(uncertainty)
        else:
            tmp_flux.append(flux[j])
            if flux_hilim[j]-flux[j]<=flux[j]-flux_lolim[j]:
                uncertainty=flux_hilim[j]-flux[j]
                uncertainty_flux.append(uncertainty)
            else:
                uncertainty=flux[j]-flux_lolim[j]
                uncertainty_flux.append(uncertainty)
            tmpstd_flux.append(flux_hilim[j])
            tmpstd_flux.append(flux[j])
            tmpstd_flux.append(flux_lolim[j])
            std_tmp=np.std(tmpstd_flux)
            std.append(std_tmp)
    #print('i',i)
    #print('j',j)
    max_flux=np.max(tmp_flux)
    index_max=tmp_flux.index(max_flux)
    min_flux=np.min(tmp_flux)
    index_min=tmp_flux.index(min_flux)
    err_min=std[index_min]
    err_max=std[index_max] #Fmax and Fmin are the maximum and minimum EPIC fluxes of a unique source, with the corresponding (statistica
    svar_tmp=(max_flux-min_flux)/np.sqrt(err_min**2+err_max**2)
    svar.append(svar_tmp)
    
    ratio=max_flux/min_flux
    fluxratio.append(ratio)
    mean_flux=np.mean(tmp_flux)

    
    max_l=max_flux*dist**2*4*np.pi
    mean_l=mean_flux*dist**2*4*np.pi
    min_l=min_flux*dist**2*4*np.pi
    
    mean_L.append(mean_l)
    log_meanl=[np.log10(n) for n in mean_L]
    max_L.append(max_l)
    min_L.append(min_l)
    tmp_excessvar=0
    for n in range(0,len(tmp_flux),1):
        tmp_excessvar+=((tmp_flux[n]-mean_flux)**2-(std[n])**2)/((len(tmp_flux))*(mean_flux)**2)
    excess_var.append(tmp_excessvar)# 
    tmp_sd1=0
    for n in range(0,len(tmp_flux),1):
        tmp_sd1+=((tmp_flux[n]-mean_flux)**2-(std[n])**2-tmp_excessvar)**2*mean_flux**2/(len(tmp_flux)-1)
    tmp_sd2=tmp_sd1/(mean_flux**2*np.sqrt(len(tmp_flux)))
    Sd2.append(tmp_sd2)
svargt3=[n for n in svar if n >=3] 

print('svar',svargt3)   
Sd=[np.sqrt(n) for n in Sd2]    
meanl_div1040=[n/(10**40) for n in mean_L]
sqrt_rms=[np.sqrt(n) for n in excess_var]
    #print(excess_var)
logmeanl=[np.log10(n) for n in mean_L]
#print('x',log_meanl)
#print('y',excess_var)
#=====change y to square
sigma1=[];sigma2=[];sigma3=[];sigma4=[];sigma5=[];sigma6=[];sigma7=[];sigma8=[];sigma9=[]
sigma10=[];sigma11=[]
meanl1=[];meanl2=[];meanl3=[];meanl4=[];meanl5=[];meanl6=[];meanl7=[];meanl8=[];meanl9=[];
meanl10=[];meanl11=[]
evar1=[];evar2=[];evar3=[];evar4=[];evar5=[];evar6=[]
for n in range(0,len(logmeanl),1):
    if 40>=logmeanl[n]>=39 :
        sigma1.append(Sd2[n])
        meanl1.append(logmeanl[n])
        evar1.append(excess_var[n])
    #if 39.5<=logmeanl[n]<=40 :
     #   sigma2.append(Sd2[n])
      #  meanl2.append(logmeanl[n])
    if 41>=logmeanl[n]>=40 :
        sigma2.append(Sd2[n])
        meanl2.append(logmeanl[n])
        evar2.append(excess_var[n])
    #if 40.5<=logmeanl[n]>=41 :
     #   sigma4.append(Sd2[n])
      #  meanl4.append(logmeanl[n])
    if 42>=logmeanl[n]>=41 :
        sigma3.append(Sd2[n])
        meanl3.append(logmeanl[n])
        evar3.append(excess_var[n])
    #if 41.5<=logmeanl[n]>=42 :
     #   sigma6.append(Sd2[n])
      #  meanl6.append(logmeanl[n])
    if 43>=logmeanl[n]>=42 :
        sigma4.append(Sd2[n])
        meanl4.append(logmeanl[n])
        evar4.append(excess_var[n])
    #if 42.5<=logmeanl[n]>=43 :
     #   sigma8.append(Sd2[n])
      #  meanl8.append(logmeanl[n])
    if 44>=logmeanl[n]>=43 :
        sigma5.append(Sd2[n])
        meanl5.append(logmeanl[n])
        evar5.append(excess_var[n])
    #if 44>=logmeanl[n]>=43.5 :
     #   sigma10.append(Sd2[n])
      #meanl10.append(logmeanl[n])
    if 45>=logmeanl[n]>=44 :
        sigma6.append(Sd2[n])
        meanl6.append(logmeanl[n])
        evar6.append(excess_var[n])
mean_sigma1=np.sum(sigma1)/len(sigma1)
mean_evar1=np.sum(evar1)/len(evar1)
tmp_err1=0;tmp_err2=0;tmp_err3=0;tmp_err4=0;tmp_err5=0;tmp_err6=0;
#tmp_err7=[];tmp_err8=[];tmp_err9=[];tmp_err10=[];tmp_err11=[]
#mean_ml1=[];mean_ml2=[];mean_ml3=[];mean_ml3=[];mean_ml4=[];mean_ml5=[];
#mean_ml6=[];mean_ml7=[];mean_ml8=[];mean_ml9=[];mean_ml10=[];mean_ml11=[]
for n in range(0,len(evar1),1):
    tmp_err1+=(evar1[n]-mean_evar1)**2/(len(evar1)*(len(evar1)-1))
print('tmp_err1',tmp_err1)
sqr_err1=np.sqrt(tmp_err1)
print('sqr_err1',sqr_err1)
mean_ml1=np.mean(meanl1)

mean_sigma2=np.sum(sigma2)/len(sigma2)
mean_evar2=np.sum(evar2)/len(evar2)
for n in range(0,len(evar2),1):
    tmp_err2+=(evar2[n]-mean_evar2)**2/(len(evar2)*(len(evar2)-1))
sqr_err2=np.sqrt(tmp_err2)
mean_ml2=np.mean(meanl2)

mean_sigma3=np.sum(sigma3)/len(sigma3)
mean_evar3=np.sum(evar3)/len(evar3)
for n in range(0,len(evar3),1):
    tmp_err3+=(evar3[n]-mean_evar3)**2/(len(evar3)*(len(evar3)-1))
sqr_err3=np.sqrt( tmp_err3)
mean_ml3=np.mean(meanl3)

mean_sigma4=np.sum(sigma4)/len(sigma4)
mean_evar4=np.sum(evar4)/len(evar4)
for n in range(0,len(evar4),1):
    tmp_err4+=(evar4[n]-mean_evar4)**2/(len(evar4)*(len(evar4)-1))
sqr_err4=np.sqrt(tmp_err4) 
mean_ml4=np.mean(meanl4)

mean_sigma5=np.sum(sigma5)/len(sigma5)
mean_evar5=np.sum(evar5)/len(evar5)
for n in range(0,len(evar5),1):
    tmp_err5+=(evar5[n]-mean_evar5)**2/(len(evar5)*(len(evar5)-1))
sqr_err5=np.sqrt(tmp_err5) 
mean_ml5=np.mean(meanl5)

mean_sigma6=np.sum(sigma6)/len(sigma6)
mean_evar6=np.sum(evar6)/len(evar6)
for n in range(0,len(evar6),1):
    tmp_err6+=(evar6[n]-mean_evar6)**2/(len(evar6)*(len(evar6)-1))
sqr_err6=np.sqrt(tmp_err6)
mean_ml6=np.mean(meanl6)

xlolim=[40-mean_ml1,41-mean_ml2,42-mean_ml3,43-mean_ml4,44-mean_ml5,45-mean_ml6]#need to change 
#print(len(xuplim))
xuplim=[mean_ml1-39,mean_ml2-40,mean_ml3-41,mean_ml4-42,mean_ml5-43,mean_ml6-44]
yuplim=[sqr_err1,sqr_err2,sqr_err3,sqr_err4,sqr_err5,sqr_err6]
#sqr_err7,sqr_err8,sqr_err9,sqr_err10,sqr_err11
ylolim=[sqr_err1,sqr_err2,sqr_err3,sqr_err4,sqr_err5,sqr_err6]

x=[mean_ml1,mean_ml2,mean_ml3,mean_ml4,mean_ml5,mean_ml6]
y=[mean_evar1,mean_evar2,mean_evar3,mean_evar4,mean_evar5,mean_evar6]
  

B=pt.ArrowStyle("Fancy", head_length=.4, head_width=.4, tail_width=.4)
a=plt.scatter(logmeanl,excess_var,s=4,facecolors='none',edgecolor='grey',label='acis_num>=2')

plt.legend(loc=0)
#plt.yscale('log')
#plt.errorbar(logmeanl,excess_var,ms=4,yerr=Sd2,lolims=True,uplims=False,markerfacecolor='none',markeredgecolor='g',fmt='o',c='grey')
#plt.errorbar(logmeanl,excess_var,ms=4,yerr=Sd2,lolims=False,uplims=True,markerfacecolor='none',markeredgecolor='g',fmt='o',c='grey')
plt.scatter(x,y,s=4,facecolor='none',edgecolor='r')
plt.errorbar(x,y,ms=4,yerr=ylolim,lolims=True,xlolims=True,xerr=xlolim,markerfacecolor='none',markeredgecolor='r',fmt='o',c='r')
plt.errorbar(x,y,ms=4,uplims=True,yerr=yuplim,xuplims=True,xerr=xuplim,markerfacecolor='none',markeredgecolor='r',fmt='o',c='r')

plt.xlabel('Mean '+'$logL_{x}$'+' '+'   '+'(0.5-7keV)'+' '+'(erg $s^{-1}$)')
plt.ylabel('$\sigma_{rms}^2$')

plt.ylim(-0.1,0.3)
plt.savefig('/media/ashley/KINGSTON/excess_varge2.pdf',dpi=400,format='pdf')  
     

ratiogt10=[]
ratiogt5=[]
ratiolt5=[]
minlgt10=[]
minlgt5=[]
minllt5=[]
for i in range(0,len(fluxratio),1):
    if fluxratio[i]>=10:
        ratiogt10.append(fluxratio[i])
        minlgt10.append(min_L[i])
    if fluxratio[i]>=5:
        ratiogt5.append(fluxratio[i])
        minlgt5.append(min_L[i])
    if fluxratio[i]<5:
        ratiolt5.append(fluxratio[i])
        minllt5.append(min_L[i])
        
logminl=[np.log10(n) for n in min_L]
logminlgt5=[np.log10(n) for n in minlgt5]
logminllt5=[np.log10(n) for n in minllt5]
sminl=pa.Series(logminl)
sminlgt5=pa.Series(logminlgt5)
sminllt5=pa.Series(logminllt5)

log_minlgt10=[np.log10(n) for n in minlgt10]
plt.figure()
Sminl=sminl.dropna()

Sminlgt5=sminlgt5.dropna()
Sminllt5=sminllt5.dropna()
bin1=[39,39.2,39.4,39.6,39.8,40,40.2,40.4,40.6,40.8,41, 41.2, 41.4, 41.6, 41.8,42, 42.2, 42.4, 42.6, 42.8,43, 43.2, 43.4, 43.6, 43.8,44, 44.2, 44.4, 44.6, 44.8]

plt.hist(Sminl,bins=bin1,facecolor='none',edgecolor='g')
plt.hist(log_minlgt10,bins=bin1,facecolor='none',edgecolor='b') 
handels=[Rectangle((0,0),1,1,color=c) for c in ['g','b']]
labels=["all","$F_{var}$>10"]
plt.legend(handels,labels)   
plt.xlabel('log_minL')
plt.title('$F_{var}$')
plt.yscale('log')
plt.savefig('/media/ashley/KINGSTON/histall.pdf',dpi=200,format='pdf')

plt.figure()
plt.hist(Sminlgt5,bins=bin1,facecolor='none',edgecolor='r')
plt.hist(Sminllt5,bins=bin1,facecolor='none',edgecolor='y')
handels1=[Rectangle((0,0),1,1,color=c) for c in ['r','y']]
labels1=["$F_{var} $>=5","$F_{var}$<5"]
plt.legend(handels1,labels1)
plt.xlabel('log_minL')
plt.title('$F_{var}$')
plt.yscale('log')
plt.savefig('/media/ashley/KINGSTON/hist_div5.pdf',dpi=200,format='pdf')





file.close()
file1.close()
"""
nameg3=[]
numg3=[]
redg3=[]
namel3=[]
numl3=[]
redl3=[]

for i in range(0,len(name3),1):
    if num_3[i]>3:
        nameg3.append(name3[i])
        numg3.append(num_3[i])
        redg3.append(red_3[i])
    if num_3[i]>=2 and num_3[i]<=3:
        namel3.append(name3[i])
        numl3.append(num_3[i])
        redl3.append(red_3[i])

print('len_nameg3',len(nameg3))
print('len_namel3',len(namel3))
"""
'''
nameg3.extend(namel3)
#print('type',type(allname))
a=list(set(name_2).difference(set(nameg3)))
print(a)
print('test:',namel3[0])
b=name2.index('CXO J121749.7+301151')
print('index_b',b)
'''
['CXO J121749.7+301151',
 'CXO J121345.5+364412',
 'CXO J141756.6+254326',
 'CXO J225715.9+074251', 
 'CXO J093037.5+495025', 
 'CXO J124203.5+063829',
 'CXO J105329.4+572104',
 'CXO J103223.5+532319',
 'CXO J121049.6+392821', 
 'CXO J122911.6+020313',
 'CXO J141953.7+541920', 
 'CXO J165258.8+022403',
 'CXO J103140.2+505436', 
 'CXO J141756.8+254430', 
 'CXO J122915.4+020529',
 'CXO J124155.0+063327',
 'CXO J121032.5+392420']
#lost data
"""
#========num<2,3>
indexl3=[]
varl3=[]
varl3_real=[]
#print('single:',name3.index(namel3[0]))
for i in range(0,len(namel3),1):
    fluxtmp=[]
    index_l3=name3.index(namel3[i])
    #index_tmp=[n for n,v in enumerate(name3) if v==namel3[i]]
    #index_l3=np.min(index_tmp)
    #print('index_l3:',index_l3)
    indexl3.append(index_l3)
    a_num=numl3[i]
    for j in range(index_l3,a_num+index_l3,1):
        if flux3[j]==0.0:
            flux_tmp=flux3_hilim[j]
        else:
            flux_tmp=flux3[j]
        fluxtmp.append(flux_tmp)
    #print('fluxtmp:',fluxtmp)
    flux_tmpmin=np.min(fluxtmp)
    flux_tmpmax=np.max(fluxtmp)
    #print
    vary=(flux_tmpmax)/(flux_tmpmin)
    vary_real=np.var(fluxtmp)/(np.mean(flux_tmp)**2)
    #print('vary_real',vary_real)

    #print('vary',vary)
    varl3.append(vary)
    varl3_real.append(vary_real)
#index=[i for i,v in enumerate(varl3) if v=='nan']
#print('index_nan:',index) 
#print('len_varl3',len(varl3))
#print('var:',var)
ymax1=np.max(varl3_real)*1.1
Dll3=[i*3E8*308568E16/67.73999999999998 for i in redl3]
Dl3=[i*3.24077929e-25 for i in Dll3]
plt.figure()
#plt.ylim(0,ymax1)
grey=plt.scatter(Dl3,varl3_real,s=4,facecolors='none',edgecolor='c',label='3>=acis_num>=2')
plt.legend(loc=0)          
plt.xscale('log')
plt.xlabel('Distance'+'  '+'Mpc')
plt.ylabel('var_real')

plt.title('var_real23')
plt.savefig('/home/ashley/Link_xray/result/var_real23',format='png',dpi=200)


#===========num>3:
indexg3=[]
varg3=[]
varg3_real=[]
#index=name3.index(nameg3[828])
#print(index)
for i in range(0,len(nameg3),1):
    fluxtmp=[]
    index_g3=name3.index(nameg3[i])
    #index_tmp=[n for n,v in enumerate(name3) if v==namel3[i]]
    #index_l3=np.min(index_tmp)
    #print('index_l3:',index_l3)
    indexg3.append(index_g3)
    #print('i',i)
    a_num=numg3[i]
    #print(index_g3)
    #print(a_num)
    for j in range(index_g3,a_num+index_g3,1):
        
        if flux3[j]==0.0:
            flux_tmp=flux3_hilim[j]
        else:
            flux_tmp=flux3[j]
        fluxtmp.append(flux_tmp)
    #print('fluxtmp:',fluxtmp)
    flux_tmpmin=np.min(fluxtmp)
    flux_tmpmax=np.max(fluxtmp)
    #print
    vary=(flux_tmpmax)/(flux_tmpmin)
    vary_real=np.var(fluxtmp)/(np.mean(flux_tmp)**2)
    #print('vary',vary)
    varg3.append(vary)
    varg3_real.append(vary_real)
#print('len_varg3',len(varg3))
#print('var:',var)
#index1=[i for i,v in enumerate(varg3) if v=='nan']
#print('index_nan:',index1) 
ymax=np.max(varg3_real)*1.1
Dgg3=[i*3E8*308568E16/67.73999999999998 for i in redg3]
Dg3=[i*3.24077929e-25 for i in Dgg3]
plt.figure()
#plt.ylim(0,ymax)
glow=plt.scatter(Dg3,varg3_real,s=4,facecolors='none',edgecolor='m',label='acis_num>3')
plt.legend(loc=0)          
plt.xscale('log')
plt.xlabel('Distance'+'  '+'Mpc')
plt.ylabel('var_real')

plt.title('var3_real')
plt.savefig('/home/ashley/Link_xray/result/var_real',format='png',dpi=200)

file3.close()

 
"""
