import numpy as np
import pylab as pl
import pandas as pd
import scipy.integrate as si

f1,ax1 = pl.subplots(1,3,figsize=[12,4])
f2,ax2 = pl.subplots(1,3,figsize=[12,4])


def model(u,t,Sv,phigv,phibv,phigb,deltv,mub,deltb,epsgb,epsgv,epsb,deltg):
    V,B,G = u[0],u[1],u[2]
    dVdt = Sv -phigv*G*V -phibv*B*V - deltv*V
    dBdt = mub*B + epsb*phibv*B*V - phigb*B*G - deltb*B**2
    dGdt = epsgv*phigv*G*V + epsgb*phigb*G*B - deltg*G**2
    return np.r_[[dVdt,dBdt,dGdt]]

# params
Sv = 5e+7
phigv = 1/1e+3
phibv = 1/1e+6
phigb = 1/1e+3
deltv = 0.1
mub = 0.8
deltb = 1e-6
epsgb = 0.1
epsgv = 1e-6
epsb = 1e-4
deltg = 0.1
params = (Sv,phigv,phibv,phigb,deltv,mub,deltb,epsgb,epsgv,epsb,deltg)

# time array
ndays,delt = 10,0.01
times = np.linspace(0,ndays,int(ndays/delt))

# spin-up
inits = np.r_[[1e+3,1e+3,1e+3]]
vs = si.odeint(model,inits,times,args=params).T
S3dicta = {'time':times,'virusa':vs[0],'bacteriaa':vs[1],'grazera':vs[2]}
ax1[0].plot(times,vs[0],label='virus')
ax1[0].plot(times,vs[1],label='bacteria')
ax1[0].plot(times,vs[2],label='grazer')
ax1[0].semilogy()

ndays,delt = 8,0.01
times = np.linspace(0,ndays,int(ndays/delt))

inits = np.r_[[1e+5,1e+5,1e+5]]
vs = si.odeint(model,inits,times,args=params).T
S3dictb = {'time':times,'virusb':vs[0],'bacteriab':vs[1],'grazerb':vs[2]}
ax1[1].plot(times,vs[0],label='virus')
ax1[1].plot(times,vs[1],label='bacteria')
ax1[1].plot(times,vs[2],label='grazer')
ax1[1].semilogy()

inits = np.r_[[1e+7,1e+7,1e+7]]
vs = si.odeint(model,inits,times,args=params).T
S3dictc = {'time':times,'virusc':vs[0],'bacteriac':vs[1],'grazerc':vs[2]}
ax1[2].plot(times,vs[0],label='virus')
ax1[2].plot(times,vs[1],label='bacteria')
ax1[2].plot(times,vs[2],label='grazer')
ax1[2].semilogy()

# save data
pd.DataFrame(S3dicta).to_csv('../data/FigS3panela.csv')
pd.DataFrame(S3dictb).to_csv('../data/FigS3panelb.csv')
pd.DataFrame(S3dictc).to_csv('../data/FigS3panelc.csv')

# labels
lab1 = '<100 micron'
lab2 = '<1.2 micron'
lab3 = '<0.1 micron'

# turn off virus supply
params = (0.0,phigv,phibv,phigb,deltv,mub,deltb,epsgb,epsgv,epsb,deltg)

# <100mum
inits = vs.T[-1]
m1 = si.odeint(model,inits,times,args=params).T
F5dictsmall = {'time':times,'virussmall':m1[0],'bacteriasmall':m1[1],'grazersmall':m1[2]}
ax2[0].plot(times,m1[0]/1e+7,label=lab1,c=str(0.8))
ax2[1].plot(times,m1[1]/1e+5,label=lab1,c=str(0.8))
ax2[2].plot(times,m1[2],label=lab1,c=str(0.8))

# <1.2mum
inits = vs.T[-1]
inits[2] = 0.0
m1 = si.odeint(model,inits,times,args=params).T
F5dictmed = {'virusmed':m1[0],'bacteriamed':m1[1],'grazermed':m1[2]}
ax2[0].plot(times,m1[0]/1e+7,label=lab2,c=str(0.5))
ax2[1].plot(times,m1[1]/1e+5,label=lab2,c=str(0.5))

# <0.1mum
inits = vs.T[-1]
inits[2] = 0.0
inits[1] = inits[1] / 2
m1 = si.odeint(model,inits,times,args=params).T
F5dictlarge = {'viruslarge':m1[0],'bacterialarge':m1[1],'grazerlarge':m1[2]}
ax2[0].plot(times,m1[0]/1e+7,label=lab3,c=str(0.1))
ax2[1].plot(times,m1[1]/1e+5,label=lab3,c=str(0.1))

# save model data
pd.DataFrame({**F5dictsmall,**F5dictmed,**F5dictlarge}).to_csv('../data/figure5all.csv')

for axes in (ax1,ax2):
    for (ax,let) in zip(axes,'abc'):
        ax.text(0.1,0.9,let,transform=ax.transAxes)
        ax.set_xlabel('Day')

ax2[0].set_ylabel(r'Virus (x10$^{7}$  ml$^{-1}$)')
ax2[1].set_ylabel(r'Bacteria (x10$^{5}$  ml$^{-1})$')
ax2[2].set_ylabel(r'Grazer ml$^{-1}$')
ax2[1].set_ylim(ymin=0.0)

ax2[2].set_ylim([0,800])

ax1[0].set_ylabel(r'Abundance (ml$^{-1}$)')

f2.subplots_adjust(wspace=0.3)

l1 = ax1[0].legend()
l2 = ax2[0].legend()
l1.draw_frame(False)
l2.draw_frame(False)

f1.savefig('../figures/spinup',dpi=300,bbox_inches='tight')
f2.savefig('../figures/filtration',dpi=300,bbox_inches='tight')

pl.show()

pl.close(f1)
pl.close(f2)

