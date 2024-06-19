import numpy as np
import pylab as pl

f,ax = pl.subplots()
times = np.linspace(0,8,1000)

m100,m02 = 0.5,0.2
V0 = 30
Vs100 = V0*np.exp(-m100*times)
Vs02 = V0*np.exp(-m02*times)

#ax.plot(times,Vs100,c=str(0.5),label=r'$<$100 $\mu$m')
#ax.plot(times,Vs02,c=str(0.0),label=r'$<$0.2 $\mu$m')

ax.text(0.3,0.8,r'simple decay model:'+'\n'+ r'$\frac{dV}{dt}=-m_vV$',transform=ax.transAxes)

ax.plot(times,Vs100,c=str(0.5),label=r'$m_v$='+str(m100))
ax.plot(times,Vs02,c=str(0.0),label=r'$m_v$='+str(m02))

l = ax.legend(loc='best')
l.draw_frame(False)

ax.set_ylabel(r'EhVs ml$^{-1}$ x10$^{6}$')
ax.set_xlabel(r'Day')

pl.show()
f.savefig('../figures/exp_decay_model_exp2')

pl.close(f)
