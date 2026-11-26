# Compare Migdal Kadanoff ladder running with SU(2) string tension data (Wilson action).
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
beta=np.array([2.2,2.3,2.4,2.5,2.635]); asig=np.array([0.4690,0.3690,0.2660,0.1871,0.1207])
N=2; b0=11*N/(48*np.pi**2); b1=34*N**2/(3*(16*np.pi**2)**2)
def two_loop(b): g2=4*N/(2*b) if False else 2*N/b; return (b0*g2)**(-b1/(2*b0**2))*np.exp(-1/(2*b0*g2))
one=lambda b: np.exp(-b/(4*N*b0)*1.0)  # exp(-1/(2 b0 g^2)) with g^2=2N/b
mk=lambda b: 2.0**(-b/0.25)            # one halving of a per Delta beta = 1/4
slope_data=np.polyfit(beta,np.log(asig),1)[0]
loc=np.diff(np.log(asig))/np.diff(beta)
print("fit slope",slope_data,"local",loc)
print("MK slope",-np.log(2)/0.25,"one loop",-1/(4*N*b0))
bb=np.linspace(2.2,2.635,50)
for f,n in [(one,"one"),(two_loop,"two")]:
    print(n,"local slope at 2.5",(np.log(f(2.51))-np.log(f(2.49)))/0.02)
fig,ax=plt.subplots(figsize=(5.4,3.6))
ax.semilogy(beta,asig,'o',color='#c0392b',label='lattice data')
for f,c,l in [(mk,'#1f4e79','Migdal Kadanoff ladder'),(one,'#7f8c8d','one loop'),(two_loop,'#e67e22','two loop')]:
    ax.semilogy(bb,f(bb)*asig[3]/f(2.5),color=c,label=l)
ax.set_xlabel(r'$\beta$'); ax.set_ylabel(r'$a\sqrt{\sigma}$'); ax.legend(frameon=False,fontsize=8)
ax.set_title('Running of the SU(2) lattice spacing',fontsize=10)
fig.tight_layout(); fig.savefig('figures/scaling_test.pdf')
