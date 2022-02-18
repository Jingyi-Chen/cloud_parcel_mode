import matplotlib.pyplot as plt
import numpy as np
def set_font_func(ftsz):
    # set the fontsize for the plots
    plt.rcParams['xtick.labelsize']=ftsz
    plt.rcParams['ytick.labelsize']=ftsz
    plt.rcParams['axes.labelsize'] = ftsz
    plt.rcParams['legend.fontsize'] = ftsz
    plt.rcParams['axes.titlesize'] = ftsz
    
cpmv = np.loadtxt('cpmv.txt')
step = cpmv[:,0]
height = cpmv[:,2]

var = {}
var['temp'] = cpmv[:,3]
var['press'] = cpmv[:,4]*1.0e-3 #hPa
var['sat'] = cpmv[:,6]
var['vmr'] = cpmv[:,11]*1.0e3
var['lmr'] = cpmv[:,12]*1.0e3
var['lwc'] = cpmv[:,14]*1.0e6
var['cdconc'] = cpmv[:,15]
var['rmean'] = cpmv[:,16]*1.0e4

units = ['K','1','g $kg^{-1}$','g $kg^{-1}$','kg $m^{-3}$','$cm^{-3}$','$\mu$m']
ii=0
for var_name in ['temp','sat','vmr','lmr','lwc','cdconc','rmean']:
    
    set_font_func(20)
    (fig,axes) = plt.subplots(1,1,figsize=(8,6),sharey=True)
    axes.plot(var[var_name],height*0.01)
    axes.set_xlabel(var_name + f' (unit: {units[ii]})')
    axes.set_ylabel('Height (m)')
    ii += 1

    fig.savefig(f'./{var_name}',dpi=fig.dpi,bbox_inches='tight', pad_inches=0.5)
    plt.close()
    
    
    
