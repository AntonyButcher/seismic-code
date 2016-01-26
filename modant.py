# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 12:56:49 2015

@author: Antony Butcher

This contains a series of function for displaying modelled synthetics. It has
a series of function to cut out the p and s arrivals, and display them 
in order to examine changes in the time and frequency domain, along with 
looking at showing the particle motion. 
"""

from obspy.fdsn import Client
from obspy import UTCDateTime
from obspy import read, Trace, Stream

import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.mlab as ml
import scipy
from scipy.interpolate import griddata

from math import sqrt

import os as os

def fft(trace):
    length = len(trace)
    samp_rate=trace.stats.sampling_rate
        
    data=trace.data
    spec_amp=np.fft.fft(data,n=length)
    spec_abs=abs(spec_amp)
        
    freq = np.fft.fftfreq(length, d=1./samp_rate)
    freq_abs=abs(freq)

    return freq_abs,spec_abs
    
def pscut(vp,vs,depth,spacing,stz_sim=stz_sim,stz_ns=stz_ns,stx_sim=stx_sim,stx_ns=stx_ns):
    
    start= stz_sim[0].stats.starttime
    length = len(stz_ns)
    
    stz_ns_cutp=stz_ns.slice(start+0.3,start+0.5)
    stx_ns_cutp=stx_ns.slice(start+0.3,start+0.5)
    
    stz_sim_cutp=stz_sim.slice(start+0.3,start+0.5)
    stx_sim_cutp=stx_sim.slice(start+0.3,start+0.5)
    
    stz_ns_cuts=stz_ns.slice(start+0.3,start+0.5)
    stx_ns_cuts=stx_ns.slice(start+0.3,start+0.5)
    
    stz_sim_cuts=stz_sim.slice(start+0.3,start+0.5)
    stx_sim_cuts=stx_sim.slice(start+0.3,start+0.5)
    
    
    
    for i in range(0,length):
        
        # Calculates the onset time of the p and s wave for an 
        # array with a 100m station spacing
        st_pwave=(sqrt(depth**2+(25*i)**2))/vp
        st_swave=(sqrt(depth**2+(25*i)**2))/vs
        
        # Cuts the trace for the p-wave for the near surface model
        stz_ns_cutp[i]=stz_ns[i].slice(start+st_pwave,start+st_pwave+0.1)
        stx_ns_cutp[i]=stx_ns[i].slice(start+st_pwave,start+st_pwave+0.1)
    
        # Cuts the trace for the p-wave for the simple model
        stz_sim_cutp[i]=stz_sim[i].slice(start+st_pwave,start+st_pwave+0.1)
        stx_sim_cutp[i]=stx_sim[i].slice(start+st_pwave,start+st_pwave+0.1)
        
        # Cuts the trace for the s-wave for the near surface model
        stz_ns_cuts[i]=stz_ns[i].slice(start+st_swave,start+st_swave+0.1)
        stx_ns_cuts[i]=stx_ns[i].slice(start+st_swave,start+st_swave+0.1)
    
        # Cuts the trace for the s-wave for the simple model
        stz_sim_cuts[i]=stz_sim[i].slice(start+st_swave,start+st_swave+0.1)
        stx_sim_cuts[i]=stx_sim[i].slice(start+st_swave,start+st_swave+0.1)
        
    return stz_ns_cutp,stx_ns_cutp,stz_sim_cutp,stx_sim_cutp,stz_ns_cuts,stx_ns_cuts,stz_sim_cuts,stx_sim_cuts
    
    
def ratio(spec1,spec2):
    spec_tmp1=np.array(spec1)
    spec_tmp2=np.array(spec2)
    
    ratio=spec_tmp1/spec_tmp2
    
    return ratio

def stnplt(stream1,stream2,spacing=1):
    """Overlays traces from modelled data. Assumes you have a simple model that 
    you wish to plot with another more complicate model. Requires a station 
    spacing for labelling, otherwise defaults to 1"""
    f, ax = plt.subplots(len(stream1),sharex=True, sharey=True,figsize=(12,50))
    
    plt.subplots_adjust(hspace=0.001)
    
    plot=0
    length = len(stream1)
    
    for j in range(0,length):
        ax[plot].plot(stream1[j].times(),stream1[j].data,'g', )
        ax[plot].plot(stream2[j].times(),stream2[j].data,'b')
        ax[plot].text(0.6, 0.9, 'Station at %dm' %(plot*spacing),transform=ax[plot].transAxes, fontsize=12)
                
        plot=plot+1
        
    
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    plt.tight_layout()
    
    if not os.path.exists('images'):
        os.makedirs('images')
    
    plt.savefig('images/traces.pdf')
    
    
def overplt(vp,vs,depth,spacing,stz_sim=stz_sim,stz_ns=stz_ns,stx_sim=stx_sim,stx_ns=stx_ns):
    
    length = len(stz_ns)    
    size= len(stz_ns[0])
    end= size*stz_ns[0].stats.delta    
    
    freq_simple_z=np.array(stz_ns)
    spec_simple_z=np.array(stz_ns)
    freq_ns_z=np.array(stz_ns)
    spec_ns_z=np.array(stz_ns)
    
    freq_simple_x=np.array(stz_ns)
    spec_simple_x=np.array(stz_ns)
    freq_ns_x=np.array(stz_ns)
    spec_ns_x=np.array(stz_ns)
    
    ratio_x=np.array(stz_ns)
    ratio_z=np.array(stz_ns)
    
    ns_zxratio=np.array(stz_ns)
    
    
    for i in range(0,length):
        freq_simple_z[i],spec_simple_z[i]=fft(stz_sim[i])
        freq_ns_z[i],spec_ns_z[i]=fft(stz_ns[i])
    
        ratio_z[i] = ratio(spec_ns_z[i],spec_simple_z[i])
    
        freq_simple_x[i],spec_simple_x[i]=fft(stx_sim[i])
        freq_ns_x[i],spec_ns_x[i]=fft(stx_ns[i])
    
        ratio_x[i] = ratio(spec_ns_x[i],spec_simple_x[i])
    
        ns_zxratio[i] = ratio(spec_ns_x[i],spec_ns_z[i])
        
# this doesn't work as it's not incrementing        
    stz_ns_cutp,stx_ns_cutp,stz_sim_cutp,stx_sim_cutp,stz_ns_cuts,stx_ns_cuts,stz_sim_cuts,stx_sim_cuts=pscut(vp=vp,vs=vs,depth=depth,spacing=spacing)


    
    for j in range(0,length):
    
        fig = plt.figure(figsize=[12,9])
        fig.subplots_adjust(hspace=0.3, wspace=0)
                
        ax2 = fig.add_subplot(3, 2, 2)
    #    ax2.set_ylim([-2e-7,2e-7])
        ax2.set_xlim([0.1,end])
        ax2.xaxis.tick_top()
        ax2.grid()
        p = plt.axvspan(stz_sim_cutp[j].stats.starttime,stz_sim_cutp[j].stats.starttime+0.1, facecolor='0.5', alpha=0.2)
        p = plt.axvspan(stz_sim_cuts[j].stats.starttime,stz_sim_cuts[j].stats.starttime+0.1, facecolor='0.5', alpha=0.2)
        ax2.text(0.65, 0.02, 'Horizontal Component',transform=ax2.transAxes, fontsize=10)
        ax2.plot(stx_sim[j].times(),stx_sim[j].data,label='Simple')
        ax2.plot(stx_ns[j].times(),stx_ns[j].data,label='NS Layer')
        
        plt.legend(fontsize=10)
        
        ax1 = fig.add_subplot(3, 2, 1,sharey=ax2)
    #    ax1.set_ylim([-2e-7,2e-7])
        ax1.set_xlim([0.1,end]) 
        ax1.xaxis.tick_top()
        ax1.grid()
        p = plt.axvspan(stz_sim_cutp[j].stats.starttime,stz_sim_cutp[j].stats.starttime+0.1, facecolor='0.5', alpha=0.2)
        p = plt.axvspan(stz_sim_cuts[j].stats.starttime,stz_sim_cuts[j].stats.starttime+0.1, facecolor='0.5', alpha=0.2)
        ax1.text(0.7, 0.9, 'Station at %dm' %(j*25),transform=ax1.transAxes, fontsize=12)
        ax1.text(0.7, 0.02, 'Vertical Component',transform=ax1.transAxes, fontsize=10)
        ax1.plot(stz_sim[j].times(),stz_sim[j].data,label='Simple')
        ax1.plot(stz_ns[j].times(),stz_ns[j].data,label='NS Layer')
        
        
        yticklabels = ax2.get_yticklabels()
        plt.setp(yticklabels, visible=False)
        
    #    fig = plt.figure(figsize=[13,3])
    #    fig.subplots_adjust(hspace=None, wspace=0.3)
        
        ax4 = fig.add_subplot(3, 4, 6)
    #    ax4.set_ylim([-4e-8,4e-8])
    #    ax4.set_xlim([-4e-8,4e-8])
        ax4.grid()
        ax4.set_title('NS-P Wave',fontsize=10)
    #    ax4.set_ylabel('z-component')
        ax4.plot(stx_ns_cutp[j].data,stz_ns_cutp[j].data, 'g',label='NS Layer')
        
        ax3 = fig.add_subplot(3, 4, 5,sharey=ax4,sharex=ax4)
    #    ax3.set_ylim([-4e-8,4e-8])
    #    ax3.set_xlim([-4e-8,4e-8])
        ax3.grid()
        ax3.set_title('Sim-P Wave',fontsize=10)
        ax3.set_ylabel('z-component')
        ax3.plot(stx_sim_cutp[j].data,stz_sim_cutp[j].data, 'b',label='NS Layer')
        
        ax6 = fig.add_subplot(3, 4, 8)
    #    ax6.set_ylim([-3e-7,3e-7])
    #    ax6.set_xlim([-3e-7,3e-7])
        ax6.grid()
        ax6.set_title('NS-S Wave',fontsize=10)
    #    ax6.set_ylabel('z-component')
        ax6.plot(stx_ns_cuts[j].data,stz_ns_cuts[j].data, 'g',label='NS Layer')
        
        ax5 = fig.add_subplot(3, 4, 7,sharey=ax6,sharex=ax6)
    #    ax5.set_ylim([-3e-7,3e-7])
    #    ax5.set_xlim([-3e-7,3e-7])
        ax5.grid()
        ax5.set_title('Sim-S Wave',fontsize=10)
        ax5.set_xlabel('x-component')
        ax5.plot(stx_sim_cuts[j].data,stz_sim_cuts[j].data, 'b',label='NS Layer')
        
        
        yticklabels = ax6.get_yticklabels()+ax4.get_yticklabels()
        plt.setp(yticklabels, visible=False)
        
        ax7 = fig.add_subplot(3, 4, 9)
        ax7.plot(freq_simple_z[j], spec_simple_z[j],label='Simple')
        ax7.plot(freq_ns_z[j], spec_ns_z[j],label='NS Layer')
    #    ax7.title('z-component spectrum')
        ax7.set_ylabel('Amplitude')
        ax7.set_xlabel('Frequency(Hz)')
        ax7.set_xlim([0,50])
        ax7.set_ylim([0,8e-5])
        ax7.legend(fontsize=8)
    
        ax8 = fig.add_subplot(3, 4, 10)
        ax8.plot(freq_simple_x[j], spec_simple_x[j],label='Simple')
        ax8.plot(freq_ns_x[j], spec_ns_x[j],label='NS Layer')
        ax8.set_xlim([0,50])
        ax8.set_ylim([0,8e-5])
        
        ax9 = fig.add_subplot(3, 4, 11)
        ax9.plot(freq_simple_x[j], ratio_z[j],label='Simple')
        ax9.plot(freq_ns_x[j], ratio_x[j],label='NS Layer')
        ax9.set_xlim([0,50])
        ax9.set_ylim([0,10])
        
        ax10 = fig.add_subplot(3, 4, 12)
        ax10.plot(freq_simple_x[j], ns_zxratio[j],label='Simple')
        ax10.set_xlim([0,50])
        ax10.set_ylim([0,30])
    
        yticklabels = ax8.get_yticklabels()+ax10.get_yticklabels()
        plt.setp(yticklabels, visible=False)
        
        if not os.path.exists('images'):
            os.makedirs('images')
        plt.savefig('images/m14_4m_stn%dm.pdf' %(j*25))
        plt.show
    
    
    
if __name__=="__main__":
    stz_ns = read('data/*ns*.z',format='SAC')
    stz_sim = read('data/*sim*.z',format='SAC')
    stx_ns = read('data/*ns*.x',format='SAC')
    stx_sim = read('data/*sim*.x',format='SAC')
    
    stnplt(stz_sim,stz_ns,25)
    
    overplt(vp=3000,vs=1500,depth=1000,spacing=25)
    #stz_ns_cutp,stx_ns_cutp,stz_sim_cutp,stx_sim_cutp,stz_ns_cuts,stx_ns_cuts,stz_sim_cuts,stx_sim_cuts = pscut(vp=vp,vs=vs,depth=depth,spacing=spacing)
    
    
