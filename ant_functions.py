# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 15:46:46 2015

@author: ab14785
"""

from obspy.fdsn import Client
from obspy import UTCDateTime
from obspy import read, Trace, Stream

import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp

#client = Client()
#
#def getWave(network, station, number, channel, UTC, dur):
#    """
#    Downloads miniseed datasets through the obspy.fdsn function.     
#    """
#    t = UTCDateTime(UTC)
#    st = client.get_waveforms(network, station, number, channel, t, t + dur, attach_response=True)
#    print st
#    return st

def preprocess(stream):
    """Carries out simple preprocessing of trace, by first merging the stream, 
    removing instrumetn response, highpass filtering at 0.2 Hz then tapering"""
#    stream.merge()
    stream.remove_response(output="disp")
    stream.filter('highpass',freq=0.02,corners=2,zerophase=True)
    stream.taper(max_percentage=0.02,type='cosine')
    return stream

def FRF(stream,comp):

    freq_all=[]
    spec_all=[]

    t1=stream[0].stats.starttime
    t2=stream[0].stats.endtime
    samp_rate=stream[0].stats.sampling_rate 

    t=t1
    t_end = t2-60*5

    while t<t_end:
        cut=stream[comp].slice(t,t+60*5)

        spec_pow=abs(np.fft.fft(cut))
        freq_abs= abs(np.fft.fftfreq(len(cut), d=1./samp_rate))


        w=50
        step=50
        freq_par=[freq_abs[i] for i in range(0,len(freq_abs),step)]

        spec_par=[sum(spec_pow[i-(w-step):i+step]) \
                  if i>(w-1) else sum(spec_pow[:i+step]) \
                  for i in range(0,len(spec_pow),step)]
#         print len(spec_pow)

    #     spec_all.append(spec)

#         plt.plot(freq_abs,spec_pow,color = '0.75')  
#         plt.show()
#         plt.plot(freq_par,spec_par,color = '0.75')
#         plt.show
        
        freq_all.append(freq_par)
        spec_all.append(spec_par)
#         print t
        t=t+60*5
#         print t


    spec_mean=[]
    for i in range(0,len(spec_all[0])):
        value=0
        for j in range(0,len(spec_all)):
            value=value+spec_all[j][i]
        average=value/len(spec_all)
        spec_mean.append(average) 

    print len(spec_mean)
    plt.plot(freq_all,spec_all,color='0.75')
    plt.plot(freq_par,spec_mean,color = 'k')
    plt.show
    
    return freq_par,spec_mean

def FRF_all(st_pre):
    freq_all = []
    spec_all = []

    freq,spec = FRF(st_pre,0)
    freq_all.append(freq)
    spec_all.append(spec)
    plt.show()

    freq,spec =FRF(st_pre,1)
    freq_all.append(freq)
    spec_all.append(spec)
    plt.show()

    freq,spec =FRF(st_pre,2)
    freq_all.append(freq)
    spec_all.append(spec)
    plt.show()

    return freq_all,spec_all

def ratio(spec1,spec2):
    spec_ratio=[]
    print len(spec1)
    for i in range(0,len(spec1)):
        ratio = spec1[i]/spec2[i]
        spec_ratio.append(ratio)
    
    return spec_ratio

def Amplification_plot(freq,spec):
    reg= sum(spec[2][30:40])/10

    spec_reg_E=spec[0]/reg
    spec_reg_N=spec[1]/reg
    spec_reg_Z=spec[2]/reg
    spec_EZ=ratio(spec[0],spec[2])
    spec_NZ=ratio(spec[1],spec[2])

    spec_trace=spec_reg_E

    plt.figure(figsize=[16,6])
    plt.plot(freq[0][0:len(spec_trace)/2],spec_reg_E[0:len(spec_trace)/2],color = '0.5',linestyle = '-',label='East')
    plt.plot(freq[1][0:len(spec_trace)/2],spec_reg_N[0:len(spec_trace)/2],color = 'k',linestyle = '--',label='North')
    plt.plot(freq[2][0:len(spec_trace)/2],spec_reg_Z[0:len(spec_trace)/2],color = 'k',linestyle = '-',label='Vertical')
    plt.plot(freq[2][0:len(spec_trace)/2],spec_EZ[0:len(spec_trace)/2],color = 'k',linestyle = ':',label='EZ Ratio')
    plt.plot(freq[2][0:len(spec_trace)/2],spec_NZ[0:len(spec_trace)/2],color = 'k',linestyle = '-.',label='NZ Ratio')
    plt.axhline(y=1,color='0.3',linestyle = ':')
    plt.legend()
    plt.xlabel('Frequency (Hz)',size=14)
    plt.ylabel('Amplification Factor',size=14)

    plt.ylim(0,15)
    plt.show()
    
    return freq[0][0:len(spec_trace)/2], spec_reg_E[0:len(spec_trace)/2]

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), sp.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m, m-h, m+h
    
def latlong_distn(lat1, long1, lat2, long2):
 
    # Convert latitude and longitude to 
    # spherical coordinates in radians.
    degrees_to_radians = np.pi/180.0
         
    # phi = 90 - latitude
    phi1 = (90.0 - lat1)*degrees_to_radians
    phi2 = (90.0 - lat2)*degrees_to_radians
         
    # theta = longitude
    theta1 = long1*degrees_to_radians
    theta2 = long2*degrees_to_radians
         
    # Compute spherical distance from spherical coordinates.
         
    # For two locations in spherical coordinates 
    # (1, theta, phi) and (1, theta', phi')
    # cosine( arc length ) = 
    #    sin phi sin phi' cos(theta-theta') + cos phi cos phi'
    # distance = rho * arc length
     
    cos = (np.sin(phi1)*np.sin(phi2)*np.cos(theta1 - theta2) + 
           np.cos(phi1)*np.cos(phi2))
    arc = np.arccos( cos )
 
    # Remember to multiply arc by the radius of the earth 
    # in your favorite set of units to get length.
    return (arc*6371)

def ml(a,r):
    """equation to calculte local magnitude using the UK equation"""
    mag=(np.log10(a))+(0.95*np.log10(r))+(0.00183*r)-1.76
    
    return mag

def ml_corr(a,r,stat_corr):
    """equation to calculte local magnitude using the UK equation"""
    mag=(np.log10(a))+(0.95*np.log10(r))+(0.00183*r)+stat_corr-1.76
    
    return mag

def ml_amp(mag,distance):
    amplitude = -(mag+1.76)+(0.95*np.log10(distance)+(0.00183*distance))
    return amplitude

def fft(trace,samp_rate):
    """Enhanced fft using numpy, function takes a trace and returns frequency, amplitude and phase arrays"""
    
    spec=np.fft.fft(trace)
    freq=np.fft.fftfreq(len(trace),d=1./samp_rate)
    
    spec_real=spec.real
    spec_imag=spec.imag
    
    amplitude =  np.sqrt(spec_real**2 + spec_imag**2)
    phase = np.arctan2(spec_imag,spec_real)
    
    return freq, amplitude, phase

def invfft(amplitude,phase):
    """Inverse fft with take the amplitude and phase and returns a time trace"""
    re = amplitude * np.cos(phase)
    im = amplitude * np.sin(phase)
    F = re + 1j*im;
    inv_F=np.fft.ifft(F)
    
    return inv_F
    

def parseval(frequency,amplitude,phase,bin_hz):
    """Bins data at user defined intervals"""
    
    nyquist = int(np.max(frequency))
    bin_size=int((len(amplitude)/2*bin_hz)/nyquist)

    freq_par=[]
    spec_par=[]
    phase_par=[]
    for i in range(0,len(amplitude)/2,bin_size):
        if i+bin_size<len(amplitude):
            freq_par.append(abs(frequency[i+bin_size/2]))
            spec_par.append(sum(amplitude[i:i+bin_size]))
            phase_par.append(sum(phase[i:i+bin_size]))
        else:
            pass
        
    return freq_par,spec_par,phase_par

def cft(st):
    cft_en=[]
    # Creates the characteristic function which is the absolute value of all three components
    for i in range(0,len(st[0])):
        en_temp=np.sqrt(st[0][i]**2+st[1][i]**2+st[2][i]**2)
        cft_en.append(en_temp)
    return cft_en
        
def noise_filter(st):
    cft_en=cft(st)
    # Determine the sample rate
    dt=st[0].stats.sampling_rate
    # length of characteristic function
    end=len(cft_en)
    # sample rate x 60sec x number of minutes. The characteristic function is slice into these lengths.
    length=dt*60*4

    # creates a dataframe of the standard deviations and start times.
    data={'start':[],'std':[]}
    for j in range(0,int(end),int(length)):
        data['start'].append(j)
        cft_temp=cft_en[j:j+int(length)]
        data['std'].append(np.std(cft_temp))
    else:
        # Required as the last section of the trace is unlikely to be the correct length
        pass

    frame=pd.DataFrame(data)
    plt.plot(frame['std'])
    sort_frame=frame.sort_index(by='std')
    sort_frame2=pd.DataFrame(sort_frame, columns=['start','std','rank'])
    sort_frame2['rank']=range(0,len(sort_frame2))
    sort_frame2=sort_frame2.set_index(['rank']) 
    plt.plot(sort_frame2['std'])
    plt.show()

    diff=np.gradient(sort_frame2['std'])
    plt.plot(diff)
    plt.show()
    return sort_frame2
    
def noise_chop(sort_frame2,n1,n2):
    
    plt.plot(sort_frame2['std'])
    plt.plot(sort_frame2['std'][n1:n2])
    plt.show()

    filtered = pd.DataFrame(sort_frame2[n1:n2], columns=['start','std','rank'])
    filtered_sort=filtered.sort_index(by='start')
    filtered_sort['rank']=range(0,len(filtered))
    filtered_sort=filtered_sort.set_index(['rank']) 

    plt.plot(filtered_sort['start'],filtered_sort['std'])
    
    return filtered_sort

def noise_fft(tr,filtered_sort):
    
#    cft_en=cft(st)
    
    dt=tr.stats.sampling_rate
    # length of characteristic function
    end=len(cft_en)
    # sample rate x 60sec x number of minutes. The characteristic function is slice into these lengths.
    length=dt*60*4
    
    freq_all = []
    spec_all = []
    flt_trace=[]

    plt.plot(cft_en)
    plt.show()

    for i in range(0,len(filtered_sort)):
#     for i in range(0,5):
        start=filtered_sort['start'][i]

        trace_temp=tr.data[start:int(start+length)]
        plt.plot(trace_temp)
    #     plt.ylim(-0.0000003,0.0000003)
        plt.show()
        freq, amplitude, phase=fft(trace_temp,dt)
        freq_par,spec_par,phase_par=parseval(freq, amplitude, phase,0.1)

        plt.plot(freq_par,spec_par)
        plt.show()
        freq_all.append(freq_par)
        spec_all.append(spec_par)
        flt_trace.append(trace_temp)

    flt_trace=np.hstack(flt_trace)
    plt.plot(cft_en)
    plt.plot(flt_trace)
    plt.show()

    spec_mean=[]
    five_all=[]
    ninetyfive_all=[]

    for i in range(0,len(spec_all[0])):
        value=[]
        for j in range(0,len(spec_all)):
            value.append(spec_all[j][i])
    #     print len(value)
        average,five,ninetyfive = mean_confidence_interval(value,0.95)

        spec_mean.append(average)
        five_all.append(five)
        ninetyfive_all.append(ninetyfive)
    print len(spec_mean)

    plt.plot(freq_all[0],spec_mean,'k')
    plt.plot(freq_all[0],five_all,'k:')
    plt.plot(freq_all[0],ninetyfive_all,'k:')
    plt.ylim(0,0.0001)
    plt.show()
    
    return freq_all,spec_mean

def sac_write(i):  
    station=st[i].stats.station
    time=st[i].stats.starttime
    channel=st[i].stats.channel
    resp_file='RESP/RESP.GB.%s.00.%s'%(station,channel)
    
    if not os.path.exists(resp_file):
        pass
    
    else:

        resp_file='RESP/RESP.GB.%s.00.%s'%(station,channel)
        date=st[i].stats.starttime
        pre_filt = (0.2,0.5,49,50)
        seedresp = {'filename': resp_file, 'date':date, 'units': 'DIS'}    
        st[i].simulate(paz_remove=None, pre_filt=pre_filt, seedresp=seedresp)
        
        directory=str(time.date)
        if not os.path.exists('SAC/%s'%directory):
            os.makedirs('SAC/%s'%directory)
        name='SAC/%s/GB.%s_%s_%02d-%02d_%s.sac'%(time.date,station,time.date,time.hour,time.minute,channel)
        print(name)
        st[i].write(name, format="SAC")
