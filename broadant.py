"""
@author: Antony Butcher

This contains a series of functions relating to broadband seismic stations. 

"""

from obspy.fdsn import Client
from obspy import UTCDateTime

import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.mlab as ml
import scipy
from scipy.interpolate import griddata

import math

client = Client()

def getWave(network, station, number, channel, UTC, dur):
    """
    Downloads miniseed datasets through the obspy.fdsn function.     
    """
    t = UTCDateTime(UTC)
    st = client.get_waveforms(network, station, number, channel, \
    t, t + dur, attach_response=True)
    print st
    return st
    
def preprocess(stream):
    """Carries out simple preprocessing of trace, by first merging the stream, 
    removing instrument response, highpass filtering at 0.2 Hz then tapering"""
    stream.merge()
    stream.remove_response(output="vel")
    stream.filter('highpass',freq=0.02,corners=2,zerophase=True)
    stream.taper(max_percentage=0.01,type='cosine')
    return stream
    
def spec_amp(stream):
    """This produces three 1d arrays (time,data, freq) constructed after splitting 
    the noise trace into 20min intervals, with the intention that they will be used to create a 
    spectrogram. Using a while loop the routine cuts the trace, creates the time
    array using the mid point of the cut trace, carries out an fft on the cut 
    trace and creates an array of frequencies. These are then stacked onto 
    previously constructed datasets. Note, at the current time I haven't really 
    got my head around structuring np arrays, hence the fudges for reordering arrays."""
    
    data=np.array((0))
    time=np.array((0))
    freq=np.array((0))
    
    t1=stream[0].stats.starttime
    t2=stream[0].stats.endtime
    samp_rate=stream[0].stats.sampling_rate 
    
    t=t1
    
    while t < t2:
        cut=stream[0].slice(t,t+1200)
        tmid=t+600
        t=t+1200
        length = len(cut)-1
        
        time_temp=np.zeros(shape=(length,1))
        
        for a in range(0,length):
            time_temp[a]=[tmid.timestamp]

        time=np.vstack((time,time_temp))
        
        cut_data=cut.data
        cut_spec=np.fft.fft(cut_data,n=length)
    
        data_temp=np.zeros(shape=(length,1),dtype=np.complex)
        
        for b in range(0,length):
            data_temp[b]=cut_spec[b]
    
        data=np.vstack((data,data_temp))
        
        freq_temp = np.fft.fftfreq(length, d=1./samp_rate)
    
        freq_temp2=np.zeros(shape=(length,1))
        
        for c in range(0,length):
            freq_temp2[c]=freq_temp[c]   
    
    
        freq=np.vstack((freq,freq_temp2))
        
    time=np.delete(time,0,0)
    data=np.delete(data,0,0)
    freq=np.delete(freq,0,0)
    
    time=time-time.min()+600

    
    return time,data,freq
    
def Contour(time,freq,data):
    """this function creates a grid of point and interpolates the real data 
    onto the grid"""
    x=np.array(())
    y=np.array(())
    z=np.array(())
    
    x = np.append(x,time)
    y = np.append(y,abs(freq))
    z = np.append(z,abs(data))
    
    points = x,y
    values = z
    
    ndata = len(data)
    ny, nx = 1, 1200
    xmin, xmax = x.min(), x.max()
    ymin, ymax = y.min(), y.max()
    
    xi, yi = np.mgrid[x.min():x.max():nx, ymin:ymax:ny]
    
    zi = griddata(points, values, (xi, yi), method='nearest') 

    
    return xi,yi,zi
    
def fft(trace):
    """fft carried out on the nominated trace and returns both frequency 
    and spectral array"""
    length = len(trace)
    samp_rate=trace.stats.sampling_rate
        
    data=trace.data
    spec_amp=np.fft.fft(data,n=length)
    spec_abs=abs(spec_amp)
        
    freq = np.fft.fftfreq(length, d=1./samp_rate)
    freq_abs=abs(freq)

    return freq_abs,spec_abs
    
    
def HVSR(trace):
    """Simple HV ratio function which sums, according to Parseval's theorem, to 
    smooth the trace"""
    length = len(trace)
    samp_rate=trace.stats.sampling_rate
        
    data=trace.data
    spec_amp=np.fft.fft(data,n=length)
    spec_pow=abs(spec_amp)*abs(spec_amp)
        
    freq = np.fft.fftfreq(length, d=1./samp_rate)
    freq_abs=abs(freq)


    w=200
    step=200
    freq_par=[freq_abs[i] for i in range(0,len(freq_abs),step)]
    spec_par=[sum(spec_pow[i-(w-step):i+step]) if i>(w-1) else sum(spec_pow[:i+step]) for i in range(0,len(spec_pow),step)]

    return freq_par,spec_par    
    
def ml(a,r):
    """equation to calculte local magnitude using the UK equation"""
    mag=math.log(a)+0.95*math.log(r)+0.00183*r-1.76
    
    return mag
    
def snell(th1,v1,v2):
    th_temp=math.asin(v1/v2.math.sin(th1*math.pi/180))
    return th_temp*180/math.pi
    
def f0(hz,vel):
    """Calculate the Fundamental frequency"""
    h=vel/(4.*hz)
    
    return h
    
def tt_diff(depth,vel1,vel2):
    """ Calculates the travel time difference as a result of velocity differences.
    eg the difference between a p- and s-wave arrival."""
    tt_diff=depth*(1./vel1-1./vel2)
    return tt_diff
    
def peak_pick(stream,level):
    """ picks the events above a certain amplitude level. 'stream' should be /
    three E,N,Z traces and level should be the require amplitude threshold"""
# The trigger list is created from the east component, which is normally the
# first trace in the stream. 
    pick_tr=stream[0]
    trig_list = pick_tr.times()[pick_tr.data>level]
    trig_ed=[]
    
# This rationalises the list so that it just contains the first trigger time. 
# Any trigger time within 5sec of the first is removed from the list, change elif statement
# to increase or decrease this criteria
    for i in range(0,len(trig_list)):
        if i == 0:
            trig_ed.append(trig_list[0])
        elif (trig_list[i]-trig_list[i-1])>5:
            trig_ed.append(trig_list[i])
        else:
            pass
     
# This loop creates a 2.5sec trace of the event, with 0.5sec before trigger and 2sec after.
# The sac file is then saved in a pick folder with the date and time.
    for tr in stream:
        tr_cut=tr.copy()
        start=tr_cut.stats.starttime
        for trig in trig_ed:
            t1=start+trig-0.5
            t2=start+trig+2
            tr_cut=tr.slice(t1,t2)
            tr_cut.plot()
            date=str(t1.date)
            time=str(t1.time)
            channel=tr_cut.stats.channel
            name='picks/'+date+'_'+time+'_'+channel+'.sac'
            tr_cut.write(name, format="SAC")    
    
if __name__=="__main__":    
    
    st = getWave('GB','HTL','*', 'BHZ', '2010-02-15T00:00:00.000',7200)
    st_pre = preprocess(st)
    time,data,freq=spec_amp(st_pre)
    xi,yi,zi=Contour(time,freq,data)
    
    fig = plt.figure()
    plt.contourf(xi, yi, zi,1000) 
    plt.clim(0,0.000005)
    plt.show()