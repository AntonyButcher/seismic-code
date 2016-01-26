# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 09:53:24 2015

@author: Antony Butcher

These functions relate to seismic refraction data. 
"""
from obspy import read, Trace, Stream
import numpy as np 
import matplotlib.pyplot as plt


import glob 
import os as os


def refrac_plot(shot_file):
    """ This function creates a frequency plot for each trace in the shot 
    gather and creates a avi file. Shot_file should be in the form 'raw/2*.dat'"""
    traces = glob.glob(shot_file)
    inc = 1
    
    for f in traces:
        
        filelist = glob.glob("refrac_images/*.png")
        for r in filelist:
            os.remove(r)
            
        x= read(f)
        print f
        len_x=len(x)
        print len_x
    
        t1=x[0].stats.starttime
        t2=x[0].stats.endtime
        samp_rate=x[0].stats.sampling_rate 
    
        for i in range(0,len_x):
            spec=np.fft.fft(x[i])
            freq= np.fft.fftfreq(len(x[i]), d=1./samp_rate)
            plt.figure(figsize=[12,4])
            plt.plot(abs(freq),abs(spec))
            plt.xlim(0,150)
#             plt.ylim(0,8e8)
            plt.title(i*2)
            plt.savefig('refrac_images/tmp_tr%02d.png'%(i),format='png')
            plt.show()
        os.system('mencoder mf://refrac_images/*.png -mf w=800:h=600:fps=2:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o refrac_images/10%02d.avi'%inc)
        
        inc = inc+1
"""
strips out refraction tomography data from plotrefa files
"""

def plotrefa_strip(path):
    files = glob.glob(path)
    print files            
                
    for file in files:
    
    # this first section reads the file in a line at a time, and strips out the cell coordinates and the velocity values
    # The tomo data is separated by 2 white spaces, and there is a '0' between the cell and velocity values.
    
        tomocell=[]
        tomodata=[]
        flag=0
        with open(file,'rt') as f:
            for line in f:
                data = [x.strip() for x in line.split(' ')]
                if data[0]=='':
                    data2 = [x.strip() for x in line.split('  ')]
                    if data2[0]!='0':
                        if flag==0:
                            tomocell.append(data2)
                        else:
                            tomodata.append(data2)
                    elif flag==0:
                        flag=1
                    else:
                        break
        f.close()
    
# This cleans up the arrays and removes unwanted values from the cell and data array
        for i in tomocell:
            i.remove('')
    
        tomocell_np=np.asarray(tomocell)
        tomocell_np=tomocell_np.astype(float)
    
        tomodata_2=[]
    
        for i in range (3,len(tomodata),3):
            tomodata_2.append(tomodata[i])
    
        for i in tomodata_2:
            i.remove('')
    
        tomodata_np=np.asarray(tomodata_2)
        tomodata_np=tomodata_np.astype(float)
        
# separates out the x and y values, then stacks the coordinates with the velocities.
        x=[]
        y=[]
    
        tomo=[]
    
        for i in range(0,len(tomodata_np)):
            x.append(tomocell_np[i*2])
            y.append(tomocell_np[i*2+1])
    
        for i in range(0,len(x)):
            for f in range(0,len(x[i])):
                a= x[i][f]
                b= y[i][f]
                c= tomodata_np[i][f]
                tomo.append(np.hstack((a, b,c)))
        tomo_np=np.asarray(tomo)  
        
        base=os.path.basename(file)
        np.savetxt("output/"+base+"_out", tomo_np, delimiter=",")
        
def s_stack(path):
    traces = glob.glob("S_Data/LineG/segy/raw/*.SGY")
    
    for t in range(0,len(traces),2):
        x1 = read(traces[t],format='SEGY')
        x2 = read(traces[t+1],format='SEGY')
        
        x_1 = x1.copy()
        x_2 = x2.copy()
        x_s = x1.copy()
        
        for i in range(0,len(x1)):

            d1= x_1[i].data
            d2=x_2[i].data
            ds=x_s[i].data
            for j in range(0,len(d1)):
                ds[j]=d1[j]+(d2[j]*-1)
                
            x_s[i].data=ds
        
        print t
        x_s.write('S_Data/LineG/segy/combined/LineG_%02d_ed.SGY'%(1+t/2),format="segy")
        
if __name__=="__main__":    
    
    
    
