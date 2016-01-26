# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 11:37:32 2016

@author: ab14785
"""
import pandas as pd
from datetime import *
import glob


magtype_map={"L":"ML","b":"mb","B":"mB","s":"Ms","S":"MS","W":"MW","G":"mbLg","C":"Md"}

def lntype(line):
  if len(line) == 80:
    type=line[79]
  else:
    type=''
  return type

def read_file(filename):
    try:
        rea_fh=open(filename,"r")
        raw=rea_fh.read()
        lines=raw.split('\n')
        rea_fh.close()
        return lines
    except IOError:
        print "Error: The specified file does not exist - %s" % (filename)

def seisan_parsing(files):
    event_cat={'ID':[],'Date':[],'Time':[],'Latitude':[],'Longitude':[],
               'Depth':[],'mag':[],'mag_type':[],'mag_agency':[]}
    picks_cat={'ID':[],'Station':[],'Date':[],'Pick_Time':[],'Amplitude':[],'Period':[],'Distance':[],'Azimuth':[],'BGS_Mag':[]}
    for file in files:
        lines=read_file(file)
        for line in lines:
            if lntype(line) != '':
            # Process Type 1 line
                if lntype(line) == '1':

                    year=int(line[1:5])
                    month=int(line[6:8])
                    day=int(line[8:10])
                    hour=int(line[11:13])
                    min=int(line[13:15])
                    isec=int(line[16:18])
                    if isec == 60: # 
                        min+=1
                        isec=isec-60
                    origin=datetime(year,month,day,hour,min,isec)

                    if line[23:30]!="       ":
                        latitude=float(line[23:30])
                    if line[30:38]!="        ":
                        longitude=float(line[30:38])
                    if line[38:43]!="     ":
                        depth=float(line[38:43])

                    if line[55:59] != '    ':
                        magnitude=float(line[55:59])
                        mag_type = magtype_map[line[59]]
                        mag_agency = line[60:63]
            if lntype(line) == 'I':
                event_id=int(line[60:74])
            
            if lntype(line) == ' ' or lntype(line) == '4':
                if line[10:14] == 'IAML':
                    if line[33:40]!='       ':
#                         tmp={}
#                         #print line[1:6]
#                         sta=line[1:6]
#                         comp=line[6:8]

                        hour=int(line[18:20])
                        minute=int(line[20:22])
                        sec=float(line[22:28])
            
                        pick_time='%02d:%02d:%s'%(hour,minute,sec)
#                         print pick_time

#                         amplitude=float(line[33:40])
#                         period=float(line[41:45])
#                         distance=float(line[70:75])
#                         azimuth=int(line[76:79])
                        
                        picks_cat['ID'].append(event_id)
                        picks_cat['Station'].append(line[1:5])
                        picks_cat['Date'].append(origin.date())
                        picks_cat['Pick_Time'].append(pick_time)
                        picks_cat['Amplitude'].append(float(line[33:40]))
                        picks_cat['Period'].append(float(line[41:45]))
                        print line[70:75]
#                         try:
#                             distance = float(line[70:75])
#                             print distance
#                         except:
#                             print event_id
                        picks_cat['Distance'].append(float(line[70:75]))

                        picks_cat['Azimuth'].append(float(line[76:79]))
                        picks_cat['BGS_Mag'].append(round(magnitude,2))

                        
                    
                
                

        event_cat['ID'].append(event_id)
        event_cat['Date'].append(origin.date())
        event_cat['Time'].append(origin.time())
        event_cat['Latitude'].append(latitude)
        event_cat['Longitude'].append(longitude)
        event_cat['Depth'].append(depth)
        event_cat['mag'].append(magnitude)
        event_cat['mag_type'].append(mag_type)
        event_cat['mag_agency'].append(mag_agency)
        
        
        
    return event_cat,picks_cat
    
files = glob.glob('data/New_Ollerton/*.S2014*')
event_cat,picks_cat = seisan_parsing(files)

events=pd.DataFrame(event_cat)
events.to_csv("output/events_1001.csv",index=False,cols=['ID','Date','Time','Latitude','Longitude','Depth',
                                                         'mag','mag_type','mag_agency'])

picks=pd.DataFrame(picks_cat)
picks.to_csv("output/BGS-picks_1001.csv",index=False,cols=['ID','Station','Date','Pick_Time','Amplitude',
                                                           'Period','Distance','Azimuth','BGS_Mag'])