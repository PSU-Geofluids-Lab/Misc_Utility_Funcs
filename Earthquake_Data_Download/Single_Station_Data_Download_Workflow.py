#!/usr/bin/env python
# coding: utf-8

# In[4]:


import matplotlib.pyplot as plt
from obspy.taup import TauPyModel
from future.builtins import * 
from datetime import timedelta
from obspy.core import read, Stream
from obspy.core.event import read_events
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.inventory import read_inventory
from obspy.geodetics import gps2dist_azimuth, locations2degrees
import numpy as np
from obspy.clients.fdsn.client import Client
import pandas as pd
import glob
from obspy.io.sac import SACTrace
import time
import sys
from obspy.clients.fdsn.mass_downloader import CircularDomain, Restrictions, MassDownloader
import os
import time


# In[5]:


## Define All Constants
continent= sys.argv[1]           ## update the continent in which station is
data_center = sys.argv[2]        ## update data center name
network = sys.argv[3]            ## update network name
station = sys.argv[4]            ## update station name
minmagnitude = 5.5               ## update minimum nad maximum magnitude from which to download data
maxmagnitude = 7
event_t_before = 300             ## update how many seconds before arrival data needs to be downloaded
event_t_after = 300              ## update how many seconds after arrival data needs to be downloaded
epicentral_dis_min = 300         ## update minimum and maximum epicentral distance, change to 0 and 180 if not restricted
epicentral_dis_max = 90

root_path = f'datadir/data'
station_inv_path = f'datadir/inv'
logging_path = f'datadir/logging'
err_path = f'datadir/err'


# In[6]:


def get_mseed_storage(network, station, location, channel, starttime, endtime):
    filename = network + '.' + station + '.' + str(event_count) + '.' + location + '.' + channel + '.mseed'
    filepath = os.path.join(root_path, network, station, filename)
    
    if os.path.exists(filepath):
        return True

    return filepath


# In[7]:


client_data = Client(data_center)
inv = client_data.get_stations(loc="*", channel="*", network=network, station=station, level='response')
total_networks = len(inv.get_contents()['networks'])

station_catalog_list = []
for i_net in range(0, total_networks):
    total_stations = len(inv[i_net].get_contents()['stations'])
    
    for i_sta in range(0, total_stations):
        station_code = inv[i_net][i_sta].__getattribute__('code')
        station_start_date = inv[i_net][i_sta].__getattribute__('start_date')
        if station_start_date==None:
            station_start_date = inv[i_net][i_sta].__getattribute__('creation_date')
        station_end_date = inv[i_net][i_sta].__getattribute__('end_date')
        if station_start_date!=None and station_end_date==None:
            station_end_date = '2023-12-31T00:00:00.000000Z'
        station_lat = inv[i_net][i_sta].__getattribute__('latitude')
        station_lon = inv[i_net][i_sta].__getattribute__('longitude')
        station_ele = inv[i_net][i_sta].__getattribute__('elevation')
        station_restricted_status = inv[i_net][i_sta].__getattribute__('restricted_status')
        network_code = inv[i_net].__getattribute__('code')
        station_catalog = [network_code,station_code,station_start_date,
                          station_end_date,station_lat,station_lon,station_ele,station_restricted_status]
        station_catalog_list.append(station_catalog)

station_metadata = pd.DataFrame(station_catalog_list, columns = ['network','station','station_start_date',
                                                                   'station_end_date','station_lat','station_lon',
                                                                   'station_ele','station_restricted_status'])


# In[7]:


event_catalog_list = []
download_err_list = []
client_events = Client('IRIS')
mdl = MassDownloader(providers=[data_center])
model = TauPyModel(model='iasp91')
event_count = 0

for idx in station_metadata.index:
    row = station_metadata.iloc[idx]
    station_start_date = UTCDateTime(row.station_start_date)
    station_end_date  = UTCDateTime(row.station_end_date)
    network = row.network
    station = row.station
    station_lat = row.station_lat
    station_lon = row.station_lon
    station_ele = row.station_ele
    
    events_list = client_events.get_events(starttime=station_start_date, endtime=station_end_date, 
                                        minmagnitude=minmagnitude, maxmagnitude=maxmagnitude, catalog = 'NEIC PDE')
    
    for events in events_list:
        tbefore = 0
        event_time = events.origins[0].time - tbefore
        evlat = events.origins[0].latitude
        evlon = events.origins[0].longitude
        mag = events.magnitudes[0].mag
        event_name = events.event_descriptions[0].text
        event_region = events.event_descriptions[0].type
        event_type = events.event_type    
        depth = events.origins[0].depth
        if (depth is not None and (depth >= 0)):
            depth_in_km = depth/1000
            distance_details  = gps2dist_azimuth(lat1=station_lat, lon1=station_lon, lat2=evlat, lon2=evlon)
            gcarc = locations2degrees(lat1=station_lat, long1=station_lon, lat2=evlat, long2=evlon)
            azi = distance_details[1]
            baz = distance_details[2]

            if (gcarc>=epicentral_dis_min and gcarc<=epicentral_dis_max):
                arrivals = model.get_travel_times(source_depth_in_km=depth_in_km,distance_in_degree=gcarc,
                                                  phase_list=["P"])
                if len(arrivals)>0:
                    first_arrival = arrivals[0].time
                    travel_time = event_time + first_arrival                    
                    event_count+=1

                    domain = CircularDomain(latitude=evlat, longitude=evlon,
                                            minradius=epicentral_dis_min, maxradius=epicentral_dis_max)
                    
                    restrictions = Restrictions(
                                    starttime = travel_time-event_t_before,
                                    endtime = travel_time+event_t_after,
                                    network=network, station=station,
                                    minimum_length=0.95)
                    
                    try:
                        mdl.download(domain, restrictions, 
                                     mseed_storage=get_mseed_storage, stationxml_storage=station_inv_path)
                        time.sleep(0.1)
                        event_catalog_list.append([event_count,network,station,station_lat,station_lon,station_ele,station_start_date,station_end_date,
                                               evlat,evlon,mag,depth_in_km,event_name,event_region,event_type,azi,baz,gcarc,event_time,travel_time])
                    except Exception as err:
                        download_err_list.append([network,station,event_count,str(err),evlat,evlon,mag,depth_in_km,event_name,event_region,event_type,azi,baz,gcarc,event_time,travel_time])
                        continue
                        
                    ##rotation and sac
                    file_list = glob.glob(f'{root_path}/{network}/{station}/{network}.{station}.{event_count}*')
                    stn_inv_file = f'{station_inv_path}/{network}.{station}.xml'
                    channel_list = []
                    npts_list = []
                    startime_list = []
                    endtime_list = []
                    sampling_rate = []
                    stream=Stream()
                    
                    if os.path.exists(stn_inv_file) and len(file_list)==3:
                        for file in file_list:
                            trace = read(file)
                            channel_list.append(trace[0].stats.channel)
                            npts_list.append(trace[0].stats.npts)
                            startime_list.append(trace[0].stats.starttime)
                            endtime_list.append(trace[0].stats.endtime)
                            sampling_rate.append(trace[0].stats.sampling_rate)
                            stream+=trace
                            
                        channel_list_str = ''.join(channel_list)
                        startime_list_str = [str(stm) for stm in startime_list]
                        endtime_list_str = [str(etm) for etm in endtime_list]
                        if len(set(sampling_rate))==1 and (len(set(startime_list_str))>1 or len(set(endtime_list_str))>1):
                            try: 
                                stream.interpolate(sampling_rate=sampling_rate[0],npts=int(600*sampling_rate[0]),
                                                   starttime=max(startime_list)+0.5, 
                                                   method='lanczos', a=10)
                            except Exception as err:
                                download_err_list.append([network,station,event_count,str(err),evlat,evlon,mag,depth_in_km,event_name,event_region,event_type,azi,baz,gcarc,event_time,travel_time])
                                continue                                
                        
                        try:
                            if channel_list_str.find('1')>=0:
                                stn_inv_data = inv.select(location=trace[0].stats.location, time=travel_time, sampling_rate=trace[0].stats.sampling_rate)
                                stream_updated = stream.select(location=trace[0].stats.location, channel='*Z', sampling_rate=trace[0].stats.sampling_rate)+ \
                                                    stream.select(location=trace[0].stats.location, channel='*1', sampling_rate=trace[0].stats.sampling_rate)+ \
                                                    stream.select(location=trace[0].stats.location, channel='*2', sampling_rate=trace[0].stats.sampling_rate)

                                stream_updated._rotate_to_zne(stn_inv_data)
                                stream_updated.rotate(method='NE->RT', back_azimuth=azi)
                            else:
                                stream_updated = stream.select(location=trace[0].stats.location, channel='*N', sampling_rate=trace[0].stats.sampling_rate)+ \
                                                    stream.select(location=trace[0].stats.location, channel='*E', sampling_rate=trace[0].stats.sampling_rate)+ \
                                                    stream.select(location=trace[0].stats.location, channel='*Z', sampling_rate=trace[0].stats.sampling_rate)

                                stream_updated.rotate(method='NE->RT', back_azimuth=azi)
                        except Exception as err:
                            download_err_list.append([network,station,event_count,str(err),evlat,evlon,mag,depth_in_km,event_name,event_region,event_type,azi,baz,gcarc,event_time,travel_time])
                            continue
                        
                        for tr in stream_updated:
                            sac = SACTrace.from_obspy_trace(tr)
                            sac.stla = station_lat
                            sac.stlo = station_lon
                            sac.o = event_time
                            sac.t0 = travel_time
                            sac.evla = evlat
                            sac.evlo = evlon
                            sac.evdp = depth_in_km
                            sac.mag = mag
                            sac.baz = azi
                            sac.gcarc = gcarc

                            kcmpnm = sac.kcmpnm
                            channel_name = kcmpnm[:-1]
                            kcmpnm = kcmpnm[-1].upper()
                            
                            sac_out = f'{root_path}/{network}/{station}/{network}.{station}.{event_count}.{trace[0].stats.location}.{kcmpnm}'
                            sac.write(sac_out)
                        
                    else:
                        download_err_list.append([network,station,event_count,'No Three Channels Or No Inventory',
                                                  evlat,evlon,mag,depth_in_km,event_name,event_region,event_type,azi,baz,gcarc,event_time,travel_time])

                        
pd.DataFrame(event_catalog_list, 
            columns=['event_id','network','station','station_lat','station_lon','station_ele',
                       'station_start_date','station_end_date',
                       'event_lat','event_lon','event_magnitude','event_depth_km','event_name','event_region',
                       'event_type','azimuth','backazimuth','event_gcarc','event_origin','event_travel_time_P'
                    ]).to_csv(f'{logging_path}/{network}_{station}_event_catalog_list.csv',index=None)

pd.DataFrame(download_err_list, 
            columns=['network','station','event_id','error_message',
                       'event_lat','event_lon','event_magnitude','event_depth_km','event_name','event_region',
                       'event_type','azimuth','backazimuth','event_gcarc','event_origin','event_travel_time_P'
                    ]).to_csv(f'{err_path}/{network}_{station}_download_err_list.csv',index=None)


# In[ ]:




