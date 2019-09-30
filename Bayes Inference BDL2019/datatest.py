# -*- coding: utf-8 -*-
"""
Created on Sat Sep 14 17:53:46 2019

@author: user
"""
from obspy import read
import obspy
import matplotlib.pyplot as plt
import numpy as np

st = read('http://examples.obspy.org/RJOB_061005_072159.ehz.new')
print(st)
tr = st[0]
print(tr)
fs = 200
N = len(tr)
print(tr.stats) 
data = tr.data
st.plot()

tr_filt = tr.copy()
tr_filt.filter('lowpass', freq=1.0, corners=2, zerophase=True)

# plot the raw and filtered data...
t = np.arange(0, tr.stats.npts / tr.stats.sampling_rate, tr.stats.delta)
plt.subplot(211)
plt.plot(t, tr.data, 'k')
plt.ylabel('Raw Data')
plt.subplot(212)
plt.plot(t, tr_filt.data, 'k')
plt.ylabel('Lowpassed Data')
plt.xlabel('Time [s]')
plt.title(tr.stats.starttime)
plt.show()

data = st[0].data
npts = st[0].stats.npts
samprate = st[0].stats.sampling_rate

# Filtering the Stream object
st_filt = st.copy()
st_filt.filter('bandpass', freqmin=1, freqmax=3, corners=2, zerophase=True)
# Envelope of filtered data
data_envelope = obspy.signal.filter.envelope(st_filt[0].data)

# The plot the raw and envelop
t = np.arange(0, npts / samprate, 1 / samprate)
plt.plot(t, st_filt[0].data, 'k')
plt.plot(t, data_envelope, 'k:')
plt.title(st[0].stats.starttime)
plt.ylabel('Filtered Data w/ Envelope')
plt.xlabel('Time [s]')
plt.xlim(80, 90)
plt.show()

from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy import read, read_inventory

t1 = UTCDateTime("2010-09-3T16:30:00.000")
t2 = UTCDateTime("2010-09-3T17:00:00.000")
fdsn_client = Client('IRIS')
# Fetch waveform from IRIS FDSN web service into a ObsPy stream object
# and automatically attach correct response
st = fdsn_client.get_waveforms(network='NZ', station='BFZ', location='10',
                               channel='HHZ', starttime=t1, endtime=t2,
                               attach_response=True)
# define a filter band to prevent amplifying noise during the deconvolution
pre_filt = (0.005, 0.006, 30.0, 35.0)
st.remove_response(output='DISP', pre_filt=pre_filt)

#tr = st[0]
#pre_filt = [0.001, 0.005, 10, 20]
#inv = read_inventory()
#st.remove_response(inventory=inv, pre_filt=pre_filt, output="DISP",plot=True)
st = read("/path/to/IU_ULN_00_LH1_2015-07-18T02.mseed")
tr = st[0]
inv = read_inventory("/path/to/IU_ULN_00_LH1.xml")
pre_filt = [0.001, 0.005, 10, 20]
tr.remove_response(inventory=inv, pre_filt=pre_filt, output="DISP",
                   water_level=60, plot=True)


from obspy.io.xseed import Parser
from obspy.signal import PPSD
st = read("https://examples.obspy.org/BW.KW1..EHZ.D.2011.038")
tr = st.select(id="BW.KW1..EHZ")[0]

inv = read_inventory("https://examples.obspy.org/BW_KW1.xml")
ppsd = PPSD(tr.stats, metadata=inv)
ppsd.add(st)
ppsd.plot()
ppsd.plot_temporal([0.1, 1, 10])
ppsd.plot_spectrogram()

from obspy.imaging.cm import obspy_sequential
from obspy.signal.tf_misfit import cwt


st = obspy.read()
tr = st[0]
npts = tr.stats.npts
dt = tr.stats.delta
t = np.linspace(0, dt * npts, npts)
f_min = 1
f_max = 50

scalogram = cwt(tr.data, dt, 8, f_min, f_max)

fig = plt.figure()
ax = fig.add_subplot(111)

x, y = np.meshgrid(
    t,
    np.logspace(np.log10(f_min), np.log10(f_max), scalogram.shape[0]))

ax.pcolormesh(x, y, np.abs(scalogram), cmap=obspy_sequential)
ax.set_xlabel("Time after %s [s]" % tr.stats.starttime)
ax.set_ylabel("Frequency [Hz]")
ax.set_yscale('log')
ax.set_ylim(f_min, f_max)
plt.show()

#transfer to matlab
for i, tr in enumerate(st):
    mdict = {k: str(v) for k, v in tr.stats.iteritems()}
    mdict['data'] = tr.data
    savemat("data-%d.mat" % i, mdict)