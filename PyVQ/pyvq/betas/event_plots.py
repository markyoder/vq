#!/usr/bin/env python

# ---------------------- IMPORTS -------------------
from __future__ import print_function
import math
import sys
import argparse
import quakelib
import gc
import operator
scipy_available = True
#
try:
    import scipy.stats
except ImportError:
    scipy_available = False
matplotlib_available = True
try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
except ImportError:
    matplotlib_available = False
numpy_available = True
try:
    import numpy as np
    import numpy		# i often use numpy as just numpy, so can we keep the alias?
except ImportError:
    numpy_available = False
h5py_available = True
try:
    import h5py
except ImportError:
    h5py_available = False
# ----------------------         -------------------

# ======= h5py I/O ============================================
def read_events_h5(sim_file, event_numbers=None):
    # TODO: Add event filters
    with h5py.File(sim_file) as vq_data:
        events = vq_data['events'][()]
    # If event_numbers specified, only return those events
    if event_numbers is not None:
        # Handle single events separately
        if isinstance(event_numbers, int): 
            events = np.core.records.fromarrays(zip(*filter(lambda x: x['event_number'] == event_numbers, events)), dtype=events.dtype)
        else:
            events = np.core.records.fromarrays(zip(*filter(lambda x: x['event_number'] in event_numbers, events)), dtype=events.dtype)	
	return events

def read_sweeps_h5(sim_file, event_number=0, block_ids=None):
	# Read sweeps sequence for multiple blocks (unless block_id specified) in a single event.
	with h5py.File(sim_file) as vq_data:
		sweep_range = [vq_data['events'][event_number]['start_sweep_rec'],
                       vq_data['events'][event_number]['end_sweep_rec']]
		sweeps = vq_data['sweeps'][sweep_range[0]:sweep_range[1]][()]
	# If block_id specified, only return those sweeps for that block
	if block_ids is not None:
		d_type = sweeps.dtype
		sweeps = np.core.records.fromarrays(zip(*filter(lambda x: x['block_id'] in block_ids, sweeps)), dtype=d_type)	
	return sweeps

def parse_sweeps_h5(sim_file=None, block_id=None, event_number=0, do_print=True, sweeps=None):
    # Read sweep data if not provided
	if sweeps is None: sweeps = read_sweeps_h5(sim_file, block_id=block_id, event_number=event_number)
	# Grab data
	data = [[rw['sweep_number'], rw['block_id'], rw['block_slip'], rw['shear_init'],
             rw['shear_final'], rw['normal_init'],rw['normal_final'], 
             (rw['shear_final']-rw['shear_init'])/rw['shear_init'], 
             (rw['normal_final']-rw['normal_init'])/rw['normal_init']] for rw in sweeps]
	if do_print:
		for rw in data: print(rw)
	cols = ['sweep_number', 'block_id', 'block_slip', 'shear_init', 
            'shear_final', 'normal_init', 'normal_final', 'shear_change', 'normal_change']
	return np.core.records.fromarrays(zip(*data), names=cols, formats = [type(x).__name__ for x in data[0]])


def hist_file(file_in=None, hist_cols=[0], comment_char='#', delim='\t', **hist_kwargs):
	'''
	# simple wrapper to plot histogram(s) from csv type data.
	# plot a histogram for each col in hist_cols. 
	# pass any histogram arguemtns in hist_kwargs.
	'''
	if hist_cols==None:
		# assume one histogram for each column...
		with open(file_in, 'r') as f:
			for rw in f:
				if rw[0]==comment_char: continue
				hist_cols = range(len(rw.split()))
				break
			#
		#
	#
	
	#
	with open(file_in, 'r') as f:			
		datas = zip(*[rw.split() for rw in f if rw[0]!=comment_char[0]])
	#
	hist_kwargs['bins']=hist_kwargs.get('bins', len(datas[0])/100)
	#
	for fnum,col in enumerate(hist_cols):
		plt.figure(fnum)
		plt.clf()
		#
		plt.hist([float(x) for x in datas[col]], **hist_kwargs)
	

# ======= SIM DATA CLASSES ===========================================
class Events(object):
    def __init__(self, sim_file, min_detectable_slip=.001):
    	# block_size should probably not be a parameter; we should read it from a corresponding fault geometry... won't really matter at first,
    	# since it won't change scaling.
        # TODO: Add event filters
        #self.events = read_events_h5(sim_file)
        #
        print("Reading events from {}".format(sim_file))
        #
        # now, derived values from sweeps: area, slip, ??
        #
        new_cols = ['mean_slip', 'area', 'area_prime']
        events_2 = []	# we'll wrap these together with events later
        #
        # in the event that we encounter large sweeps files, load sweep data one (or a few -- we'll get to that later) at a time.
        with h5py.File(sim_file) as vc_data:
        	self.events = vc_data['events'][()]
        	#return self.events
        	
        	sweeps=vc_data['sweeps']
        	for j, rw in enumerate(self.events):
        		#"event-sweep"
        		
        		#j0=
        		#j1=rw
        		#print("range: ", j0,j1)
        		ev_sw = sweeps[rw['start_sweep_rec']:rw['end_sweep_rec']][()]
        		#n_elements = len(set(ev_sw['block_id']))
        		n_elements = rw['end_sweep_rec'] - rw['start_sweep_rec']
        		#block_areas = list(set([[rw['block_id'], rw['block_area']] for rw in ev_sw]))
        		block_areas_unique = {key:val for key,val in [[rw['block_id'], rw['block_area']] for rw in ev_sw]}		# this is a bit more stable, in the event that there are small
        																										# maybe floating point differences in block area. taking the mean
        																										# area would probabu be better...
        		area_total = numpy.sum(block_areas_unique.values())
        		mean_slip=numpy.sum([rw2['block_slip']*rw2['block_area'] for rw2 in ev_sw])/area_total
        		
        		block_slips = {b_id:0. for b_id in ev_sw['block_id']}
        		for rw in ev_sw: block_slips[rw['block_id']]+=rw['block_slip']
        		#area_slip_weighted = numpy.sum([block_areas_unique[key]*block_slips[key] for key in block_slips.iterkeys()])/mean_slip
        		area_prime = sum([val for key,val in block_areas_unique.iteritems() if abs(block_slips[key]) > min_detectable_slip])
        		#
        		events_2 += [[mean_slip, area_total, area_prime]]
        	#
        #
        #self.events_2=events_2
        self.events = np.core.records.fromarrays(list(zip(*self.events)) + list(zip(*events_2)), names = list(self.events.dtype.names) + new_cols, formats = [rw[1] for rw in self.events.dtype.descr] + [type(x).__name__ for x in events_2[0]])
        #
    def plot_slip_mag(self, fignum=0, event_ids=None, block_ids=None):
        # TODO: eventually allow event, block_id filters.
        plt.figure(fignum)
        plt.clf()
        ax=plt.gca()
        ax.set_yscale('log')
        #
        try:
            # intercept, slope
            lS = numpy.log10(numpy.abs(self.events['mean_slip']))
            a,b = scipy.optimize.curve_fit(lambda x,a,b: a + b*x, self.events['event_magnitude'], lS)
            print ("fits: %f, %f" % (a,b))
            inv_log = lambda x: 10.**(a + b*x)
            X = [[0], self.events['event_magnitude'][-1]]
            #ax.plot(X, [inv_log(x) for x in X], '.-', lw=2, label='scaling: a=%.3f, b=%.3f' % (a,b))
            ax.plot(X, [10.**(a+b*x) for x in X], '.-', lw=2, label='scaling: a=%.3f, b=%.3f' % (a,b))
        except:
        	print("fit to data failed...")
        #
        #ax.plot(self.events['event_magnitude'], numpy.abs(self.events['mean_slip']), '.')
        ax.plot(*zip(*[[m,abs(s)] for m,s in zip(self.events['event_magnitude'], self.events['mean_slip']) if s<0.]), marker='.', ls='', label='neg')
        ax.plot(*zip(*[[m,abs(s)] for m,s in zip(self.events['event_magnitude'], self.events['mean_slip']) if s>0.]), marker='.', ls='', label='pos')
        
        plt.title('Slip-magnitude scaling')
        plt.xlabel('Event magnitude $m$', size=18)
        plt.ylabel('Event slip $s$', size=18)
        plt.legend(loc=0, numpoints=1)
    
    def plot_area_mag(self, fignum=0, event_ids=None, block_ids=None, min_slip=.001):
    	# @min_slip: minimum slip on block to consider for area calc; if slip<min_slip, exclude from area calculation.
        # TODO: eventually allow event, block_id filters.
        plt.figure(fignum)
        plt.clf()
        ax=plt.gca()
        ax.set_yscale('log')
        #
        try:
            # fit these data to get a scaling expon/slope:
            #
            a_0, b_0 = scipy.optimize.curve_fit(lambda x,a,b: a+b*x, self.events['event_magnitude'], numpy.log10(self.events['area']))[0]
            a_1, b_1 = scipy.optimize.curve_fit(lambda x,a,b: a+b*x, self.events['event_magnitude'], numpy.log10(self.events['area_prime']))[0]
            print ("fits_0: %f, %f" % (a_0,b_0))
            print ("fits_1: %f, %f" % (a_1,b_1))
            #
            X = [self.events['event_magnitude'][0], self.events['event_magnitude'][-1]]
            inv_log = lambda x: 10.**(a_0 + b_0*x)
            ax.plot(X, [inv_log(x) for x in X], '.-', lw=2, label='full_area: a=%.3f, b=%.3f' % (a_0,b_0), color='b')
            inv_log = lambda x: 10.**(a_1 + b_1*x)
            ax.plot(X, [inv_log(x) for x in X], '.-', lw=2, label='full_area: a=%.3f, b=%.3f' % (a_1,b_1), colog='g')
        except:
			print("fit to area scaling data failed.")
        
        
        #
        ax.plot(self.events['event_magnitude'], self.events['area'], '.', color='b')
        ax.plot(self.events['event_magnitude'], self.events['area_prime'], '.', color='g')
        
        
        
        #area_prime = [[rw['event_magnitude'], rw['area']]
        
        plt.title('Area-magnitude scaling')
        plt.xlabel('Event magnitude $m$', size=18)
        plt.ylabel('Event area $A$', size=18)


class Sweeps:
    def __init__(self, sim_file, event_number=0, block_ids=None):
        self.sweeps = read_sweeps_h5(sim_file, event_number=event_number, block_ids=block_ids)
        self.sweep_data = parse_sweeps_h5(sweeps=self.sweeps, do_print=False)					# not sure what the difference is here...
        self.block_ids = self.sweep_data['block_id'].tolist()
        self.mag = read_events_h5(sim_file,event_numbers=event_number)['event_magnitude'][0]
        self.event_number = event_number
        print("Read event {} sweeps from {}".format(event_number,sim_file))
        # we could also, at this point, parse out the individual block sequences, maybe make a class Block().
    #
    def block_slips(self):
        # return total slip on each block.
        slips = {rw['block_id']:0. for rw in self.sweeps}
        for rw in self.sweeps: slips[rw['block_id']]+=rw['block_slip']
        #
        return numpy.core.records.fromarrays(zip(*[[key,val] for key,val in slips.iteritems()]), names=('block_id','slip'), formats=('u8', 'f8'))
    #
    def plot_event_block_slips(self, block_ids=None, fignum=0):
        block_ids = self.check_block_ids_list(block_ids)
        plt.figure(fignum)
        plt.clf()
        for block_id in block_ids:
            rws = np.core.records.fromarrays(zip(*filter(lambda x: x['block_id']==block_id, self.sweep_data)), dtype=self.sweep_data.dtype)
            #plt.semilogy(rws['sweep_number'], rws['block_slip'], '.-', label='blk: %d' % block_id)
            plt.plot(rws['sweep_number'], rws['block_slip'], '.-', label='blk: %d' % block_id)
        if len(block_ids) <= 10:
            plt.legend(loc='best', numpoints=1,fontsize=8,ncol=3,handlelength=2,handletextpad=1)
        plt.title('Event {} (M={:.2f}) slips for {} blocks'.format(self.event_number,self.mag,len(block_ids)))
        plt.xlabel('sweep number')
        plt.ylabel('slip [m]')
        min_sweep = 0
        max_sweep = int(max(self.sweep_data['sweep_number']))
        if max(self.sweep_data['sweep_number']) < 3:
            max_sweep += 1
        ticks = range(max_sweep+1)
        plt.xticks(ticks,[str(tick) for tick in ticks])
        plt.xlim(min_sweep, max_sweep)
    #
    def plot_stress_changes(self, block_ids=None, fignum=0, shear=True,log=False,max_val=None):
        block_ids = self.check_block_ids_list(block_ids)
        #
        plt.figure(fignum)
        plt.clf()
        ax = plt.gca()
        #
        # we should move all plotting parameters into some sort of **kwargs argument. for now, go with what we've got... but clean up a little bit.
        ax.set_yscale('log' if log else 'linear')
        #
        for block_id in block_ids:
            rws = np.core.records.fromarrays(zip(*filter(lambda x: x['block_id']==block_id, self.sweep_data)), dtype=self.sweep_data.dtype)
            if shear: 
                plt.plot(rws['sweep_number'], rws['shear_change'], '.-', label=block_id)
                #if not log:
                #    plt.plot(rws['sweep_number'], rws['shear_change'], '.-', label=block_id)
                #else:
                #    plt.semilogy(rws['sweep_number'], rws['shear_change'], '.-', label=block_id)
            else: 
                plt.plot(rws['sweep_number'], rws['shear_change'], '.-', label=block_id)
                #if not log:
                #    plt.plot(rws['sweep_number'], rws['shear_change'], '.-', label=block_id)
                #else:
                #    plt.semilogy(rws['sweep_number'], rws['shear_change'], '.-', label=block_id)
        #
		plt.plot([min(self.sweep_data['sweep_number']), max(self.sweep_data['sweep_number'])], [0., 0.], 'k-')
        if len(block_ids) <= 10:
            plt.legend(loc='best', numpoints=1,fontsize=8,ncol=3,handlelength=2,handletextpad=1)
        if shear: 
            plt.title('Event {} (M={:.2f}) shear stress changes for {} blocks'.format(self.event_number,self.mag,len(block_ids)))
        else: 
            plt.title('Event {} (M={:.2f}) normal stress changes for {} blocks'.format(self.event_number,self.mag,len(block_ids)))
        plt.xlabel('sweep number')
        plt.ylabel('fractional stress change')
        min_sweep = 0
        max_sweep = int(max(self.sweep_data['sweep_number']))
        if max(self.sweep_data['sweep_number']) < 3:
            max_sweep += 1
        ticks = range(max_sweep+1)
        plt.xticks(ticks,[str(tick) for tick in ticks])
        plt.xlim(min_sweep, max_sweep)
        if max_val is not None: plt.ylim(-max_val,max_val)
    #    
    def check_block_ids_list(self, block_ids):
        # Make sure the block_ids are a list
        if block_ids is None: block_ids=self.block_ids
        if isinstance(block_ids, float): block_ids=[int(block_ids)]
        if isinstance(block_ids, int): block_ids = [block_ids]
        return block_ids

if __name__=='__main__':
	# probably use 'Agg' backend for mpl...
	pass
else:
	plt.ion()

# ============================ TEMP. RUNNER ===================
# Example usage below
"""
SIM_FILE = "../../../../Desktop/RUNNING/UCERF2/events_ALLCAL2_VQmeshed_3km_EQSim_StressDrops_1600yr_22June2015.h5"
EVENT_NUM = 1541 #13, 948, 504, 1541
BLOCK_IDS = None
sim_sweeps = Sweeps(SIM_FILE, event_number=EVENT_NUM, block_ids=BLOCK_IDS)
# ---- plot slips ---------

sim_sweeps.plot_event_block_slips()
savename = "../../../../VQScripts/event_{}_slips{}.png".format(sim_sweeps.event_number,SIM_FILE.split("/")[-1].split(".")[0].split("events")[-1])
plt.savefig(savename,dpi=100)

# ---- plot stresses ---------
sim_sweeps.plot_stress_changes()
savename = "../../../../VQScripts/event_{}_shear_changes{}.png".format(sim_sweeps.event_number,SIM_FILE.split("/")[-1].split(".")[0].split("events")[-1])
plt.savefig(savename,dpi=100)
"""
