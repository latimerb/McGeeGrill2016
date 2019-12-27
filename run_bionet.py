import os, sys
from bmtk.simulator import bionet
from bmtk.simulator.bionet.default_setters.cell_models import loadHOC
from bmtk.simulator.bionet.modules.sim_module import SimulatorMod
from bmtk.utils.reports.spike_trains import SpikeTrains
from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator
import numpy as np
import matplotlib.pyplot as plt
from neuron import h
import pdb

np.random.seed(2315)

# Import the synaptic depression/facilitation model
import synapses
synapses.load()

# Create lists for firing rates, bladder volumes, and bladder pressures
times = []
b_vols = []
b_pres = []

press_thres = 40 # cm H20
                 # Lingala, et al. 2016
change_thres = 10 # cm H20 
                 # Need biological value for this

bionet.pyfunction_cache.add_cell_model(loadHOC, directive='hoc', model_type='biophysical')

# Huge thank you to Kael Dai, Allen Institute 2019 for the template code
# we used to create the feedback loop below
class FeedbackLoop(SimulatorMod):
    def __init__(self):
        self._synapses = {}
        self._netcon = None
        self._spike_events = None

        self._high_level_neurons = []
        self._low_level_neurons = []

        self._synapses = {}
        self._netcons = {}

        self._spike_records = {}
        self._glob_press = 0
        self._prev_glob_press = 0

    def _set_spike_detector(self, sim):
        for gid in self._low_level_neurons:
            tvec = h.Vector()
            nc = self._netcons[gid]
            nc.record(tvec)
            self._spike_records[gid] = tvec

    def _activate_hln(self, sim, block_interval, firing_rate):
        next_block_tstart = (block_interval[1] + 1) * sim.dt  # The next time-step
        next_block_tstop = next_block_tstart + sim.nsteps_block*sim.dt  # The time-step when the next block ends

        # This is where you can use the firing-rate of the low-level neurons to generate a set of spike times for the
        # next block
        psg = PoissonSpikeGenerator()
        psg.add(node_ids=[0], firing_rate=firing_rate/1000, times=(next_block_tstart, next_block_tstop))
        self._spike_events = psg.get_times(0)
        
        if self._spike_events.shape[0]==0:
            self._spike_events = np.array((next_block_tstart+0.2,))

        # Update firing rate of bladder afferent neurons
        for gid in self._high_level_neurons:
            nc = self._netcons[gid]
            for t in self._spike_events:
                nc.event(t)

        # If pressure is maxxed, update firing rate of EUS motor neurons 
        # Guarding reflex
        #press_change = self._prev_glob_press - self._glob_press
        #if self._glob_press > press_thres or press_change > change_thres:
        #    psg = PoissonSpikeGenerator()
        #    eus_fr = self._glob_press*10 + press_change*10 # Assumption: guarding reflex
        #                                                   # depends on current pressure
        #                                                   # and change from last pressure
        #    psg.add(node_ids=[0], firing_rate=eus_fr, times=(next_block_tstart, next_block_tstop))
        #    self._spike_events = psg.get_times(0)
        #    for gid in self._spn_neurons:
        #        nc = self._netcons[gid]
        #        for t in self._spike_events:
        #            nc.event(t)

    def initialize(self, sim):
        # Attach a NetCon/synapse on the high-level neuron(s) soma. We can use the NetCon.event(time) method to send
        # a spike to the synapse. Which, is a spike occurs, the high-level neuron will inhibit the low-level neuron.
        spikes = h.Vector([])  # start off with empty input
        vec_stim = h.VecStim()
        vec_stim.play(spikes)
        self._high_level_neurons = list(sim.net.get_node_set('high_level_neurons').gids())
        for gid in self._high_level_neurons:
            cell = sim.net.get_cell_gid(gid)
            # Create synapse
            # These values will determine how the high-level neuron behaves with the input
            new_syn = h.Exp2Syn(0.5, cell.hobj.soma[0])
            new_syn.e = 0.0
            new_syn.tau1 = 1.0
            new_syn.tau2 = 8.0

            nc = h.NetCon(vec_stim, new_syn)
            nc.weight[0] = 0.09
            nc.delay = 1.0
            nc.threshold = -50

            self._synapses[gid] = new_syn
            self._netcons[gid] = nc


        # Attach another netcon to the low-level neuron(s) that will record
        self._low_level_neurons = list(sim.net.get_node_set('low_level_neurons').gids())
        for gid in self._low_level_neurons:
            cell = sim.net.get_cell_gid(gid)
            # NOTE: If you set the segment to 0.5 it will interfere with the record_spikes module. Neuron won't let
            # us record from the same place multiple times, so for now just set to 0.0. TODO: read and save spikes in biosimulator.
            nc = h.NetCon(cell.hobj.soma[0](0.0)._ref_v, None, sec=cell.hobj.soma[0])
            nc.threshold = -50
            self._netcons[gid] = nc

        self._set_spike_detector(sim)

    def step(self, sim, tstep):
        pass

    def block(self, sim, block_interval):
        """This function is called every n steps during the simulation, as set in the config.json file (run/nsteps_block).

        We can use this to get the firing rate of PGN during the last block and use it to calculate
        firing rate for bladder afferent neuron
        """

        # Calculate the avg number of spikes per neuron
        block_length = sim.nsteps_block*sim.dt/1000.0  #  time length of previous block of simulation TODO: precalcualte
        n_gids = 0
        n_spikes = 0
        for gid, tvec in self._spike_records.items():
            n_gids += 1
            n_spikes += len(list(tvec))  # use tvec generator. Calling this deletes the values in tvec
        avg_spikes = n_spikes/(float(n_gids)*3.0)

        # Calculate the firing rate the the low-level neuron(s)
        fr = avg_spikes/float(block_length)
        print("SPN firing rate: {}".format(fr))

        # Grill function for polynomial fit according to PGN firing rate
	    # Grill, et al. 2016
        def pgn_fire_rate(x):
            f = 2.0E-03*x**3 - 3.3E-02*x**2 + 1.8*x - 0.5
            
            # Take absolute value of firing rate if negative
            if f < 0:
                f *= -1

            return f

        # Grill function for polynomial fit according to bladder volume
	    # Grill, et al. 2016
        def blad_vol(vol):
            f = 1.5*20*vol - 10

            return f

        # Grill function returning pressure in units of cm H20
	    # Grill, et al. 2016
        def pressure(fr,v):
            p = fr + v

            # Round negative pressure up to 0
            if p < 0:
                p = 0

            return p 

        # Grill function returning bladder afferent firing rate in units of Hz
	    # Grill, et al. 2016
        def blad_aff_fr(p):
            fr = -3.0E-08*p**5 + 1.0E-5*p**4 - 1.5E-03*p**3 + 7.9E-02*p**2 - 0.6*p
			
            # Take absolute value of negative firing rates
            if fr < 0:
                fr *= -1 

            return fr # Using scaling factor of 5 here to get the correct firing rate range 

        # Calculate bladder volume using Grill's polynomial fit equation
        v_init = 0.33       # TODO: get biological value for initial bladder volume
        fill = 0.05 	 	# ml/min (Asselt et al. 2017)
        fill /= (1000 * 60) # Scale from ml/min to ml/ms
        void = 4.6 	 		# ml/min (Streng et al. 2002)
        void /= (1000 * 60) # Scale from ml/min to ml/ms
        max_v = 0.76 		# ml (Grill et al. 2019)
        vol = v_init

        # Update blad aff firing rate
        if len(tvec) > 0:
            PGN_fr = pgn_fire_rate(fr)

            # Filling: 0 - 50000 ms
            if tvec[0] < 50000 and vol < max_v:
                vol = fill*tvec[0]*150 + v_init
                
            # Voiding: 50000 - 60,000 ms
            else:
                vol = max_v - void*(10000-tvec[0])*100
            print("volume: {}".format(vol))
            # Maintain minimum volume
            if vol < v_init:
                vol = v_init
            grill_vol = blad_vol(vol)

            # Calculate pressure using Grill equation
            p = pressure(PGN_fr, grill_vol)
            print("pressure: {}".format(p))

            # Update global pressure (to be used to determine if EUS motor
            # needs to be updated for guarding reflex)
            self._prev_glob_press = self._glob_press
            self._glob_press = p 

            # Calculate bladder afferent firing rate using Grill equation
            bladaff_fr = blad_aff_fr(p)

            print('PGN firing rate = %.2f Hz' %fr)
            print('Volume = %.2f ml' %vol)
            print('Pressure = %.2f cm H20' %p)
            print('Bladder afferent firing rate = %.2f Hz' %bladaff_fr)

            # Save values in appropriate lists
            times.append(tvec[0])
            b_vols.append(vol)
            b_pres.append(p)
            
            # Set the activity of high-level neuron
            self._activate_hln(sim, block_interval, bladaff_fr)

        else:
            self._activate_hln(sim, block_interval, fr) 
            print('Volume = %.2f ml' %vol)
            #print('Pressure = %.2f cm H20' %p)

        # NEURON requires resetting NetCon.record() every time we read the tvec.
        self._set_spike_detector(sim)

    def finalize(self, sim):
        pass

def run(config_file):
    conf = bionet.Config.from_json(config_file, validate=True)
    conf.build_env()

    fbmod = FeedbackLoop()

    graph = bionet.BioNetwork.from_config(conf)
    sim = bionet.BioSimulator.from_config(conf, network=graph)
    sim.add_mod(fbmod)  # Attach the above module to the simulator.
    sim.run()

    ############ Plotting ################
#    # Plot bladder afferent firing rate
#    spike_trains = SpikeTrains.from_sonata('output/spikes.h5')
#    
#    ba_means = np.zeros(sim.n_steps)
#    ba_stdevs = np.zeros(sim.n_steps)
#    ba_fr_conv = np.zeros((numBladaff,sim.n_steps))
#    plt.figure()
#    
#    for gid in np.arange(0, numBladaff):
#        spikes = np.zeros(sim.n_steps, dtype=np.int)
#        if len(spike_trains.get_times(gid)) > 0:
#            spikes[(spike_trains.get_times(gid)/sim.dt).astype(np.int)] = 1
#        window_size = 10000
#        window = np.zeros(window_size)
#        window[0:window_size] = 1
#
#        frs = np.convolve(spikes, window)
#
#        for n in range(len(ba_means)):
#            ba_means[n] += frs[n]
#            ba_fr_conv[gid][n] = frs[n]
#
#    # Calculate mean and standard deviation for all points
#    for n in range(len(ba_means)):
#        ba_means[n] /= numBladaff*2
#        ba_stdevs[n] = np.std(ba_fr_conv[:,n])
#
#    # Only plot one point each 1000 samples
#    plt_ba_means = []
#    plt_ba_stdevs = []
#    for n in np.arange(0,len(ba_means), 1000):
#        plt_ba_means.append(ba_means[n])   
#        plt_ba_stdevs.append(ba_stdevs[n]) 
#
#    plt.plot(np.arange(0,len(ba_means)/10,100), plt_ba_means, 
#                 color='b', marker='^', mfc='b', mec='b', label='Bladder Afferent')
#    plt.xlabel('Time (t) [ms]')
#
#    # Plot PGN firing rate
#    pgn_means = np.zeros(sim.n_steps)
#    pgn_stdevs = np.zeros(sim.n_steps)
#    pgn_fr_conv = np.zeros((numPGN,sim.n_steps))
#
#    for gid in np.arange(PGN_gids, PGN_gids + numPGN):
#        spikes = np.zeros(sim.n_steps, dtype=np.int)
#        if len(spike_trains.get_times(gid)) > 0:
#            spikes[(spike_trains.get_times(gid)/sim.dt).astype(np.int)] = 1
#        window_size = 10000
#        window = np.zeros(window_size)
#        window[0:window_size] = 1
#
#        frs = np.convolve(spikes, window)
#
#        for n in range(len(pgn_means)):
#            pgn_means[n] += frs[n]
#            pgn_fr_conv[gid % PGN_gids][n] = frs[n]
#
#    # Calculate mean and standard deviation for all points
#    for n in range(len(pgn_means)):
#        pgn_means[n] /= numPGN*2
#        pgn_stdevs[n] = np.std(pgn_fr_conv[:,n])
#
#    # Only plot one point each 1000 samples
#    plt_pgn_means = []
#    plt_pgn_stdevs = []
#    for n in np.arange(0,len(pgn_means),1000):
#        plt_pgn_means.append(pgn_means[n])   
#        plt_pgn_stdevs.append(pgn_stdevs[n]) 
#
#    plt.plot(np.arange(0,len(pgn_means)/10,100), plt_pgn_means,  
#                 color='g', marker='o', mfc='g', mec='g', label='PGN')
#
#    # Plot EUS motor neuron firing rate
#    eus_means = np.zeros(sim.n_steps)
#    eus_stdevs = np.zeros(sim.n_steps)
#    eus_fr_conv = np.zeros((numEUSmn,sim.n_steps))
#
#    for gid in np.arange(EUSmn_gids, EUSmn_gids + numEUSmn):
#        spikes = np.zeros(sim.n_steps, dtype=np.int)
#        if len(spike_trains.get_times(gid)) > 0:
#            spikes[(spike_trains.get_times(gid)/sim.dt).astype(np.int)] = 1
#        window_size = 10000
#        window = np.zeros(window_size)
#        window[0:window_size] = 1
#
#        frs = np.convolve(spikes, window)
#
#        for n in range(len(eus_means)):
#            eus_means[n] += frs[n]
#            eus_fr_conv[gid % EUSmn_gids][n] = frs[n]
#
#    # Calculate mean and standard deviation for all points
#    for n in range(len(eus_means)):
#        eus_means[n] /= numEUSmn*2
#        eus_stdevs[n] = np.std(eus_fr_conv[:,n])
#
#    # Only plot one point each 1000 samples
#    plt_eus_means = []
#    plt_eus_stdevs = []
#    for n in np.arange(0,len(eus_means),1000):
#        plt_eus_means.append(eus_means[n])   
#        plt_eus_stdevs.append(eus_stdevs[n]) 
#
#    plt.plot(np.arange(0,len(eus_means)/10,100), plt_eus_means, 
#                 color='k', marker='D', mfc='k', mec='k', label='EUS Motor Neurons')
#    plt.xlabel('Time (t) [ms]')
#    plt.ylabel('Neuron Firing Rate (FR) [Hz]')
#    plt.legend()
#
#    # Plot bladder volume and bladder pressure
#    fig1, ax1_1 = plt.subplots()
#
#    color = 'tab:red'
#    ax1_1.set_xlabel('Time (t) [ms]')
#    ax1_1.set_ylabel('Bladder Volume (V) [ml]', color=color)
#    ax1_1.plot(times, b_vols, color=color)
#    ax1_1.tick_params(axis='y', labelcolor=color)
#
#    ax2_1 = ax1_1.twinx()  # instantiate a second axes that shares the same x-axis
#
#    color = 'tab:blue'
#    ax2_1.set_ylabel('Bladder Pressure (P) [cm H20]', color=color)  # we already handled the x-label with ax1
#    ax2_1.plot(times, b_pres, color=color)
#    ax2_1.tick_params(axis='y', labelcolor=color)
#
#    fig1.tight_layout()  # otherwise the right y-label is slightly clipped
#
#    # Plot interneurons and EUS afferent firing rates
#    inmm_means = np.zeros(sim.n_steps)
#    inmm_stdevs = np.zeros(sim.n_steps)
#    inmm_fr_conv = np.zeros((numINmminus,sim.n_steps))
#    plt.figure()
#    
#    for gid in np.arange(INmminus_gids, INmminus_gids + numINmminus):
#        spikes = np.zeros(sim.n_steps, dtype=np.int)
#        if len(spike_trains.get_times(gid)) > 0:
#            spikes[(spike_trains.get_times(gid)/sim.dt).astype(np.int)] = 1
#        window_size = 10000
#        window = np.zeros(window_size)
#        window[0:window_size] = 1
#
#        frs = np.convolve(spikes, window)
#
#        for n in range(len(ba_means)):
#            inmm_means[n] += frs[n]
#            inmm_fr_conv[gid % INmminus_gids][n] = frs[n]
#
#    # Calculate mean and standard deviation for all points
#    for n in range(len(ba_means)):
#        inmm_means[n] /= numINmminus
#        inmm_stdevs[n] = np.std(inmm_fr_conv[:,n])
#
#    # Only plot one point each 1000 samples
#    plt_inmm_means = []
#    plt_inmm_stdevs = []
#    for n in np.arange(0,len(inmm_means), 1000):
#        plt_inmm_means.append(inmm_means[n])   
#        plt_inmm_stdevs.append(inmm_stdevs[n]) 
#
#    plt.plot(np.arange(0,len(inmm_means)/10,100), plt_inmm_means, 
#                 color='b', marker='^', mfc='b', mec='b', label='INm-')
#    plt.xlabel('Time (t) [ms]')
#
#    inmp_means = np.zeros(sim.n_steps)
#    inmp_stdevs = np.zeros(sim.n_steps)
#    inmp_fr_conv = np.zeros((numINmplus,sim.n_steps))
#    
#    for gid in np.arange(INmplus_gids, INmplus_gids + numINmplus):
#        spikes = np.zeros(sim.n_steps, dtype=np.int)
#        if len(spike_trains.get_times(gid)) > 0:
#            spikes[(spike_trains.get_times(gid)/sim.dt).astype(np.int)] = 1
#        window_size = 10000
#        window = np.zeros(window_size)
#        window[0:window_size] = 1
#
#        frs = np.convolve(spikes, window)
#
#        for n in range(len(inmp_means)):
#            inmp_means[n] += frs[n]
#            inmp_fr_conv[gid % INmplus_gids][n] = frs[n]
#
#    # Calculate mean and standard deviation for all points
#    for n in range(len(inmp_means)):
#        inmp_means[n] /= numINmplus
#        inmp_stdevs[n] = np.std(inmp_fr_conv[:,n])
#
#    # Only plot one point each 1000 samples
#    plt_inmp_means = []
#    plt_inmp_stdevs = []
#    for n in np.arange(0,len(inmp_means), 1000):
#        plt_inmp_means.append(inmp_means[n])   
#        plt_inmp_stdevs.append(inmp_stdevs[n]) 
#
#    plt.plot(np.arange(0,len(inmp_means)/10,100), plt_inmp_means, 
#                 color='r', marker='^', mfc='r', mec='r', label='INm+')
#    plt.xlabel('Time (t) [ms]')
#
#    eus_means = np.zeros(sim.n_steps)
#    eus_stdevs = np.zeros(sim.n_steps)
#    eus_fr_conv = np.zeros((numEUSaff,sim.n_steps))
#    
#    for gid in np.arange(EUS_gids, EUS_gids + numEUSaff):
#        spikes = np.zeros(sim.n_steps, dtype=np.int)
#        if len(spike_trains.get_times(gid)) > 0:
#            spikes[(spike_trains.get_times(gid)/sim.dt).astype(np.int)] = 1
#        window_size = 10000
#        window = np.zeros(window_size)
#        window[0:window_size] = 1
#
#        frs = np.convolve(spikes, window)
#
#        for n in range(len(eus_means)):
#            eus_means[n] += frs[n]
#            eus_fr_conv[gid % EUS_gids][n] = frs[n]
#
#    # Calculate mean and standard deviation for all points
#    for n in range(len(eus_means)):
#        eus_means[n] /= numEUSaff
#        eus_stdevs[n] = np.std(eus_fr_conv[:,n])
#
#    # Only plot one point each 1000 samples
#    plt_eus_means = []
#    plt_eus_stdevs = []
#    for n in np.arange(0,len(eus_means), 1000):
#        plt_eus_means.append(eus_means[n])   
#        plt_eus_stdevs.append(eus_stdevs[n]) 
#
#    plt.plot(np.arange(0,len(eus_means)/10,100), plt_eus_means, 
#                 color='m', marker='^', mfc='m', mec='m', label='EUS Afferent')
#    plt.xlabel('Time (t) [ms]')
#    plt.ylabel('Neuron Firing Rate (FR) [Hz]')
#    plt.legend()
#
#    plt.show()

    bionet.nrn.quit_execution()

if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        run(sys.argv[-1])
    else:
        run('config.json')
