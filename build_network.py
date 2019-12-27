import numpy as np
from bmtk.builder.networks import NetworkBuilder
import math
import random

random.seed(42)
output_dir='network'

#######################################################################
##################### Create the cells ################################
#######################################################################
print("\nCreating Cells")

# Build the main network
net = NetworkBuilder('LUT')

# Specify number of cells in each population #
numPUDaff   = 1
numPelaff   = 1
numINmplus  = 1
numINmminus = 1
numIND      = 1
numFB       = 1
numSPN      = 1

# Create the nodes ----------------------------------------
net.add_nodes(N=numPelaff, pop_name='Pelaff',model_type='biophysical',model_template='hoc:LIF_adapt',morphology=None)
net.add_nodes(N=numPUDaff, pop_name='PUDaff',model_type='biophysical',model_template='hoc:LIF_adapt',morphology=None)
net.add_nodes(N=numINmplus, pop_name='INmplus',model_type='biophysical',model_template='hoc:LIF_adapt',morphology=None)
net.add_nodes(N=numINmminus, pop_name='INmminus',model_type='biophysical',model_template='hoc:LIF_adapt',morphology=None)
net.add_nodes(N=numIND, pop_name='IND',model_type='biophysical',model_template='hoc:LIF_adapt',morphology=None)
net.add_nodes(N=numFB, pop_name='FB',model_type='biophysical',model_template='hoc:LIF_adapt',morphology=None)
net.add_nodes(N=numSPN, pop_name='SPN',model_type='biophysical',model_template='hoc:LIF_adapt',morphology=None) 

##################################################################################
####################### Connect the cells ########################################
##################################################################################
print("\nConnecting Cells")


def conn_props(source,target,mu,sigma):
    """
    Simply add connection properties as normal
    https://github.com/AllenInstitute/bmtk/blob/develop/docs/tutorial/NetworkBuilder_Intro.ipynb

    Can also adjust weights before running the network.
    See https://github.com/AllenInstitute/bmtk/blob/develop/docs/tutorial/02_single_cell_syn.ipynb
    Section 5
    If done this way the function will need to be imported in the run script, consider refactoring?
    """

    #syn_weight = np.random.lognormal(mean=mu,sigma=sigma)
    syn_weight = mu

    return syn_weight,0,0.5


def one_to_one(source,target):
    tmp_syn = 1
    return tmp_syn

# Add connections -----------------------------------------

conn = net.add_edges(source=net.nodes(pop_name='PUDaff'), target=net.nodes(pop_name='INmminus'),
                   connection_rule=one_to_one,
                   delay=0.5,
                   syn_weight = 0.07,
                   target_sections=['somatic'],
                   distance_range=[0.0, 300.0],
                   dynamics_params='AMPA_ExcToExc.json',
                   model_template='Exp2Syn')
				
conn = net.add_edges(source=net.nodes(pop_name='PUDaff'), target=net.nodes(pop_name='INmplus'),
                   connection_rule=one_to_one,
                   delay=0.5,
                   syn_weight = 0.044,
                   target_sections=['somatic'],
                   distance_range=[0.0, 300.0],
                   dynamics_params='AMPA_ExcToExc.json',
                   model_template='Exp2Syn')

conn = net.add_edges(source=net.nodes(pop_name='PUDaff'), target=net.nodes(pop_name='IND'),
                   connection_rule=one_to_one,
                   delay=0.5,
                   syn_weight = 0.06,
                   target_sections=['somatic'],
                   distance_range=[0.0, 300.0],
                   dynamics_params='AMPA_ExcToExc.json',
                   model_template='Exp2Syn')


conn = net.add_edges(source=net.nodes(pop_name='Pelaff'), target=net.nodes(pop_name='IND'),
                   connection_rule=one_to_one,
                   delay=0.5,
                   syn_weight = 0.07,
                   target_sections=['somatic'],
                   distance_range=[0.0, 300.0],
                   dynamics_params='AMPA_ExcToExc.json',
                   model_template='Exp2Syn')

conn = net.add_edges(source=net.nodes(pop_name='INmminus'), target=net.nodes(pop_name='SPN'),
                   connection_rule=one_to_one,
                   delay=0.5,
                   syn_weight = 0.065,
                   target_sections=['somatic'],
                   distance_range=[0.0, 300.0],
                   dynamics_params='GABA_InhToExc.json',
                   model_template='Exp2Syn')

conn = net.add_edges(source=net.nodes(pop_name='INmplus'), target=net.nodes(pop_name='SPN'),
                   connection_rule=one_to_one,
                   delay=0.5,
                   syn_weight = 0.06,
                   target_sections=['somatic'],
                   distance_range=[0.0, 300.0],
                   dynamics_params='AMPA_ExcToExc.json',
                   model_template='Exp2Syn')

conn = net.add_edges(source=net.nodes(pop_name='IND'), target=net.nodes(pop_name='SPN'),
                   connection_rule=one_to_one,
                   delay=0.5,
                   syn_weight = 0.07,
                   target_sections=['somatic'],
                   distance_range=[0.0, 300.0],
                   dynamics_params='AMPA_ExcToExc.json',
                   model_template='Exp2Syn')

conn = net.add_edges(source=net.nodes(pop_name='FB'), target=net.nodes(pop_name='IND'),
                   connection_rule=one_to_one,
                   delay=0.5,
                   syn_weight = 0.09,
                   target_sections=['somatic'],
                   distance_range=[0.0, 300.0],
                   dynamics_params='GABA_InhToExc.json',
                   model_template='Exp2Syn')

conn = net.add_edges(source=net.nodes(pop_name='FB'), target=net.nodes(pop_name='IND'),
                   connection_rule=one_to_one,
                   delay=0.5,
                   syn_weight = 0.06,
                   target_sections=['somatic'],
                   distance_range=[0.0, 300.0],
                   dynamics_params='GABA_InhToExc.json',
                   model_template='Exp2Syn')

conn = net.add_edges(source=net.nodes(pop_name='SPN'), target=net.nodes(pop_name='FB'),
                   connection_rule=one_to_one,
                   delay=0.5,
                   syn_weight = 0.07,
                   target_sections=['somatic'],
                   distance_range=[0.0, 300.0],
                   dynamics_params='AMPA_ExcToExc.json',
                   model_template='Exp2Syn')

conn = net.add_edges(source=net.nodes(pop_name='PMC'), target=net.nodes(pop_name='IND'),
                   connection_rule=one_to_one,
                   delay=0.5,
                   syn_weight = 0.033,
                   target_sections=['somatic'],
                   distance_range=[0.0, 300.0],
                   dynamics_params='AMPA_ExcToExc.json',
                   model_template='Exp2Syn')


net.build()
net.save_nodes(output_dir=output_dir)
net.save_edges(output_dir=output_dir)

#inp = NetworkBuilder('INPUT') # Virtual cells delivering input to PUD
#inp.add_nodes(N=1, pop_name = 'Pel_aff_virt', model_type='virtual', potential='exc')
#inp.add_edges(source=inp.nodes(pop_name='Pel_aff_virt'), target=net.nodes(pop_name='SPN'),
#                   connection_rule=1,
#                   syn_weight=0.06,
#                   target_sections=['somatic'],
#				   delay=0.5,
#                   distance_range=[0.0, 300.0],
#                   dynamics_params='AMPA_ExcToExc.json',
#                   model_template='Exp2Syn')
#
#inp.build()
#inp.save_nodes(output_dir=output_dir)
#inp.save_edges(output_dir=output_dir)
####################################################################################
########################## Build and save network ##################################
####################################################################################

print("\nBuilding network and saving to directory \"" + output_dir + "\"")

from bmtk.utils.sim_setup import build_env_bionet

#build_env_bionet(base_dir='./',      
#                 network_dir='network',
#                 tstop=3000.0, dt=0.1,
#                 report_vars=['v'],     # Record membrane potential and calcium (default soma)
#                 spikes_inputs=[('INPUT',   # Name of population which spikes will be generated for
#                                'PUD_spikes.h5')],
#                 include_examples=True,    # Copies components files
#                 compile_mechanisms=True   # Will try to compile NEURON mechanisms
#                )


#from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator
#
#psg = PoissonSpikeGenerator(population='PUDaff')
#psg.add(node_ids=range(1),  # Have nodes to match mthalamus
#        firing_rate=4.0,    # 15 Hz, we can also pass in a nonhomoegenous function/array
#        times=(0.5, 1.0))    # Firing starts at 0 s up to 3 s
#psg.to_sonata('PUD_spikes.h5')


print("Done")
