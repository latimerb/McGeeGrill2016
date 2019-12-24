from bmtk.simulator import bionet
#from bmtk.simulator.bionet.default_setters.cell_models import loadHOC

# Load hoc cell templates
#bionet.pyfunction_cache.add_cell_model(loadHOC, directive='hoc', model_type='biophysical')

# Import the synaptic depression/facilitation model
import synapses
synapses.load()


conf = bionet.Config.from_json('simulation_configSIMPLE.json')
conf.build_env()
net = bionet.BioNetwork.from_config(conf)
sim = bionet.BioSimulator.from_config(conf, network=net)
sim.run()
