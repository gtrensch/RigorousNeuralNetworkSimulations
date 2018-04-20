#  polychronizationNetwork.py
#
#  This file is part of the SpiNNaker Izhikevich polychronization model implementation.
#
#  Copyright (C) 2018, Author: G. Trensch
#
#  The SpiNNaker Izhikevich polychronization model implementation is free software:
#  you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  It is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this application. If not, see <http://www.gnu.org/licenses/>.
#

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# =   SpiNNaker Izhikevich polychronization network implementation w/o STDP.
# =   Simulation of selected network states created with the C model implementation.
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

import spynnaker8 as sim
from pyNN.utility.plotting import Figure, Panel
import pylab
import os
import matplotlib.pyplot as plt
import numpy as np
import time
import seaborn


SPIKEPLOT = True

start_time = time.time()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# =   Import the network state that was generated with the C model
# =   Parameterize according to the simulated network and state
# =
# =                            source  target  weight    delay
# =   exc_exc_connections = [ (nnn,    nnn,    ww.wwww,  d ), () ... ]
# =   exc_inh_connections = [ (), () ... ]
# =   inh_exc_connections = [ (), () ... ]
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
EXTERNAL_CURRENT_FILE_NAME = 'randomNetworkInput.dat'         # random input


from pythonConWeightDelay_after1h import *                    # 1st state: connections, weights, and delays
FIRINGS_FILE_NAME  = 'firings_after1h.dat'                    #            network acivity data output file

# from pythonConWeightDelay_after2h import *                  # 2nd state
# FIRINGS_FILE_NAME  = 'firings_after2h.dat'                  #

# from pythonConWeightDelay_after3h import *                  # 3rd state
# FIRINGS_FILE_NAME  = 'firings_after3h.dat'                  #

# from pythonConWeightDelay_after4h import *                  # 4th state
# FIRINGS_FILE_NAME  = 'firings_after4h.dat'                  #

# from pythonConWeightDelay_after5h import *                  # 5th state
# FIRINGS_FILE_NAME  = 'firings_after5h.dat'                  #


# population sizes have to correspond to the imported network state
NUM_EXCITATORY_NEURONS  = 800
NUM_INHIBITORY_NEURONS  = 200
NUM_NEURONS             = NUM_EXCITATORY_NEURONS + NUM_INHIBITORY_NEURONS
SYNAPSES_PER_NEURON     = 100
CURRENT_INJECTION_VALUE = 20.0


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# =   Simulation parameters
# =   The simulation is carried out in 60 cycles with 1000 milliseconds simulation time each
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
SIM_TIME_TOTAL  = 1000 * 60                                  # [ms]
SIM_TIME_SINGLE = 1000                                       # [ms]
SIM_CYCLES      = SIM_TIME_TOTAL / SIM_TIME_SINGLE

SIM_TIME_STEP = 1.0                                          # [ms] (1ms = real time)
SIM_MIN_DELAY = 1.0
SIM_MAX_DELAY = 144.0

SIM_NEURONS_PER_CORE = 70

sim.setup( timestep = SIM_TIME_STEP, min_delay= SIM_MIN_DELAY, max_delay = SIM_MAX_DELAY )

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# =   Neuron model and Izhikevich neuron parameters
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
NEURON_MODEL = sim.Izhikevich


# Regular spiking type Izhikevich neuron
NEURON_PARAMS_EXCITATORY = { 'a':          0.02
                           , 'b':          0.2
                           , 'c':        -65.0
                           , 'd':          8.0
                           , 'v_init':   -65.0
                           , 'u_init':   -65.0 * 0.2         # according to the polychronization model
                           , 'i_offset':   0.0
                           }

# Fast spiking type Izhikevich neuron
NEURON_PARAMS_INHIBITORY = { 'a':          0.1
                           , 'b':          0.2
                           , 'c':        -65.0
                           , 'd':          2.0
                           , 'v_init':   -65.0
                           , 'u_init':   -65.0 * 0.2        # according to the polychronization model
                           , 'i_offset':   0.0
                           }

sim.set_number_of_neurons_per_core( NEURON_MODEL, SIM_NEURONS_PER_CORE )


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# =   Populations: create the two populations of the polychronization  model
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
pop_exc = sim.Population( NUM_EXCITATORY_NEURONS, NEURON_MODEL(**NEURON_PARAMS_EXCITATORY) )
pop_inh = sim.Population( NUM_INHIBITORY_NEURONS, NEURON_MODEL(**NEURON_PARAMS_INHIBITORY) )

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# =   Populations: create two spike source arrays for emulating external current and two empty lists,
# =                filled at begin of every cycle with the spike sequence, that is, the random input data
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
spikes_inp_exc = list()
for i in range( NUM_EXCITATORY_NEURONS ):
    spikes_inp_exc.append( [] )

pop_inp_exc = sim.Population( NUM_EXCITATORY_NEURONS, sim.SpikeSourceArray( spike_times = spikes_inp_exc ) )

spikes_inp_inh = list()
for i in range( NUM_INHIBITORY_NEURONS ):
    spikes_inp_inh.append( [] )

pop_inp_inh = sim.Population( NUM_INHIBITORY_NEURONS, sim.SpikeSourceArray( spike_times = spikes_inp_inh ) )

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# =   Projections: set up the poylchronization network
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
proj_exc_exc = sim.Projection( pop_exc, pop_exc
                             , sim.FromListConnector( exc_exc_connections )
                             , sim.StaticSynapse()
                             , receptor_type = "excitatory" )

proj_exc_inh = sim.Projection( pop_exc, pop_inh
                             , sim.FromListConnector( exc_inh_connections )
                             , sim.StaticSynapse()
                             , receptor_type = "excitatory" )

proj_inh_exc = sim.Projection( pop_inh, pop_exc
                             , sim.FromListConnector( inh_exc_connections )
                             , sim.StaticSynapse()
                             , receptor_type = "inhibitory" )

proj_inp_exc = sim.Projection( pop_inp_exc, pop_exc
                             , sim.OneToOneConnector()
                             , sim.StaticSynapse(weight = CURRENT_INJECTION_VALUE)
                             , receptor_type = "excitatory" )

proj_inp_inh = sim.Projection( pop_inp_inh, pop_inh
                             , sim.OneToOneConnector()
                             , sim.StaticSynapse(weight = CURRENT_INJECTION_VALUE)
                             , receptor_type = "excitatory" )

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# =   Set up recording
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
pop_exc.record('spikes')
pop_inh.record('spikes')
pop_inp_exc.record('spikes')
pop_inp_inh.record('spikes')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# =   Open files
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
inputCurrentFile = open( EXTERNAL_CURRENT_FILE_NAME, 'r' )

if os.path.isfile( FIRINGS_FILE_NAME ):
    os.remove( FIRINGS_FILE_NAME )

outputSpikesFile = open( FIRINGS_FILE_NAME, 'a' )


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# =   Set spike plot properties
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
if SPIKEPLOT:
    seaborn.set()
    seaborn.despine()
    seaborn.set_style("white")
    seaborn.axes_style("darkgrid")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# =   Run simulation
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
for cycle in range( SIM_CYCLES ):

    simStartTime = SIM_TIME_SINGLE * cycle

    spikes_inp_exc = list()
    for i in range(NUM_EXCITATORY_NEURONS):
        spikes_inp_exc.append([])

    spikes_inp_inh = list()
    for i in range(NUM_INHIBITORY_NEURONS):
        spikes_inp_inh.append([])

    # load random input date from file for this cycle
    for tickCount in range( SIM_TIME_SINGLE ):
        inputDataForThisTick = inputCurrentFile.readline().split()
        inputNeuronForThisTick = int(inputDataForThisTick[2])

        if( inputNeuronForThisTick < NUM_EXCITATORY_NEURONS ):
            spikes_inp_exc[inputNeuronForThisTick].append( tickCount + simStartTime )

        else:
            inputNeuronForThisTick = inputNeuronForThisTick - NUM_EXCITATORY_NEURONS
            spikes_inp_inh[inputNeuronForThisTick].append( tickCount + simStartTime )

    # set spike source array with random input data fo this cycle
    pop_inp_exc.set( spike_times = spikes_inp_exc )
    pop_inp_inh.set( spike_times = spikes_inp_inh )


    sim.run( SIM_TIME_SINGLE )

    # retrieve input spikes from populations, i.e., from spike source arrays
    spikesInpToExcPop = pop_inp_exc.get_data()
    spikesInpToInhPop = pop_inp_inh.get_data()

    # retrieve output spikes from populations and clear segment
    spikesExcPop = pop_exc.get_data('spikes', clear = True)
    spikesInhPop = pop_inh.get_data('spikes', clear = True)

    simEndTime = sim.get_current_time()
    print( simStartTime, simEndTime )

    if SPIKEPLOT:
        Figure( Panel(spikesInpToInhPop.segments[0].spiketrains
              , xticks = True, yticks = True, markersize = 2.0, xlim = (simStartTime, simEndTime))

              , Panel(spikesInpToExcPop.segments[0].spiketrains
              , xticks = True, yticks = True, markersize = 2.0, xlim = (simStartTime, simEndTime))

              , Panel(spikesInhPop.segments[0].spiketrains
              , xticks = True, yticks = True, markersize = 2.0, xlim = (simStartTime, simEndTime))

              , Panel(spikesExcPop.segments[0].spiketrains
              , xticks = True, yticks = True, markersize = 2.0, xlim = (simStartTime, simEndTime), xlabel = "Time [ms]")

              , title="Plochronization Network", annotations="Simulated with {}".format(sim.name())
              )

        plt.pause(0.1)

    # write firings to file
    firings = list()
    for i in range( SIM_TIME_SINGLE ):
        firings.append([])

    for neuron in range( NUM_EXCITATORY_NEURONS ):
        for spikeTimeIdx in range( len(spikesExcPop.segments[0].spiketrains[neuron]) ):
            firingTime = int(spikesExcPop.segments[0].spiketrains[neuron][spikeTimeIdx]) - simStartTime
            firings[firingTime].append( neuron )

    for neuron in range( NUM_INHIBITORY_NEURONS ):
        for spikeTimeIdx in range( len(spikesInhPop.segments[0].spiketrains[neuron]) ):
            firingTime = int(spikesInhPop.segments[0].spiketrains[neuron][spikeTimeIdx]) - simStartTime
            firings[firingTime].append( neuron + NUM_EXCITATORY_NEURONS )

    for i in range( SIM_TIME_SINGLE ):
        simCycle = str( '{:06d}'.format(cycle) )
        tickInCycle = str('{:06d}'.format(i))

        for n in firings[i]:
            neuronId = str('{:03d}'.format(n))
            line = ' '.join([simCycle, tickInCycle, neuronId, "\n"])
            outputSpikesFile.writelines( line )


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# =   End simulation, close files, show spike plot
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
sim.end()

inputCurrentFile.close()
outputSpikesFile.close()

print("Simulation time: %s seconds" % (time.time() - start_time))

if SPIKEPLOT:
    pylab.show()
