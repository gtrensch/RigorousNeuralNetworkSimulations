#
#  singleNeuronDynamics.py
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
# =   Plot the individual neuron dynamics, i.e., v(t), for a regular spiking (RS type) and a fast spiking
# =   (FS type) Izhikevich neuron, that is stimulated with an external current of 5.0pA.
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

import spynnaker8 as sim
from pyNN.utility.plotting import Figure, Panel
import pylab

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# =   Simulation time and resolution
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
SIM_TIME = 1000    # [ms]
sim.setup( timestep = 1.0 )

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# =   Neuron model and Izhikevich neuron parameters
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
NEURON_MODEL = sim.Izhikevich

# Regular spiking type Izhikevich neuron
NEURON_PARAMS_RS = { 'a':          0.02
                   , 'b':          0.2
                   , 'c':        -65.0
                   , 'd':          8.0
                   , 'v_init':   -75.0
                   , 'u_init':     0.0
                   , 'tau_syn_E':  0.01
                   , 'tau_syn_I':  0.01
                   , 'i_offset':   5.0
                   }

# Fast spiking type Izhikevich neuron
NEURON_PARAMS_FS = { 'a':          0.1
                   , 'b':          0.2
                   , 'c':        -65.0
                   , 'd':          2.0
                   , 'v_init':   -75.0
                   , 'u_init':     0.0
                   , 'tau_syn_E':  0.01
                   , 'tau_syn_I':  0.01
                   , 'i_offset':   5.0
                   }

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# =   Func: record membrane voltages to file in a printable ASCII format
# =
# =   t [ms]   v
# =   000     -75.0 mV
# =   001     -77.9880371094 mV
# =   002     -78.662109375 mV
# =   003     -78.6662902832 mV
# =   004     -78.4951477051 mV
# =   ...
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def writeVtoDisk( population, fileName ):
    outFile = open(fileName, 'w')
    x = population.get_data('v')
    v = x.segments[0].filter(name='v')[0]

    for i in range(len(v)):
        time = str('{:03d}'.format(i))
        voltage = str(v[i][0])
        line = ' '.join([time, voltage, "\n"])
        outFile.writelines(line)

    outFile.close()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# =   Create an RS-type and FS-type Izhikevich neuron
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
neuron_RS = sim.Population( 1, NEURON_MODEL(**NEURON_PARAMS_RS) )
neuron_FS = sim.Population( 1, NEURON_MODEL(**NEURON_PARAMS_FS) )

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# =   Set up recording of the membrane voltages for both neurons
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
neuron_RS.record('v')
neuron_FS.record('v')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# =   Run simulation
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
sim.run( SIM_TIME )

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# =   Retrieve and plot data, and write data to disk
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
v_neuron_RS = neuron_RS.get_data('v').segments[0].filter(name='v')[0]
v_neuron_FS = neuron_FS.get_data('v').segments[0].filter(name='v')[0]

writeVtoDisk( neuron_RS, "v(t)_RS_type.dat")
writeVtoDisk( neuron_FS, "v(t)_FS_type.dat")

Figure( Panel(v_neuron_RS, xticks = True, yticks = True, markersize = 2.0, xlabel = "regular spiking")
      , Panel(v_neuron_FS, xticks = True, yticks = True, markersize = 2.0, xlabel = "fast spiking")
      , annotations="Simulated with {}".format(sim.name())
      )

sim.end()
pylab.show()
