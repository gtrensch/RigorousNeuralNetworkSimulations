#
#  Copyright (C) 2021, G. Trensch, Forschungszentrum Jülich, JSC, Simulation & Data Laboratory Neuroscience
#
#  The refactored Izhikevich polychronization model application is free software:
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

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
#  PyNEST script of the two-population Izhikevich network model described in [1].
#
#  This implementation serves three purposes:
#  1) According to the findings from [2], [3] and [4] that the polychronization described in [1] is a consequence of
#     numerical inaccuracy and thus a simulation artifact, this implementation is aiming for numerical accuracy.
#  2) It is used in a validation task, i.e., the reproduction of network states exported from a
#     reference C implementation: https://github.com/gtrensch/RigorousNeuralNetworkSimulations/tree/master/C_model
#  3) It is used in a benchmarking task: Systematically varying an external input current from
#     i_ext = −3.0 pA to 100.0 pA allows to run the network through a wide range of activity,
#     from quiescence up to an average firing rate of 300 spks/s and thus simulating workload.
#
#
#  Prerequisites:
#  - NEST 2.20.1
#    The PyNEST script was developed using NEST version 2.20.1.
#    NEST was built using the following CMake options:
#    cmake -DCMAKE_INSTALL_PREFIX=$PWD/INSTALL \
#          -Dwith-python=3 \
#          -Dwith-optimize=ON \
#          -Dwith-warning=ON \
#          -Dwith-openmp=ON \
#          -Dwith-mpi=OFF \
#          -Dwith-gsl=ON \
#          -Dwith-ltdl=ON \
#          -Dwith-boost=ON \
#          -Dwith-debug=OFF \
#          -Dwith-version-suffix=IhkBuild \
#           ../../nest-simulator
#  - To achieve numerical accuracy while allowing to communication spike events in an one milliseconds interval,
#    the Izhikevich model implementation was adapted to progress neuron model dynamics in h = 0.1ms steps.
#    The source code modification applied to https://github.com/nest/nest-simulator/blob/v2.20.1/models/izhikevich.cpp
#    is shown below.
#
#    void nest::izhikevich::update( Time const& origin, const long from, const long to ) {
#      assert( to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
#      assert( from < to );
#
#      const double h = Time::get_resolution().get_ms();
#
#      bool thresholdDetected = false;
#
#      for ( long lag = from; lag < to; ++lag ) {
#
#        double I_syn = B_.spikes_.get_value( lag );
#
#        // h = 0.1ms, set nest.SetKernelStatus({"resolution": 1.0})
#        for(int i = 0; i < 10; ++i) {
#          double prev_V = S_.v_;
#          double prev_U = S_.u_;
#          S_.v_ += 0.1 * h * ( 0.04 * prev_V * prev_V + 5.0 * prev_V + 140.0 - prev_U + S_.I_ + P_.I_e_ + I_syn );
#          S_.u_ += 0.1 * h * P_.a_ * ( P_.b_ * prev_V - prev_U );
#
#          // lower bound of membrane potential
#          S_.v_ = ( S_.v_ < P_.V_min_ ? P_.V_min_ : S_.v_ );
#
#          // threshold detection inside ODE solver
#          if ( S_.v_ >= P_.V_th_ ) {
#            S_.v_ = P_.c_;
#            S_.u_ = S_.u_ + P_.d_;
#
#            thresholdDetected = true;
#          }
#        }
#
#        if ( thresholdDetected ) {
#          // compute spike time
#          set_spiketime( Time::step( origin.get_steps() + lag + 1 ) );
#
#          SpikeEvent se;
#          kernel().event_delivery_manager.send( *this, se, lag );
#        }
#
#        // set new input current
#        S_.I_ = B_.currents_.get_value( lag );
#
#        // voltage logging
#        B_.logger_.record_data( origin.get_steps() + lag );
#      }
#    }
#
# References:
# [1] Izhikevich, E. M. (2006). Polychronization: Computation with spikes. Neural Computation, 18:245–282.
#
# [2] Pauli, R., Weidel, P., Kunkel, S., and Morrison, A. (2018). Reproducing polychronization:
#     A guide to maximizing the reproducibility of spiking network models. Frontiers in Neuroinformatics, 12.
#
# [3] Trensch, G., Gutzen, R., Blundell, I., Denker, M., and Morrison, A. (2018). Rigorous neural network simulations:
#     A model substantiation methodology for increasing the correctness of simulation results in the absence of
#     experimental validation data. Frontiers in Neuroinformatics, 12:81.
#
# [4] Gutzen, R., von Papen, M., Trensch, G., Quaglio, P., Grün, S., and Denker, M. (2018). Reproducible neural network
#     simulations: Statistical methods for model validation on the level of network activity data.
#     Frontiers in Neuroinformatics, 12:90.
#
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

import nest
import nest.raster_plot
import matplotlib.pyplot as plt
import numpy as np
from networkConnectivity import *           # import network state
import time

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
#  TWO-POPULATION IZHIKEVICH NETWORK CLASS
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
class IzhkNetwork:
    def __init__(self, i_ext):
        # ---------------------------------------------------------
        # set population sizes
        self.numNrn_inh = 200
        self.numNrn_exc = 800
        self.numNrn_all = self.numNrn_exc + self.numNrn_inh
        self.indegree = 100

        # ---------------------------------------------------------
        # set regular spiking Izhikevich neuron parameters
        self.nrnParams_exc = {"a": 0.02,
                              "b": 0.2,
                              "c": -65.0,
                              "d": 8.0,
                              "U_m": 0.2 * -65.0,
                              "V_m": -65.0,
                              "I_e": i_ext}

        # ---------------------------------------------------------
        # set fast spiking Izhokevich neuron parameters
        self.nrnParams_inh = {"a": 0.1,
                              "b": 0.2,
                              "c": -65.0,
                              "d": 2.0,
                              "U_m": 0.2 * -65.0,
                              "V_m": -65.0,
                              "I_e": 0.0}

    def Build(self, sim_time):
        # ---------------------------------------------------------
        # CREATE NETWORK IN NEST,
        # RETURN SPIKE RECORDER
        # ---------------------------------------------------------
        # create populations
        self.nrnPop_exc = nest.Create("izhikevich", n = self.numNrn_exc, params = self.nrnParams_exc)
        self.nrnPop_inh = nest.Create("izhikevich", n = self.numNrn_inh, params = self.nrnParams_inh)
        self.nrnPop_all = self.nrnPop_exc + self.nrnPop_inh

        # ---------------------------------------------------------
        # create spike generators
        #
        # spike_gen 1    -> w -> neuron 1
        # spike_gen 2    -> w -> neuron 2
        # ..
        # spike_gen 1000 -> w -> neuron 1000
        #
        # each spike generator will be "loaded" with a sequence of spike times
        self.spikeGen = nest.Create("spike_generator", len(self.nrnPop_all))     # create as many spike generators as there are neurons in the network

        np.random.seed(1)                                                        # make random numbers replicable
        list_of_nrnIds = np.random.choice(self.numNrn_all, sim_time)             # as many neuronIds as milliseconds simulation; every millisecoand a new neuron is stimulated
        list_of_spike_times = np.array(np.linspace(0, sim_time, sim_time + 1))   # 0 is removed in next line; + 1 to get 1 ... sim_time
        list_of_spike_times = list_of_spike_times[list_of_spike_times > 0]       # remove zero time as it is not allowed for the spike generator
        # ---------------------------------------------------------
        # sort neuronIds with their spike times and load spike generators
        for i in np.unique(list_of_nrnIds):
            idx = list_of_nrnIds == i
            times = list_of_spike_times[idx]
            nest.SetStatus([self.spikeGen[i]], {'spike_times': times})           # NEST 2 expects a list for the generator

        # ---------------------------------------------------------
        # create recording device
        self.spikeRec = nest.Create("spike_detector")

        # ---------------------------------------------------------
        # connect populations
        # connections are imported from networkConnectivity.py; pre-trained network with static synapses
        nest.Connect(con_exc_to_all["source"], con_exc_to_all["target"], "one_to_one", {"weight": con_exc_to_all["weight"], "delay": con_exc_to_all["delay"]})
        nest.Connect(con_inh_to_exc["source"], con_inh_to_exc["target"], "one_to_one", {"weight": con_inh_to_exc["weight"], "delay": con_inh_to_exc["delay"]})

        # ---------------------------------------------------------
        # connect spike generators
        nest.Connect(self.spikeGen, self.nrnPop_all, 'one_to_one')
        nest.SetStatus(nest.GetConnections(self.spikeGen), 'weight', 20.0)

        # ---------------------------------------------------------
        # connect spike recording device
        nest.Connect(self.nrnPop_all, self.spikeRec)

        return(self.spikeRec)


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
#  HELPER FUNCTIONS
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
def WriteSpikesToDisk(outFileName, spikes, num_spikes, numNrn_all):
    # -------------------------------------------------------------
    # write spikes to disk
    # -------------------------------------------------------------
    import os
    outFile = open( outFileName, 'w' )
    # -------------------------------------------------------------
    # header
    outFile.write("#  size = {}".format(numNrn_all) + "\n"
                  "#  first_index = 1 \n"
                  "#  first_id = 1 \n"
                  "#  n = 1 \n"
                  "#  variable = spikes \n"
                  "#  last_id = {}".format(numNrn_all) + "\n"
                  "#  last_index = {}".format(numNrn_all) + "\n"
                  "#  dt = 0.1 \n"
                  "#  label = NEST \n")
    # -------------------------------------------------------------
    # write spikes to disk
    neurons = spikes['senders']
    times = spikes['times']
    for i in range(num_spikes):
        t = str('{:6.3f}'.format(times[i]))
        n = str('{:06d}'.format(neurons[i]))
        line = ' '.join([t, n, "\n"])
        outFile.writelines( line )

    outFile.close()

def WriteBenchmakResultsToDisk(outFileName, i_ext, spikesPerStep, accelFactor):
    # -------------------------------------------------------------
    # write benchmark results to disk
    # -------------------------------------------------------------
    import os
    outFile = open( outFileName, 'w' )
    outFile.write(" i_ext   spikes/dt   Accel. Factor \n")
    outFile.write("---------------------------------------- \n")
    for i in range(len(benchmark_spikesPerStep)):
        outFile.write(f"{list_of_i_ext[i]: >6.2f}    {benchmark_spikesPerStep[i]: >5.2f}        {benchmark_accelFactor[i]: >5.2f} \n")

    outFile.write("---------------------------------------- \n")
    outFile.close()

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
#  MAIN ENTRY
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
if __name__ == "__main__":

    # -------------------------------------------------------------
    # PARAMETERS:
    # select the task to perform 'VALIDATION' or 'BENCHMARK'
    RUN_MODE = "BENCHMARK"
    THREADS = 1
    OUTFILE_SPIKE_DATA = "./results/recordedSpikes.txt"
    OUTFILE_BENCHMARK_DATA = "./results/benchmark_1_thread.txt"
    # -------------------------------------------------------------

    # -------------------------------------------------------------
    # simulated time
    if RUN_MODE == "VALIDATION":
        T_SIM = 1000 * 60 * 30        # 30 minutes
    else:        # "BENCHMARK"
        T_SIM = 1000 * 60 * 5         # 5 minutes
    # -------------------------------------------------------------
    # to store benchmark statistics
    benchmark_spikesPerStep = []
    benchmark_accelFactor = []

    # -------------------------------------------------------------
    # set external current
    if RUN_MODE == "VALIDATION":
        # no external current
        # validation task: reproduce network state exported from C simulation
        list_of_i_ext = [0.0]
    else:        # "BENCHMARK"
        # increasing external current
        # benchmarking  task: simulate different workloads
        list_of_i_ext = [-3.0, -2.0, -1.0, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 1.5, 2.0, \
                          3.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, \
                          55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0]

    # -------------------------------------------------------------
    # run a simulation for each external current value in the list
    for i_ext in list_of_i_ext:
        # ---------------------------------------------------------
        # ensure clean state of NEST
        nest.ResetKernel()
        nest.SetKernelStatus({"local_num_threads": THREADS})
        # Although spikes are communicated in an one milliseconds interval, the resolution of state update is h = 0.1ms.
        # This is achieved by an adaption of the Izhikevich model implementation in NEST 2.20.1 (see the comment above).
        nest.SetKernelStatus({"resolution": 1.0})

        # ---------------------------------------------------------
        # instantiate the two-population Izhikevich network
        izhk_network = IzhkNetwork(i_ext)
        spikeRec = izhk_network.Build(T_SIM)

        # ---------------------------------------------------------
        # simulate
        startTime = time.time()
        nest.Simulate(T_SIM)
        endTime = time.time()

        # ---------------------------------------------------------
        # calculate and print statistics
        num_spikes = nest.GetStatus(spikeRec, "n_events")[0]
        spikes = nest.GetStatus(spikeRec, "events")[0]
        firingNeurons = spikes['senders']
        firingTimes = spikes['times']
        accelFactor = T_SIM / (endTime - startTime) / 1000
        average_firing_rate = (num_spikes / (T_SIM / 1000)) / izhk_network.numNrn_all
        spikesPerDt = average_firing_rate / 10

        print("")
        print("Two-population Izhikevich Network")
        print("----------------------------------------")
        print(f"Simulated time    : {T_SIM / 1000} sec.")
        print(f"Number of neurons : {izhk_network.numNrn_all}")
        print(f"External current  : {i_ext} pA")
        print(f"Number of spikes  : {num_spikes}")
        print(f"Avr. firing rate  : {average_firing_rate:.2f} spikes/sec.")
        print(f"Spikes per dt     : {spikesPerDt:.4f}")
        print(f"Duration          : {endTime - startTime:.2f} sec.")
        print(f"Accel. factor     : {accelFactor:.2f}")
        print("----------------------------------------")

        # ---------------------------------------------------------
        # store benchmark statistics
        benchmark_spikesPerStep.append(spikesPerDt)
        benchmark_accelFactor.append(accelFactor)

    if RUN_MODE == "VALIDATION":
        WriteSpikesToDisk(OUTFILE_SPIKE_DATA, spikes, num_spikes, izhk_network.numNrn_all)
        print("... spike data written to: ", OUTFILE_SPIKE_DATA)

    if RUN_MODE == "BENCHMARK":
        WriteBenchmakResultsToDisk(OUTFILE_BENCHMARK_DATA, list_of_i_ext, benchmark_spikesPerStep, benchmark_accelFactor)
        # ---------------------------------------------------------
        # print benchmark statistics
        print(" i_ext   Spikes/dt   Accel. factor")
        print("----------------------------------------")
        for i in range(len(benchmark_spikesPerStep)):
            print(f"{list_of_i_ext[i]: >6.2f}    {benchmark_spikesPerStep[i]: >5.2f}        {benchmark_accelFactor[i]: >5.2f}")

        print("----------------------------------------")

