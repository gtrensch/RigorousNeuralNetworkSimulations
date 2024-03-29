/*
 *  main.cpp
 *
 *  Copyright (C) 2018, G. Trensch, Forschungszentrum Jülich, JSC, Simulation & Data Laboratory Neuroscience
 *
 *  The refactored Izhikevich polychronization model application is free software:
 *  you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  It is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this application. If not, see <http://www.gnu.org/licenses/>.
 *
 */

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
// Refactored version of the two-population network described in [2]. This source code is based on the poly_spnet.cpp
// and spnet.cpp source code available online at: https://www.izhikevich.org/publications/spnet.htm.
//
// References:
//
// [1] Izhikevich, E. M. (2003). Simple model of spiking neurons. Trans. Neur. Netw. 14, 1569–1572.
//     doi:10.1109/TNN.2003.820440
// [2] Izhikevich, E. M. (2006). Polychronization: Computation with spikes. Neural Computation, 18:245–282.
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

#include "params.h"

#if( __RUN_REFACTORED_VERSION__ )

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <chrono>
#include <stdio.h>
#include <iostream>

#include "globals.h"
#include "utils.h"


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   C H E C K   P A R A M E T E R   S E T T I N G S
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#ifndef SIM_TIME
  #error SIM_TIME not defined!
#endif

#ifdef GENERATE_NETWORK_FROM_EXTERNAL_DATA
  #if( GENERATE_NETWORK_FROM_EXTERNAL_DATA )
    #ifndef INFILE_CONNECTION_MATRIX
      #error INFILE_CONNECTION_MATRIX not defined!
    #endif
    #ifndef INFILE_DELAY_MATRIX
      #error INFILE_DELAY_MATRIX not defined!
    #endif
    #ifndef INFILE_WEIGHT_MATRIX
      #error INFILE_WEIGHT_MATRIX not defined!
    #endif
  #else
    #ifndef OUTFILE_CONNECTIONS
      #error OUTFILE_CONNECTIONS not defined!
    #endif
    #ifndef OUTFILE_DELAYS
      #error OUTFILE_DELAYS not defined!
    #endif
    #ifndef OUTFILE_STIMULUS
      #error OUTFILE_STIMULUS not defined!
    #endif
  #endif
#else
  #error GENERATE_NETWORK_FROM_EXTERNAL_DATA not defined!
#endif

#ifdef USE_EXTERNAL_STIMULUS
  #if( USE_EXTERNAL_STIMULUS )
    #ifndef INFILE_STIMULUS
      #error INFILE_STIMULUS not defined!
    #endif
  #endif
#else
  #error USE_EXTERNAL_STIMULUS not defined!
#endif

#ifndef USE_STDP
  #error USE_STDP not defined!
#endif

#ifndef OUTFILE_FIRINGS
  #error OUTFILE_FIRINGS not defined!
#endif


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   F O R W A R D   D E C L A R A T I O N S
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void InitializeNetwork();
void InitializeSimulation();
void FinalizeSimulation();


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   M A I N   E N T R Y
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
int main() {
  unsigned long int iterCount  = 0;
  unsigned long int spikeCount = 0;
  printf( "\n\n%s\n\n", PARAM_INFO_STRING );

  printf( "[INFO]  Running simulation for %d seconds\n", SIM_TIME );

#if( ODE_SOLVER_REFINEMENT )
  for( int i = 0; i < NUM_TOTAL_NEURONS; ++i ) {
    spike[i] = 0;
  }
#endif

  InitializeSimulation();
  InitializeNetwork();

  std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();

  for( int simTimeSecond = 0; simTimeSecond < SIM_TIME; ++simTimeSecond ) {          // simulation loop [seconds]
    for( int t = 0; t < 1000; ++t ) {                                                // simulation loop [milliseconds]

      // clear I_ext[] and select a random neuron for input
      for( int i = 0; i < NUM_TOTAL_NEURONS; ++i ) {
        I_ext[i] = 0.0;
      }

      //  REFACTOR COMMENT: the original implementation does not allow networks sizes < 1000 neurons
      //                    for( int idx = 0; idx < NUM_TOTAL_NEURONS / 1000; ++idx ) { ... }
#if( USE_EXTERNAL_STIMULUS )
      int inputNeuron = GetNextExternalStimulusFromFile(INFILE_STIMULUS, simTimeSecond, t);
#else
      int inputNeuron = GET_RANDOM_INT( NUM_TOTAL_NEURONS );
#endif
      if( inputNeuron > 0 && inputNeuron != RC_NOID) {
        I_ext[inputNeuron] = I_EXT;
      }

#ifdef OUTFILE_STIMULUS
      RecordRandomStimulusToFile( OUTFILE_STIMULUS, simTimeSecond, t, inputNeuron );
#endif

#if( ODE_SOLVER_REFINEMENT )
      for( int n = 0; n < NUM_TOTAL_NEURONS; ++n ) {
        if( spike[n] > 0 ) {                                                         // Did neuron n fire?

#if( REPORT_ON_MULTIPLE_SPIKES )
          if( spike[n] > 1 ) printf( "[INFO]  Multiple spikes in simulation time step detected.\n" );
#endif
          spike[n] = 0;

          LTP[n][t + MAX_SYNAPSE_DELAY] = 0.1;
          LTD[n] = 0.12;

          for( int idx = 0; idx < numPreSynapticNeuronsOfTarget[n]; ++idx ) {        // LTP: pre-synaptic spike precedes post-synaptic spike
            int preSynNeuron = listOfPresynapticNeurons[n][idx];
            int preSynNeuron_correspondingDelay = listOfPresynapticDelays[n][idx];

            *listOfPointersToSynapticWeights_derivatives[n][idx] += LTP[preSynNeuron][t + MAX_SYNAPSE_DELAY - preSynNeuron_correspondingDelay - 1];
          }

          firings[numFirings][IDX_TIME] = t;
          firings[numFirings][IDX_NEURON] = n;
          numFirings++;
          if( numFirings == MAX_NUM_FIRINGS) {
            printf( "[WARNING]  Too many spikes at t = %d (ignoring all)\n", t );
            numFirings = 1;
          }
        }
      }
#else
      for( int n= 0; n < NUM_TOTAL_NEURONS; ++n ) {
        if( v[n] >= RS_FS_THR ) {                                                    // Did neuron n fire?

          v[n] = RS_FS_C;                                                            // threshold dynamics
          u[n] += d[n];

          LTP[n][t + MAX_SYNAPSE_DELAY] = 0.1;
          LTD[n] = 0.12;

          for ( int idx = 0; idx < numPreSynapticNeuronsOfTarget[n]; ++idx ) {       // LTP: pre-synaptic spike precedes post-synaptic spike
            int preSynNeuron = listOfPresynapticNeurons[n][idx];
            int preSynNeuron_correspondingDelay = listOfPresynapticDelays[n][idx];

            *listOfPointersToSynapticWeights_derivatives[n][idx] += LTP[preSynNeuron][t + MAX_SYNAPSE_DELAY - preSynNeuron_correspondingDelay - 1];
          }

          firings[numFirings][IDX_TIME]   = t;
          firings[numFirings][NEURON] = n;
          numFirings++;
          if( numFirings == MAX_NUM_FIRINGS ) {
            printf( "[INFO]  Two many spikes at t = %d (ignoring all).\n", t );
            numFirings = 1;
          }
        }
      }
#endif

      //  Go backwards through firings as far as MAX_SYNAPSE_DELAY lasts. Calculate I_ext for the next time step.
      //  Take the delays of the synapses into account. I_ext applies to the next simulation step.
      //
      //                          |        |                    |
      //     |                          |     |  |  |     |  |  |  <-- spikes
      //  ---+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--> simulation time
      //    -5 -4 -3 -2 -1  0  1  2  3  4  5  6  7  8  9 10 11 12
      //                                                        t
      //
      //                                         |              |
      //                                         +--------------+  MAX_SYNAPSE_DELAY (e.g., 5)

      long int idx = numFirings - 1;                                                      // array index starts with 0
      while( t - firings[idx][IDX_TIME] < MAX_SYNAPSE_DELAY ) {

        int neuronFired = firings[idx][IDX_NEURON];
        int timeSinceNeuronFired = t - firings[idx][IDX_TIME];
        int numEntriesOfDelay = numEntriesPerDelay[neuronFired][timeSinceNeuronFired];

        for( int i = 0; i < numEntriesOfDelay; ++i ) {
          int s = listOfSynsByNeuronAndDelay[neuronFired][timeSinceNeuronFired][i];
          int postSynNeuron = matrixOfPostSynapticNeurons[neuronFired][s];

          I_ext[postSynNeuron] += matrixOfSynapticWeights[neuronFired][s];

          if( neuronFired < NUM_EXCITATORY_NEURONS) {                                  // LTD: this spike is before postsynaptic spikes
            matrixOfSynapticWeights_derivatives[neuronFired][s] -= LTD[postSynNeuron];
          }
        }
        idx--;
      }

      // = = = = = = = = = = = = = = = = = = = = = =
      // = UPDATE NETWORK STATE
      // = = = = = = = = = = = = = = = = = = = = = =
      for( int n = 0; n < NUM_TOTAL_NEURONS; ++n ) {

#if( ODE_SOLVER_REFINEMENT )
        for( int intCycles = 0; intCycles < ODE_SOLVER_STEPS; ++intCycles ) {

  #if( ODE_SOLVER_NOT_SYMPLECTIC )
          double v_old = v[n];
          v[n] += (1.0 / ODE_SOLVER_STEPS) * ((0.04 * v[n] + 5) * v[n] + 140 - u[n] + I_ext[n] + i_offs[n]);
          u[n] += (1.0 / ODE_SOLVER_STEPS) * (a[n] * (RS_FS_B * v_old - u[n]));
          ++iterCount;
  #else
          v[n] += (1.0 / ODE_SOLVER_STEPS) * ((0.04 * v[n] + 5) * v[n] + 140 - u[n] + I_ext[n] + i_offs[n]);
          u[n] += (1.0 / ODE_SOLVER_STEPS) * (a[n] * (RS_FS_B * v[n] - u[n]));
          ++iterCount;
  #endif

  #if( LOG_MIN_MAX_V_U )
          // log min and max values of v(t) and u(t)
          if( v[n] > v_max ) v_max = v[n];
          if( v[n] < v_min ) v_min = v[n];
          if( u[n] > u_max ) u_max = u[n];
          if( u[n] < u_min ) u_min = u[n];
  #endif

          // threshold detection within loop for exact detection and precise spike timing
          if( v[n] >= 30.0 ) {
            v[n] = RS_FS_C;
            u[n] += d[n];
            spike[n]++;                                                              // remember the spike event, which is aligned to next grid point
            spikeCount++;
          }
        }
#else   // original Izhikevich implementation using a 1 ms integration time step
        v[n] += 0.5 * ((0.04 * v[n] + 5) * v[n] + 140 - u[n] + I_ext[n] + i_offs[n]);    // for numerical stability
        v[n] += 0.5 * ((0.04 * v[n] + 5) * v[n] + 140 - u[n] + I_ext[n] + i_offs[n]);    // time step is 0.5 ms, see [1]
        u[n] += a[n] * (RS_FS_B * v[n] - u[n]);

        #if( LOG_MIN_MAX_V_U )                                                       // log min and max values of v(t) and u(t)
          if(v[n] > v_max)  v_max = v[n];
          if(v[n] < v_min)  v_min = v[n];
          if(u[n] > u_max)  u_max = u[n];
          if(u[n] < u_min)  u_min = u[n];
          #endif
#endif

        LTP[n][t + MAX_SYNAPSE_DELAY + 1] = 0.95 * LTP[n][t + MAX_SYNAPSE_DELAY];
        LTD[n] *= 0.95;
      }
    }                                            // end of simulation loop [milliseconds]

    // printf("aver. firing rate = %f \n", (double)numFirings / (double)NUM_TOTAL_NEURONS);

    // every minute: report on the average firing rate of the excitatory and inhibitory population
    if( simTimeSecond > 0 && (simTimeSecond + 1) % 60 == 0 ) {
      unsigned long int spikesTotalExcitatory = 0;
      unsigned long int spikesTotalInhibitory = 0;
      for( unsigned long
      int idx = 0; idx < numFirings; ++idx ) {
        if( firings[idx][IDX_NEURON] < NUM_EXCITATORY_NEURONS) {
          spikesTotalExcitatory++;
        }
        else {
          spikesTotalInhibitory++;
        }
      }
      printf( "[INFO]  Simulated time: %d seconds\n", simTimeSecond + 1 );
      printf( "[INFO]  ... Firing rate (excitatory) = %f (total = %f, neurons = %f)\n"
            , (double) spikesTotalExcitatory / (double) NUM_EXCITATORY_NEURONS, (double) spikesTotalExcitatory, (double) NUM_EXCITATORY_NEURONS );
      printf( "[INFO]  ... Firing rate (inhibitory) = %f (total = %f, neurons = %f)\n"
            , (double) spikesTotalInhibitory / (double) NUM_INHIBITORY_NEURONS, (double) spikesTotalInhibitory, (double) NUM_INHIBITORY_NEURONS );

#if(LOG_MIN_MAX_V_U)
      printf( "[INFO]  ... v_max = %f\n", v_max );
      printf( "[INFO]  ... v_min = %f\n", v_min );
      printf( "[INFO]  ... u_max = %f\n", u_max );
      printf( "[INFO]  ... u_min = %f\n", u_min );
#endif
    }

    // = = = = = = = = = = = = = = = = = = = = = =
    // = SAVE FIVE SELECTED NETWORK STATES
    // = = = = = = = = = = = = = = = = = = = = = =
#ifdef SELECTED_STATE_1_AFTER_N_SECONDS
  #ifndef OUTFILE_STATE_1_WEIGHT_MATRIX
    #error OUTFILE_STATE_1_WEIGHT_MATRIX not defined!
  #endif
    if( simTimeSecond == SELECTED_STATE_1_AFTER_N_SECONDS - 1 ) {
      ExportWeightMatrixToFile( OUTFILE_STATE_1_WEIGHT_MATRIX );
    }
#endif

#ifdef SELECTED_STATE_2_AFTER_N_SECONDS
  #ifndef OUTFILE_STATE_2_WEIGHT_MATRIX
    #error OUTFILE_STATE_2_WEIGHT_MATRIX not defined!
  #endif
    if( simTimeSecond == SELECTED_STATE_2_AFTER_N_SECONDS - 1 ) {
      ExportWeightMatrixToFile( OUTFILE_STATE_2_WEIGHT_MATRIX );
    }
#endif

#ifdef SELECTED_STATE_3_AFTER_N_SECONDS
  #ifndef OUTFILE_STATE_3_WEIGHT_MATRIX
    #error OUTFILE_STATE_3_WEIGHT_MATRIX not defined!
  #endif
    if( simTimeSecond == SELECTED_STATE_3_AFTER_N_SECONDS - 1 ) {
      ExportWeightMatrixToFile( OUTFILE_STATE_3_WEIGHT_MATRIX );
    }
#endif

#ifdef SELECTED_STATE_4_AFTER_N_SECONDS
  #ifndef OUTFILE_STATE_4_WEIGHT_MATRIX
    #error OUTFILE_STATE_4_WEIGHT_MATRIX not defined!
  #endif
    if( simTimeSecond == SELECTED_STATE_4_AFTER_N_SECONDS - 1 ) {
      ExportWeightMatrixToFile( OUTFILE_STATE_4_WEIGHT_MATRIX );
    }
#endif

#ifdef SELECTED_STATE_5_AFTER_N_SECONDS
  #ifndef OUTFILE_STATE_5_WEIGHT_MATRIX
    #error OUTFILE_STATE_5_WEIGHT_MATRIX not defined!
  #endif
    if( simTimeSecond == SELECTED_STATE_5_AFTER_N_SECONDS - 1 ) {
      ExportWeightMatrixToFile( OUTFILE_STATE_5_WEIGHT_MATRIX );
    }
#endif

#if( RECORD_SPIKE_DATA )
    RecordNetworkActivityToFile( OUTFILE_FIRINGS, simTimeSecond, numFirings );
#endif

    // = = = = = = = = = = = = = = = = = = = = = =
    // = PREPARE NEXT SIMULATION TIME STEP
    // = = = = = = = = = = = = = = = = = = = = = =
    for( int n = 0; n < NUM_TOTAL_NEURONS; ++n ) {
      for( int i = 0; i < MAX_SYNAPSE_DELAY + 1; ++i ) {
        LTP[n][i] = LTP[n][1000 + i];
      }
    }

    long int k = numFirings - 1;
    while( 1000 - firings[k][IDX_TIME] < MAX_SYNAPSE_DELAY) {
      k--;
    }
    for( long int i = 1; i < numFirings - k; ++i ) {
      firings[i][IDX_TIME] = firings[k + i][IDX_TIME] - 1000;
      firings[i][IDX_NEURON] = firings[k + i][IDX_NEURON];
    }
    numFirings = numFirings - k;

#if( USE_STDP )
    // = = = = = = = = = = = = = = = = = = = = = =
    // = UPDATE SYNAPTIC WEIGHTS
    // = = = = = = = = = = = = = = = = = = = = = =
    // only excitatory connections are plastic
    for( int n = 0; n < NUM_EXCITATORY_NEURONS; ++n ) {
      for( int s = 0; s < NUM_SYNAPSES_PER_NEURON; ++s ) {

        // REFACTOR COMMENT:
        // The following two lines have a different order in spnet.cpp.
        matrixOfSynapticWeights_derivatives[n][s] *= 0.9;
        matrixOfSynapticWeights[n][s] += 0.01 + matrixOfSynapticWeights_derivatives[n][s];

        if( matrixOfSynapticWeights[n][s] > MAX_SYNAPTIC_STRENGTH) {
          matrixOfSynapticWeights[n][s] = MAX_SYNAPTIC_STRENGTH;
        }

        if( matrixOfSynapticWeights[n][s] < 0 ) {
          matrixOfSynapticWeights[n][s] = 0.0;
        }
      }
    }
#endif
  }                     // end of simulation loop [seconds]

  std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();

  FinalizeSimulation();

  std::cout << "Simulated time         : " << SIM_TIME << "seconds" << std::endl;
  std::cout << "Simulation duration    : " << (std::chrono::duration_cast<std::chrono::microseconds>( endTime - startTime ).count()) / 1000000.0 << " seconds" << std::endl;
  std::cout << "Acceleration factor    : " << SIM_TIME / ((std::chrono::duration_cast<std::chrono::microseconds>( endTime - startTime ).count()) / 1000000.0)  << std::endl;
  std::cout << "Total ODE solver steps : " << iterCount << std::endl;
  std::cout << "Total number of spikes : " << spikeCount << std::endl;
  std::cout << "Average firing rate    : " << static_cast<float>(spikeCount)/static_cast<float>(SIM_TIME)/static_cast<float>(NUM_TOTAL_NEURONS) << std::endl;
  float spikesPerMSec = (static_cast<float>(spikeCount)/static_cast<float>(SIM_TIME))/1000;
  std::cout << "Spikes per time step   : " << spikesPerMSec/static_cast<float>(ODE_SOLVER_STEPS) << std::endl;

  printf( "\n\nSimulation terminated normally! \n\n" );

  return (0);
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = Initialize the polychronization network
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void InitializeNetwork() {

  // initialize izhikevich neuron parameters
  for( int n = 0; n < NUM_TOTAL_NEURONS; ++n ) {
    if( n < NUM_EXCITATORY_NEURONS ) {                                               // excitatory neurons, RS type
      a[n] = RS_A;
      d[n] = RS_D;
      v[n] = RS_V_INIT;
      u[n] = RS_U_INIT;
      i_offs[n] = I_OFFS_RS;
    }
    else {                                                                           // inhibitory neurons, FS type
      a[n] = FS_A;
      d[n] = FS_D;
      v[n] = FS_V_INIT;
      u[n] = FS_U_INIT;
      i_offs[n] = I_OFFS_FS;
    }
  }


  // = = = = = = = = = = = = = = = = = = = = = =
  // = CONNECTION MATRIX
  // = = = = = = = = = = = = = = = = = = = = = =
#if( GENERATE_NETWORK_FROM_EXTERNAL_DATA )
  ImportConnectionMatrixFromFile( INFILE_CONNECTION_MATRIX );
#else
  // fill matrixOfPostSynapticNeurons[n][n] with random target neurons
  // avoid autapses and multapses

  int numOfTargets[NUM_TOTAL_NEURONS] = {'\0'};
  int numRedrawnConnections = 0;

#if( BALANCED_OUTGOING_CONNS )
  for( int n = NUM_TOTAL_NEURONS - 1; n >= 0; --n ) {                                // draw inhibitory connections first
#else
  for( int n = 0; n < NUM_TOTAL_NEURONS; ++n ) {
#endif
    for( int s = 0; s < NUM_SYNAPSES_PER_NEURON; ++s ) {

      bool duplicateTarget = false;
      bool selfAssigned = false;
      bool limitTargets = false;
      int randomTargetNeuron = 0;

      // repeat and redraw connection in case of an autapse, multapse or
      // too much outgoing connections detected (if BALANCED_OUTGOING_CONNS is set)
      do {
        duplicateTarget = false;
        selfAssigned = false;
        limitTargets = false;

        if( n < NUM_EXCITATORY_NEURONS ) {                                           // excitatory neurons can connect to all neurons
          randomTargetNeuron = GET_RANDOM_INT( NUM_TOTAL_NEURONS );
        }
        else {                                                                       // inhibitory neurons can connect to excitatory neurons only
          randomTargetNeuron = GET_RANDOM_INT( NUM_EXCITATORY_NEURONS );
        }
        if( randomTargetNeuron == n ) {
          selfAssigned = true;
        }
        for( int i = 0; i < s; ++i ) {
          if( matrixOfPostSynapticNeurons[n][i] == randomTargetNeuron ) {
            duplicateTarget = true;
          }
        }

        if( BALANCED_OUTGOING_CONNS && (numOfTargets[randomTargetNeuron] > LIMIT_OUTGOING_CONNS)) {
          limitTargets = true;
          numRedrawnConnections++;

          // loop detection
          // if LIMIT_OUTGOING_CONNS is set too small, this my end in an infinite loop as there might be no free entry in any
          // neuron's target list
          if(numRedrawnConnections > 20000) {
            // if LIMIT_OUTGOING_CONNS is set too small, we may never
            printf( "[ERROR]  Number of redrawn connections exceeds 20000. Possible loop at %s %d\n", __FILE__, __LINE__);
            fflush(stdout);
            exit(-1);
          }
        }
        else {
          numOfTargets[randomTargetNeuron]++;
        }

      } while( duplicateTarget || selfAssigned || limitTargets );

      matrixOfPostSynapticNeurons[n][s] = randomTargetNeuron;
    }
  }
#endif

#ifdef OUTFILE_CONNECTIONS
  ExportConnectionMatrixToFile( OUTFILE_CONNECTIONS );
#endif
#if __DEBUG__
  PrintMatrixOfPostSynapticNeurons();
#endif


  // = = = = = = = = = = = = = = = = = = = = = =
  // = WEIGHT MATRIX
  // = = = = = = = = = = = = = = = = = = = = = =
#if(GENERATE_NETWORK_FROM_EXTERNAL_DATA)
  ImportWeightMatrixFromFile(INFILE_WEIGHT_MATRIX);
#else
  for( int n = 0; n < NUM_TOTAL_NEURONS; ++n ) {
    for( int s = 0; s < NUM_SYNAPSES_PER_NEURON; ++s ) {

      if( n < NUM_EXCITATORY_NEURONS) {
  #if(INIT_WITH_RANDOM_WEIGHTS)
        double weight = GET_RANDOM_FLOAT(INIT_EXC_SYNAPTIC_WEIGHT);
        matrixOfSynapticWeights[n][s] = weight;
  #else
        matrixOfSynapticWeights[n][s] = INIT_EXC_SYNAPTIC_WEIGHT;
  #endif
      }
      else {
  #if(INIT_WITH_RANDOM_WEIGHTS)
        double weight = GET_RANDOM_FLOAT(INIT_INH_SYNAPTIC_WEIGHT);
        matrixOfSynapticWeights[n][s] = weight;
  #else
        matrixOfSynapticWeights[n][s] = INIT_INH_SYNAPTIC_WEIGHT;
  #endif
      }
      matrixOfSynapticWeights_derivatives[n][s] = 0.0;
    }
  }
#endif

#ifdef OUTFILE_WEIGHTS_INITIAL
  ExportWeightMatrixToFile(OUTFILE_WEIGHTS_INITIAL);
#endif
#if __DEBUG__
  PrintMatrixOfSynapticWeights();
#endif


  // = = = = = = = = = = = = = = = = = = = = = =
  // = DELAY MATRIX
  // = = = = = = = = = = = = = = = = = = = = = =
#if(GENERATE_NETWORK_FROM_EXTERNAL_DATA)
  ImportDelayMatrixFromFile(INFILE_DELAY_MATRIX);
#else
  for( int n = 0; n < NUM_TOTAL_NEURONS; ++n ) {
    int synapseCorrespondingToDelay = 0;
    for( int delayIdx = 0; delayIdx < MAX_SYNAPSE_DELAY; ++delayIdx ) {
      numEntriesPerDelay[n][delayIdx] = 0;
    }

    if( n < NUM_EXCITATORY_NEURONS ) {
      // For the excitatory neurons the delays are drawn from a uniform integer distribution.
      // e.g.:
      // delayIdx:     0      1      2      3      4      5       6 ..........      each slot corresponds to a delay value, i.e., 0 -> 1 ms, 1 -> 2 ms  etc.
      //               0, 1   2, 3   4, 5   6, 7   8, 9   9, 10   11, 12   ...      each slot has space for NUM_SYNAPSES_PER_NEURON (which is not needed) and
      //                                                                            holds the synapse numbers corresponding to that delay
      // REFACTOR COMMENT: In the Matlab implementation this is random.
      //                   This algorithm implicitly creates as many entries in delays as NUM_SYNAPSES_PER_NEURON.
      //                   This only works if NUM_SYNAPSES_PER_NEURON is a multiple of MAX_SYNAPSE_DELAY!
      for( int delayIdx = 0; delayIdx < MAX_SYNAPSE_DELAY; ++delayIdx ) {
        for( int synListIdx = 0; synListIdx < (NUM_SYNAPSES_PER_NEURON / MAX_SYNAPSE_DELAY); ++synListIdx ) {
          listOfSynsByNeuronAndDelay[n][delayIdx][synListIdx] = synapseCorrespondingToDelay++;
          numEntriesPerDelay[n][delayIdx]++;
        }
      }
    }
    else {
      // the inhibitory synapse delay is always 1 ms, thus, all synapses are added to delayIdx 0.
      for( int synListIdx = 0; synListIdx < NUM_SYNAPSES_PER_NEURON; ++synListIdx ) {
        listOfSynsByNeuronAndDelay[n][0][synListIdx] = synapseCorrespondingToDelay++;
        numEntriesPerDelay[n][0]++;
      }
    }
  }
#endif

  // for all neurons: find the excitatory presynaptic neurons and delays and store them in lists
  for( int n = 0; n < NUM_TOTAL_NEURONS; ++n ) {
    numPreSynapticNeuronsOfTarget[n] = 0;
  }

  for( int n = 0; n < NUM_TOTAL_NEURONS; ++n ) {
    for( int excPreSynNeuron = 0; excPreSynNeuron < NUM_EXCITATORY_NEURONS; ++excPreSynNeuron ) {
      for( int s = 0; s < NUM_SYNAPSES_PER_NEURON; ++s ) {
        int targetNeuron = matrixOfPostSynapticNeurons[excPreSynNeuron][s];

        if( targetNeuron == n ) {
          int idx = numPreSynapticNeuronsOfTarget[n];
          listOfPresynapticNeurons[n][idx] = excPreSynNeuron;

          // determine the delay of a connection and store it in a list
          for( int delayIdx = 0; delayIdx < MAX_SYNAPSE_DELAY; ++delayIdx ) {
            int excPreSynNeuron_numEntriesOfDelay = numEntriesPerDelay[excPreSynNeuron][delayIdx];

            for( int i = 0; i < excPreSynNeuron_numEntriesOfDelay; ++i ) {
              int synapse = listOfSynsByNeuronAndDelay[excPreSynNeuron][delayIdx][i];
              int postSynNeuron = matrixOfPostSynapticNeurons[excPreSynNeuron][synapse];

              if( postSynNeuron == n ) {
                int idx = numPreSynapticNeuronsOfTarget[n];
                listOfPresynapticDelays[n][idx] = delayIdx;
                // REFACTOR COMMENT: delayIdx represents the delay value, where 0 represents 1ms, thus,
                //                   one needs to distinguish between a 1ms delay and an empty entry.
                //                   This is done by tracking the number of entries in delays_numEntriesPerDelay.
              }
            }
          }

          // maintain list of pointers into the weight matrices; used for weight update
          {
            int idx = numPreSynapticNeuronsOfTarget[n];
            // REFACTOR COMMENT: dead code removed
            // listOfPointersToSynapticWeights[n][idx] = &matrixOfSynapticWeights[excPreSynNeuron][s];
            listOfPointersToSynapticWeights_derivatives[n][idx] = &matrixOfSynapticWeights_derivatives[excPreSynNeuron][s];
            idx++;
            numPreSynapticNeuronsOfTarget[n] = idx;
          }
        }
      }
    }
  }

#if __DEBUG__
  PrintMatrixOfSynapticDelays();
#endif

#ifdef OUTFILE_DELAYS
  ExportDelayMatrixToFile( OUTFILE_DELAYS );
#endif

  // REFACTOR COMMENT: the array is much larger and remains uninitialized
  for( int n = 0; n < NUM_TOTAL_NEURONS; ++n ) {
    for( int d = 0; d < 1 + MAX_SYNAPSE_DELAY; ++d ) {
      LTP[n][d] = 0.0;
    }
  }

  for( int n = 0; n < NUM_TOTAL_NEURONS; ++n ) {
    LTD[n] = 0.0;
  }

  // REFACTOR COMMENT: just for the algorithm; does not contribute to the network activity
  numFirings = 1;                                                                        // a dummy spike ...
  firings[0][IDX_TIME] = -MAX_SYNAPSE_DELAY;                                             // ... at -MAX_SYNAPSE_DELAY for ...
  firings[0][IDX_NEURON] = 0;                                                            // ... neuron n = 0

#ifdef OUTFILE_CON_WEIGHT_DELAY_INITIAL
  ExportConnectionMatrixWeightAndDelay(OUTFILE_CON_WEIGHT_DELAY_INITIAL);
#endif

#ifdef OUTFILE_HNC_NODE_CONNECT_CALLS
  ExportHNCNodeConnectCalls(OUTFILE_HNC_NODE_CONNECT_CALLS);
#endif

#ifdef OUTFILE_NEST_CONNECTIVITY_DICTS
  ExportNESTConnectionDicts(OUTFILE_NEST_CONNECTIVITY_DICTS);
#endif
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = Initialize simulation: set seed, delete previously created files, open files
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void InitializeSimulation() {
  srand( 0 );                                                                            // for repeatability

#if( LOG_MIN_MAX_V_U )
  v_max = 0.0;
  v_min = 0.0;
  u_max = 0.0;
  u_min = 0.0;
#endif

#ifdef OUTFILE_STIMULUS
  DeleteFile( OUTFILE_STIMULUS );
#endif

#ifdef OUTFILE_CONNECTIONS
  DeleteFile( OUTFILE_CONNECTIONS );
#endif

#ifdef OUTFILE_DELAYS
  DeleteFile( OUTFILE_DELAYS );
#endif

#ifdef OUTFILE_CON_WEIGHT_DELAY_INITIAL
  DeleteFile(OUTFILE_CON_WEIGHT_DELAY_INITIAL);
#endif

#ifdef OUTFILE_HNC_NODE_CONNECT_CALLS
  DeleteFile(OUTFILE_HNC_NODE_CONNECT_CALLS);
#endif

#ifdef OUTFILE_NEST_CONNECTIVITY_DICTS
  DeleteFile(OUTFILE_NEST_CONNECTIVITY_DICTS);
#endif

#ifdef OUTFILE_WEIGHTS_INITIAL
  DeleteFile(OUTFILE_WEIGHTS_INITIAL);
#endif

  DeleteFile( OUTFILE_FIRINGS );

#ifdef OUTFILE_STIMULUS
  pFileStimulusOutput = fopen( OUTFILE_STIMULUS, "w" );
  if( pFileStimulusOutput == nullptr ) {
    printf( "[ERROR] Failed to open stimulus output file. Filename: %s\n", OUTFILE_STIMULUS );
    exit( RC_ERROR_EXIT );
  }
#endif

#if(USE_EXTERNAL_STIMULUS)
  pFileStimulusInput = fopen(INFILE_STIMULUS,"r" );
  if(pFileStimulusInput == nullptr) {
    printf( "[ERROR] Failed to open stimulus input file. Filename: %s\n", INFILE_STIMULUS );
    exit( RC_ERROR_EXIT );
  }
#endif
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = Finalize simulation: free resources
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void FinalizeSimulation() {
  if( pFileStimulusInput != nullptr ) fclose( pFileStimulusInput );
  if( pFileStimulusOutput != nullptr ) fclose( pFileStimulusOutput );
}


#else                                                                                // original version
  // original version of the Izhikevich polychronization model available for download at
  // https://www.izhikevich.org/publications/spnet.htm
  #include "poly_spnet.cpp"
#endif
