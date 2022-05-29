/*
 *  utils.cpp
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

#include "params.h"

#if( __RUN_REFACTORED_VERSION__ )

#define __IMPORT_GLOBAL_VARIABLES__ true

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <cstring>

#include "globals.h"
#include "utils.h"


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   E X P O R T   N E T W O R K   S T A T E S
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = Export connection matrix to file in printable ASCII format.
// =
// =                      001 002                     100         <- synapse
// = Source neuron 0001   nnn nnn nnn nnn nnn nnn ... nnn
// = Source neuron 0002   nnn nnn nnn nnn nnn nnn ... nnn
// = ...
// = Source neuron 1000   nnn nnn nnn nnn nnn nnn ... nnn
// =
// = nnn ... target neuron id
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void ExportConnectionMatrixToFile( const char *pFileName ) {
  FILE *pFile = fopen( pFileName, "w" );
  if( pFile == nullptr ) {
    printf( "[ERROR] Failed to open file: %s   (errno: %d)\n", pFileName, errno);
    exit(RC_ERROR_EXIT);
  }

  printf( "[INFO]  Export connection matrix to file: %s\n", pFileName );

  for( int n = 0; n < NUM_TOTAL_NEURONS; ++n ) {
    for( int s = 0; s < NUM_SYNAPSES_PER_NEURON; ++s ) {
      fprintf( pFile, "%03d ", matrixOfPostSynapticNeurons[n][s] );
    }
    fprintf( pFile, "\n" );
  }

  fclose( pFile );
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = Export weight matrix to file in printable ASCII format.
// =
// =                      001     002                         100        <- synapse
// = Source neuron 0001   ww.wwww ww.wwww ww.wwww ww.wwww ... ww.wwww
// = Source neuron 0002   ww.wwww ww.wwww ww.wwww ww.wwww ... ww.wwww
// = ...
// = Source neuron 1000   ww.wwww ww.wwww ww.wwww ww.wwww ... ww.wwww
// =
// = ww.wwww ... synaptic strength
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void ExportWeightMatrixToFile( const char *pFileName ) {
  FILE *pFile = fopen( pFileName, "w" );
  if( pFile == nullptr ) {
    printf( "[ERROR] Failed to open file: %s   (errno: %d)\n", pFileName, errno);
    exit(RC_ERROR_EXIT);
  }

  printf( "[INFO]  Export weight matrix to file: %s\n", pFileName );

  for( int n = 0; n < NUM_TOTAL_NEURONS; ++n ) {
    for( int s = 0; s < NUM_SYNAPSES_PER_NEURON; ++s ) {
      fprintf( pFile, "%f ", matrixOfSynapticWeights[n][s] );
    }
    fprintf( pFile, "\n" );
  }

  fclose( pFile );
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = Export delay matrix to file in printable ASCII format.
// =
// =                      001 002                     100   <- synapse
// = Source neuron 0001   ddd ddd ddd ddd ddd ddd ... ddd
// = Source neuron 0002   ddd ddd ddd ddd ddd ddd ... ddd
// = ...
// = Source neuron 1000   ddd ddd ddd ddd ddd ddd ... ddd
// =
// = ddd ... conduction delay
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void ExportDelayMatrixToFile( const char *pFileName ) {
  FILE *pFile = fopen( pFileName, "w" );
  if( pFile == nullptr ) {
    printf( "[ERROR] Failed to open file: %s   (errno: %d)\n", pFileName, errno);
    exit(RC_ERROR_EXIT);
  }

  printf( "[INFO]  Export delay matrix to file: %s\n", pFileName );

  for( int n = 0; n < NUM_TOTAL_NEURONS; ++n ) {
    for( int s = 0; s < NUM_SYNAPSES_PER_NEURON; ++s ) {
      int delayValue = GetDelayOfConnection( n, s );
      fprintf( pFile, "%03d ", delayValue );
    }
    fprintf( pFile, "\n" );
  }

  fclose( pFile );
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = Export connection, weight, and delay matrix to file for import into PyNN.
// =
// =
// = exc_exc_connections = [
// =   (   source neuron, target neuron, weight, delay ), ...,
// =   ... , (   source neuron, target neuron, weight, delay )
// =                      ]
// =
// = exc_inh_connections = [
// =   (   source neuron, target neuron, weight, delay ), ...,
// =   ... , (   source neuron, target neuron, weight, delay )
// =                      ]
// = inh_exc_connections = [
// =   (   source neuron, target neuron, weight, delay ), ...,
// =   ... , (   source neuron, target neuron, weight, delay )
// =                      ]
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void ExportConnectionMatrixWeightAndDelay( const char *pFileName ) {
  const int CONST_ENTRIES_PER_LINE = 8;

  FILE *pFile = fopen( pFileName, "w" );
  if( pFile == nullptr ) {
    printf( "[ERROR] Failed to open file: %s   (errno: %d)\n", pFileName, errno);
    exit(RC_ERROR_EXIT);
  }

  printf( "[INFO]  Export connection, weight and delay matrix for import into PyNN to file: %s\n", pFileName );

  // excitatory to excitatory connections
  fprintf( pFile, "exc_exc_connections = [\n" );
  int entriesCount = 0;
  for( int excPreSynNeuron = 0; excPreSynNeuron < NUM_EXCITATORY_NEURONS; ++excPreSynNeuron ) {
    for( int synapse = 0; synapse < NUM_SYNAPSES_PER_NEURON; ++synapse ) {
      if( matrixOfPostSynapticNeurons[excPreSynNeuron][synapse] < NUM_EXCITATORY_NEURONS ) {

        int excPostSynNeuron = matrixOfPostSynapticNeurons[excPreSynNeuron][synapse];
        double weight = fabs( matrixOfSynapticWeights[excPreSynNeuron][synapse] );   // PyNN expects positive weight values for inhibitory synapses
        double delay = GetDelayOfConnection( excPreSynNeuron, synapse );

        fprintf( pFile, "( %3d, %3d, %f, %f ), ", excPreSynNeuron, excPostSynNeuron, weight, delay );

        entriesCount++;
        if( entriesCount == CONST_ENTRIES_PER_LINE ) {
          fprintf( pFile, "\n" );
          entriesCount = 0;
        }
      }
    }
  }
  fprintf( pFile, "]\n\n" );

  // excitatory to inhibitory connections
  fprintf( pFile, "exc_inh_connections = [\n" );
  entriesCount = 0;
  for( int excPreSynNeuron = 0; excPreSynNeuron < NUM_EXCITATORY_NEURONS; ++excPreSynNeuron ) {
    for( int synapse = 0; synapse < NUM_SYNAPSES_PER_NEURON; ++synapse ) {
      if( matrixOfPostSynapticNeurons[excPreSynNeuron][synapse] >= NUM_EXCITATORY_NEURONS) {

        int inhPostSynNeuron = matrixOfPostSynapticNeurons[excPreSynNeuron][synapse];
        double weight = fabs( matrixOfSynapticWeights[excPreSynNeuron][synapse] );   // PyNN expects positive weight values for inhibitory synapses
        double delay = GetDelayOfConnection( excPreSynNeuron, synapse );

        fprintf( pFile, "( %3d, %3d, %f, %f ), "
               , excPreSynNeuron, inhPostSynNeuron - NUM_EXCITATORY_NEURONS, weight, delay );

        entriesCount++;
        if( entriesCount == CONST_ENTRIES_PER_LINE ) {
          fprintf( pFile, "\n" );
          entriesCount = 0;
        }
      }
    }
  }
  fprintf( pFile, "]\n\n" );

  // inhibitory to excitatory connections
  entriesCount = 0;
  fprintf( pFile, "inh_exc_connections = [\n" );
  for( int inhPreSynNeuron = NUM_EXCITATORY_NEURONS; inhPreSynNeuron < NUM_TOTAL_NEURONS; ++inhPreSynNeuron ) {
    for( int synapse = 0; synapse < NUM_SYNAPSES_PER_NEURON; ++synapse ) {
      if( matrixOfPostSynapticNeurons[inhPreSynNeuron][synapse] < NUM_EXCITATORY_NEURONS) {

        int excPostSynNeuron = matrixOfPostSynapticNeurons[inhPreSynNeuron][synapse];
        double weight = fabs( matrixOfSynapticWeights[inhPreSynNeuron][synapse] );   // PyNN expects positive weight values for inhibitory synapses
        double delay = GetDelayOfConnection( inhPreSynNeuron, synapse );

        fprintf( pFile, "( %3d, %3d, %6.3f, %6.3f ), "
               , inhPreSynNeuron - NUM_EXCITATORY_NEURONS, excPostSynNeuron, weight, delay );

        entriesCount++;
        if( entriesCount == CONST_ENTRIES_PER_LINE ) {
          fprintf( pFile, "\n" );
          entriesCount = 0;
        }
      }
    }
  }
  fprintf( pFile, "]\n\n" );

  fclose( pFile );
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = Export connectivity in HNC node format.
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
//
//   Create a list of C-API function calls.
//
//   Connect(preSynNeuronId, postSynNeuronId, weight, delay(in units of 0.1 ms steps));
//   Connect( ... );
//   ...
//   Connect( ... );
//
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void ExportHNCNodeConnectCalls( const char *pFileName) {
  FILE *pFile = fopen( pFileName, "w" );
  if( pFile == nullptr ) {
    printf( "[ERROR] Failed to open file: %s   (errno: %d)\n", pFileName, errno);
    exit(RC_ERROR_EXIT);
  }
  printf( "[INFO]  Export network connectivity data in HNC node format: %s\n", pFileName );

  fprintf( pFile, "#ifndef __ACP_BUILD_NETWORK_H__\n");
  fprintf( pFile, "#define __ACP_BUILD_NETWORK_H__\n");
  fprintf( pFile, "#define RS_Iext_EX  (float)(0.0)\n");
  fprintf( pFile, "#define FS_Iext_INH (float)(0.0)\n");
  fprintf( pFile, "// + + + Two-population Izhikevich Network + + +\n");
  fprintf( pFile, "  void BuildNetwork() {\n");
  fprintf( pFile, "  Divider_Enable();\n");
  fprintf( pFile, "  // Neuron Setup\n");
  // neurons
  for( int neuron = 0; neuron < NUM_EXCITATORY_NEURONS; ++neuron ) {
    fprintf( pFile, "  Izhk_SetNeuronStateVars(%d, IZHK_NEURON_TYPE_RS, %+09.4f, %+09.4f, RS_Iext_EX);\n", neuron, RS_V_INIT, RS_U_INIT);
  }
  fprintf( pFile, "\n");
  for( int neuron = NUM_EXCITATORY_NEURONS; neuron < NUM_EXCITATORY_NEURONS + NUM_INHIBITORY_NEURONS; ++neuron ) {
    fprintf( pFile, "  Izhk_SetNeuronStateVars(%d, IZHK_NEURON_TYPE_FS, %+09.4f, %+09.4f, FS_Iext_INH);\n", neuron, FS_V_INIT, FS_U_INIT);
  }

  fprintf( pFile, "\n");

  fprintf( pFile, "  // Connection Setup\n");
  // connections
  for( int preSynNeuron = 0; preSynNeuron < (NUM_EXCITATORY_NEURONS + NUM_INHIBITORY_NEURONS); ++preSynNeuron ) {
    for( int synapse = 0; synapse < NUM_SYNAPSES_PER_NEURON; ++synapse ) {
      int    postSynNeuron = matrixOfPostSynapticNeurons[preSynNeuron][synapse];
      double weight        = matrixOfSynapticWeights[preSynNeuron][synapse];
      double delay         = GetDelayOfConnection( preSynNeuron, synapse );        // is in ms
      unsigned int delay_01msResolution = (unsigned int)( delay * 10);             // in 0.1 ms
      if(weight == 0) {
        fprintf( pFile, "  // Connect( %d, %d, %+09.4f, %d );\n", preSynNeuron, postSynNeuron, weight, delay_01msResolution );
      }
      else {
        fprintf( pFile, "  Connect( %d, %d, %+09.4f, %d );\n", preSynNeuron, postSynNeuron, weight, delay_01msResolution );
      }
    }
    fprintf( pFile, "\n");
  }

  fprintf( pFile, "  NeuronStates_BramLoad();\n");
  fprintf( pFile, "  return;\n");
  fprintf( pFile, "}\n");
  fprintf( pFile, "#endif // __ACP_BUILD_NETWORK_H__\n");

  fclose( pFile );
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = Export connectivity for import into NEST.
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// export as dictionaries, all lists in a dictionary have to have same length
// con_exc_to_all = {'source': [1,2,3,4,..,n], 'target': [1,2,3,4,..,n], 'weight': [0.1,0.2,..,n], 'delay': [0.1,0.2,..,n]}
// con_inh_to_exc = {'source': [1,2,3,4,..,n], 'target': [1,2,3,4,..,n], 'weight': [-0.1,-0.2,..,n], 'delay': [0.1,0.2,..,n]}
//
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void ExportNESTConnectionDicts(const char* pFileName) {

  typedef struct connection_ connection_t;
  struct     connection_ {
      int    preSynNeuron;
      int    postSynNeuron;
      double weight;
      double delay;
  };

  printf( "[INFO]  Export network connectivity data for import into NEST: %s\n", pFileName );

  // neurons
  // ... nothing to do ...     // no export of state variable initial values needed, this is done in NEST

  // collect the connections in two array
  connection_t* connMatrix_exc_to_all = (connection_t*)malloc(NUM_EXCITATORY_NEURONS * NUM_SYNAPSES_PER_NEURON * sizeof(connection_t));
  connection_t* connMatrix_inh_to_exc = (connection_t*)malloc(NUM_INHIBITORY_NEURONS * NUM_SYNAPSES_PER_NEURON * sizeof(connection_t));

  int connIdx = 0;
  for( int preSynNeuron = 0; preSynNeuron < NUM_EXCITATORY_NEURONS; ++preSynNeuron ) {
    for( int synapse = 0; synapse < NUM_SYNAPSES_PER_NEURON; ++synapse ) {
      connMatrix_exc_to_all[connIdx].preSynNeuron = preSynNeuron + 1;                                            // NEST numbers neurons starting from 1
      connMatrix_exc_to_all[connIdx].postSynNeuron = matrixOfPostSynapticNeurons[preSynNeuron][synapse] + 1;     // NEST numbers neurons starting from 1
      connMatrix_exc_to_all[connIdx].weight = matrixOfSynapticWeights[preSynNeuron][synapse];
      connMatrix_exc_to_all[connIdx].delay = GetDelayOfConnection( preSynNeuron, synapse );                      // in ms
      connIdx++;
    }
  }
  connIdx = 0;
  for( int preSynNeuron = NUM_EXCITATORY_NEURONS; preSynNeuron < (NUM_EXCITATORY_NEURONS + NUM_INHIBITORY_NEURONS); ++preSynNeuron ) {
    for( int synapse = 0; synapse < NUM_SYNAPSES_PER_NEURON; ++synapse ) {
      connMatrix_inh_to_exc[connIdx].preSynNeuron = preSynNeuron + 1;                                            // NEST numbers neurons starting from 1
      connMatrix_inh_to_exc[connIdx].postSynNeuron = matrixOfPostSynapticNeurons[preSynNeuron][synapse] + 1;     // NEST numbers neurons starting from 1
      connMatrix_inh_to_exc[connIdx].weight = matrixOfSynapticWeights[preSynNeuron][synapse];
      connMatrix_inh_to_exc[connIdx].delay = GetDelayOfConnection( preSynNeuron, synapse );                      // in ms
      connIdx++;
    }
  }

  // export to file

  FILE *pFile = fopen( pFileName, "w" );
  if( pFile == nullptr ) {
    printf( "[ERROR] Failed to open file: %s   (errno: %d)\n", pFileName, errno);
    exit(RC_ERROR_EXIT);
  }

  // exc to all
  fprintf( pFile, "con_exc_to_all = { \\\n");
  fprintf( pFile, "'source': [");
  for( int i = 0; i < NUM_EXCITATORY_NEURONS * NUM_SYNAPSES_PER_NEURON - 1; ++i) {
    fprintf(pFile, "%d, ", connMatrix_exc_to_all[i].preSynNeuron);
    if(i > 0 && i % 20 == 0) fprintf(pFile, " \\\n");
  }
  fprintf(pFile, "%d], \\\n", connMatrix_exc_to_all[NUM_EXCITATORY_NEURONS * NUM_SYNAPSES_PER_NEURON - 1].preSynNeuron);

  fprintf( pFile, "'target': [");
  for( int i = 0; i < NUM_EXCITATORY_NEURONS * NUM_SYNAPSES_PER_NEURON - 1; ++i) {
    fprintf(pFile, "%d, ", connMatrix_exc_to_all[i].postSynNeuron);
    if(i > 0 && i % 20 == 0) fprintf(pFile, " \\\n");
  }
  fprintf(pFile, "%d], \\\n", connMatrix_exc_to_all[NUM_EXCITATORY_NEURONS * NUM_SYNAPSES_PER_NEURON - 1].postSynNeuron);

  fprintf( pFile, "'weight': [");
  for( int i = 0; i < NUM_EXCITATORY_NEURONS * NUM_SYNAPSES_PER_NEURON - 1; ++i) {
    fprintf(pFile, "%.6f, ", connMatrix_exc_to_all[i].weight);
    if(i > 0 && i % 20 == 0) fprintf(pFile, " \\\n");
  }
  fprintf(pFile, "%.6f], \\\n", connMatrix_exc_to_all[NUM_EXCITATORY_NEURONS * NUM_SYNAPSES_PER_NEURON - 1].weight);

  fprintf( pFile, "'delay': [");
  for( int i = 0; i < NUM_EXCITATORY_NEURONS * NUM_SYNAPSES_PER_NEURON - 1; ++i) {
    fprintf(pFile, "%.6f, ", connMatrix_exc_to_all[i].delay);
    if(i > 0 && i % 20 == 0) fprintf(pFile, " \\\n");
  }
  fprintf(pFile, "%.6f]} \n", connMatrix_exc_to_all[NUM_EXCITATORY_NEURONS * NUM_SYNAPSES_PER_NEURON - 1].delay);

  fprintf(pFile, "\n");

  // inh to exc
  fprintf( pFile, "con_inh_to_exc = { \\\n");
  fprintf( pFile, "'source': [");
  for( int i = 0; i < NUM_INHIBITORY_NEURONS * NUM_SYNAPSES_PER_NEURON - 1; ++i) {
    fprintf(pFile, "%d, ", connMatrix_inh_to_exc[i].preSynNeuron);
    if(i > 0 && i % 20 == 0) fprintf(pFile, " \\\n");
  }
  fprintf(pFile, "%d], \\\n", connMatrix_inh_to_exc[NUM_INHIBITORY_NEURONS * NUM_SYNAPSES_PER_NEURON - 1].preSynNeuron);

  fprintf( pFile, "'target': [");
  for( int i = 0; i < NUM_INHIBITORY_NEURONS * NUM_SYNAPSES_PER_NEURON - 1; ++i) {
    fprintf(pFile, "%d, ", connMatrix_inh_to_exc[i].postSynNeuron);
    if(i > 0 && i % 20 == 0) fprintf(pFile, " \\\n");
  }
  fprintf(pFile, "%d], \\\n", connMatrix_inh_to_exc[NUM_INHIBITORY_NEURONS * NUM_SYNAPSES_PER_NEURON - 1].postSynNeuron);

  fprintf( pFile, "'weight': [");
  for( int i = 0; i < NUM_INHIBITORY_NEURONS * NUM_SYNAPSES_PER_NEURON - 1; ++i) {
    fprintf(pFile, "%.6f, ", connMatrix_inh_to_exc[i].weight);
    if(i > 0 && i % 20 == 0) fprintf(pFile, " \\\n");
  }
  fprintf(pFile, "%.6f], \\\n", connMatrix_inh_to_exc[NUM_INHIBITORY_NEURONS * NUM_SYNAPSES_PER_NEURON - 1].weight);

  fprintf( pFile, "'delay': [");
  for( int i = 0; i < NUM_INHIBITORY_NEURONS * NUM_SYNAPSES_PER_NEURON - 1; ++i) {
    fprintf(pFile, "%.6f, ", connMatrix_inh_to_exc[i].delay);
    if(i > 0 && i % 20 == 0) fprintf(pFile, " \\\n");
  }
  fprintf(pFile, "%.6f]} \n", connMatrix_inh_to_exc[NUM_INHIBITORY_NEURONS * NUM_SYNAPSES_PER_NEURON - 1].delay);

  fclose( pFile );
  free(connMatrix_exc_to_all);
  free(connMatrix_inh_to_exc);
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   I M P O R T   N E T W O R K   S T A T E S
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = Import the connection matrix from file generated with
// = ExportConnectionMatrixToFile().
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void ImportConnectionMatrixFromFile( const char *pFileName ) {
  FILE *pFile = fopen( pFileName, "r" );
  if( pFile == nullptr ) {
    printf( "[ERROR] Failed to open file: %s   (errno: %d)\n", pFileName, errno);
    exit(RC_ERROR_EXIT);
  }

  printf( "[INFO]  Import connection matrix from file: %s\n", pFileName );

  char line[NUM_SYNAPSES_PER_NEURON * SIZE_OF_ENTRY_NEURON] = {};
  int targetNeuron = 0;
  char *pToken = nullptr;

  for( int n = 0; n < NUM_TOTAL_NEURONS; ++n ) {
    fgets( line, sizeof( line ), pFile );

    // first token, i.e. synapse 0 target neuron
    pToken = strtok( line, " " );
    targetNeuron = atoi( pToken );
    matrixOfPostSynapticNeurons[n][0] = targetNeuron;

    // subsequent tokens, i.e. synapses 1 .. s
    for( int s = 1; s < NUM_SYNAPSES_PER_NEURON; ++s ) {
      pToken = strtok( NULL, " " );
      targetNeuron = atoi( pToken );
      matrixOfPostSynapticNeurons[n][s] = targetNeuron;
    }
  }

  fclose( pFile );
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = Import the weight matrix from file generated with
// = ExportWeightMatrixToFile().
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void ImportWeightMatrixFromFile( const char *pFileName ) {
  FILE *pFile = fopen( pFileName, "r" );
  if( pFile == nullptr ) {
    printf( "[ERROR] Failed to open file: %s   (errno: %d)\n", pFileName, errno );
    exit( RC_ERROR_EXIT );
  }

  printf( "[INFO]  Import weight matrix from file: %s\n", pFileName );

  char line[NUM_SYNAPSES_PER_NEURON * SIZE_OF_ENTRY_WEIGHT] = {};
  double weight = 0.0;
  char *pToken = nullptr;

  for( int n = 0; n < NUM_TOTAL_NEURONS; ++n ) {
    fgets( line, sizeof( line ), pFile );

    // first token, i.e. weight of connection synapse 0 of target neuron n
    pToken = strtok( line, " " );
    weight = atof( pToken );
    matrixOfSynapticWeights[n][0] = weight;

    // subsequent tokens, i.e. weights of connections synapses s of target neuron n
    for( int s = 1; s < NUM_SYNAPSES_PER_NEURON; ++s ) {
      pToken = strtok( NULL, " " );
      weight = atof( pToken );
      matrixOfSynapticWeights[n][s] = weight;
    }
  }

  fclose( pFile );
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = Import the delay matrix from file generated with
// = ExportDelayMatrixToFile().
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void ImportDelayMatrixFromFile( const char *pFileName ) {
  FILE *pFile = fopen( pFileName, "r" );
  if( pFile == nullptr ) {
    printf( "[ERROR] Failed to open file: %s   (errno: %d)\n", pFileName, errno );
    exit( RC_ERROR_EXIT );
  }

  printf( "[INFO]  Import delay matrix from file: %s\n", pFileName );

  char line[NUM_SYNAPSES_PER_NEURON * SIZE_OF_ENTRY_DELAY] = {};
  char *pToken = nullptr;

  for( int n = 0; n < NUM_TOTAL_NEURONS; ++n ) {
    fgets( line, sizeof( line ), pFile );

    // first token, i.e. delay of synapse 0
    pToken = strtok( line, " " );
    int delayIdx = atoi( pToken ) - 1;

    int idx = numEntriesPerDelay[n][delayIdx];
    listOfSynsByNeuronAndDelay[n][delayIdx][idx] = 0;                        // synapse 0
    numEntriesPerDelay[n][delayIdx]++;

    // subsequent tokens, i.e. delays of synapses s
    for( int s = 1; s < NUM_SYNAPSES_PER_NEURON; ++s ) {
      pToken = strtok( NULL, " " );
      delayIdx = atoi( pToken ) - 1;

      idx = numEntriesPerDelay[n][delayIdx];
      listOfSynsByNeuronAndDelay[n][delayIdx][idx] = s;
      numEntriesPerDelay[n][delayIdx]++;
    }
  }

  fclose( pFile );
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   R E C O R D   E X T E R N A L   S T I M U L U S
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = Record the random input data the network is stimulated with.
// =
// = second  millisecond  id of the neuron that receives the input
// = ssssss  mmm          nnn
// = ...
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void RecordRandomStimulusToFile( const char *pFileName, int simTimeSecond, int simTimeMillisecond, int inputNeuron ) {
  if( pFileStimulusOutput == nullptr ) {
    printf( "[ERROR] File not open: %s   (errno: %d)\n", pFileName, errno );
    exit( RC_ERROR_EXIT );
  }

  if( inputNeuron > 0 ) {
    fprintf( pFileStimulusOutput, "%06d %03d %03d\n", simTimeSecond, simTimeMillisecond, inputNeuron );
  }
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   R E A D   E X T E R N A L   S T I M U L U S
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = Read next record from external stimulus file and return the neuron id.
// =
// = second  millisecond  id of the neuron that receives the input
// = ssssss  mmm          nnn
// = ...
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
int GetNextExternalStimulusFromFile( const char *pFileName, int simTimeSecond, int t ) {
  char line[32] = {};
  char *pToken = nullptr;

  if( pFileStimulusInput == nullptr ) {
    printf( "[ERROR] File not open: %s   (errno: %d)\n", pFileName, errno );
    exit( RC_ERROR_EXIT );
  }

  if( !fgets( line, sizeof( line ), pFileStimulusInput )) {
    printf( "[ERROR] Unexpected end of file while reading stimulus data. Simulation time: %ds %dms  Filename: %s\n"
          , simTimeSecond, t, pFileName );
    exit( RC_ERROR_EXIT );
  }

  pToken = strtok( line, " " );
  if( pToken == nullptr ) {
    return RC_NOID;
  }

  pToken = strtok( nullptr, " " );
  if( pToken == nullptr ) {
    return RC_NOID;
  }

  pToken = strtok( nullptr, " " );
  if( pToken == nullptr ) {
    return RC_NOID;
  }
  int neuronId = atoi( pToken );

  return neuronId;
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   R E C O R D   N E T W O R K   A C T I V I T Y   D A T A
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = Write the spike times to file in a printable ASCII format.
// =
// = second  millisecond  id of the neuron that has fired
// = ssssss  mmm          nnn
// = ...
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void RecordNetworkActivityToFile( const char *pFileName, int simulationSecond, int numFirings ) {
  FILE *pFile = fopen( pFileName, "a+" );
  if( pFile == nullptr ) {
    printf( "[ERROR] Failed to open file: %s   (errno: %d)\n", pFileName, errno );
    exit( RC_ERROR_EXIT );
  }

  // skip negative times
  int idx = 0;
  while( firings[idx][IDX_TIME] < 0 ) {
    idx++;
  }

  for( ; idx < numFirings; ++idx ) {
    fprintf( pFile, "%06d %03d %03d\n", simulationSecond, firings[idx][IDX_TIME], firings[idx][IDX_NEURON] );
  }

  fclose( pFile );
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   H E L P E R   F U N C T I O N S
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = Return the delay value of a connection.
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
int GetDelayOfConnection( int preSynNeuron, int synapse ) {
  int delay = 0;

  for( int delayIdx = 0; delayIdx < MAX_SYNAPSE_DELAY; ++delayIdx ) {
    int preSynNeuron_numEntriesOfDelay = numEntriesPerDelay[preSynNeuron][delayIdx];
    for( int i = 0; i < preSynNeuron_numEntriesOfDelay; ++i ) {
      if( synapse == listOfSynsByNeuronAndDelay[preSynNeuron][delayIdx][i] ) {
        delay = delayIdx + 1;                                                  // index 0 corresponds to 1ms delay
      }
    }
  }

  return delay;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = Delete a file from disk.
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void DeleteFile( const char *fileName ) {
  if( remove( fileName ) == 0 ) {
    printf( "[INFO]  File deleted: %s\n", fileName );
  }
  else {
    printf( "[INFO]  File could not be deleted: %s\n", fileName );
  }
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   D E B U G   F U N C T I O N S
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = Debugprint: connection matrix
// =
// = Matrix of Post Synaptic Neurons
// =
// =          0    1    2    3    4    5    6    7    8    9  ...  <- synapse
// = -------------------------------------------------------- ...
// =  0 :   840  394  783  798  911  197  335  768  277  553  ...
// =  1 :   950  920  147  881  641  431  619  281  786  307  ...  <- target neuron id
// =  2 :    76  649  248  629  229  700  316  328  231   74  ...
// =  3 :   291  180  684  727  139  603  492  838  724  178  ...
// =  4 :    98  923  169  481  225  826  290  357  878  344  ...
// =  5 :   324  874  589  637  759  775  794  262  604  470  ...
// =  6 :   887  933  173  447  487  795  639  965  155  292  ...
// =  7 :   347  205  522  400  307  679  645  443  269  703  ...
// =  8 :   452  160  308  433    5  649  126  461   84  780  ...
// =  9 :   805  749  398  366  394  272  599   68  901  432  ...
// = 10 :    20   53  897  899   39  419  183  219  778  622  ...
// = ....
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void PrintMatrixOfPostSynapticNeurons() {
  printf( "     Matrix of Post Synaptic Neurons\n\n" );
  printf( "        " );
  for( int s = 0; s < NUM_SYNAPSES_PER_NEURON; ++s ) {
    printf( " %3d ", s );
  }
  printf( "\n" );
  printf( "-------" );
  for( int s = 0; s < NUM_SYNAPSES_PER_NEURON; ++s ) {
    printf( "-----" );
  }
  printf( " \n" );
  for( int n = 0; n < NUM_TOTAL_NEURONS; ++n ) {
    printf( "%3d :   ", n );
    for( int s = 0; s < NUM_SYNAPSES_PER_NEURON; ++s ) {
      printf( " %3d ", matrixOfPostSynapticNeurons[n][s] );
    }
    printf( "\n" );
  }
  printf( "Press enter ...\n" );
  getchar();
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = Debugprint: weight matrix
// =
// = Matrix of Synaptic Weights
// =
// =                0         1         2         3         4  ...  <- synapse
// = --------------------------------------------------------- ...
// =  0 :    6.000000  6.000000  6.000000  6.000000  6.000000  ...
// =  1 :    6.000000  6.000000  6.000000  6.000000  6.000000  ...  <- connection strength
// =  2 :    6.000000  6.000000  6.000000  6.000000  6.000000  ...
// =  3 :    6.000000  6.000000  6.000000  6.000000  6.000000  ...
// =  4 :    6.000000  6.000000  6.000000  6.000000  6.000000  ...
// =  5 :    6.000000  6.000000  6.000000  6.000000  6.000000  ...
// =  6 :    6.000000  6.000000  6.000000  6.000000  6.000000  ...
// =  7 :    6.000000  6.000000  6.000000  6.000000  6.000000  ...
// =  8 :    6.000000  6.000000  6.000000  6.000000  6.000000  ...
// =  9 :    6.000000  6.000000  6.000000  6.000000  6.000000  ...
// = 10 :    6.000000  6.000000  6.000000  6.000000  6.000000  ...
// = ....
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void PrintMatrixOfSynapticWeights() {
  printf( "     Matrix of Synaptic Weights\n\n" );
  printf( "        " );
  for( int s = 0; s < NUM_SYNAPSES_PER_NEURON; ++s ) {
    printf( " %8d ", s );
  }
  printf( "\n" );
  printf( "-------" );
  for( int s = 0; s < NUM_SYNAPSES_PER_NEURON; ++s ) {
    printf( "----------" );
  }
  printf( " \n" );
  for( int n = 0; n < NUM_TOTAL_NEURONS; ++n ) {
    printf( "%3d :   ", n );
    for( int s = 0; s < NUM_SYNAPSES_PER_NEURON; ++s ) {
      printf( " %f ", matrixOfSynapticWeights[n][s] );
    }
    printf( "\n" );
  }
  printf( "Press enter ...\n" );
  getchar();
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = Debugprint: delay matrix
// =
// = Matrix of Synaptic Delays
// =
// =           0    1    2    3    4    5    6    7    8    9   10  ...  <- synapse
// = -------------------------------------------------------------- ...
// =  0 :    001  001  001  001  001  002  002  002  002  002  003  ...
// =  1 :    001  001  001  001  001  002  002  002  002  002  003  ...  <- conduction delay
// =  2 :    001  001  001  001  001  002  002  002  002  002  003  ...
// =  3 :    001  001  001  001  001  002  002  002  002  002  003  ...
// =  4 :    001  001  001  001  001  002  002  002  002  002  003  ...
// =  5 :    001  001  001  001  001  002  002  002  002  002  003  ...
// =  6 :    001  001  001  001  001  002  002  002  002  002  003  ...
// =  7 :    001  001  001  001  001  002  002  002  002  002  003  ...
// =  8 :    001  001  001  001  001  002  002  002  002  002  003  ...
// =  9 :    001  001  001  001  001  002  002  002  002  002  003  ...
// = 10 :    001  001  001  001  001  002  002  002  002  002  003  ...
// = ....
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void PrintMatrixOfSynapticDelays() {
  printf( "     Matrix of Synaptic Delays\n\n" );
  printf( "        " );
  for( int s = 0; s < NUM_SYNAPSES_PER_NEURON; ++s ) {
    printf( " %8d ", s );
  }
  printf( "\n" );
  printf( "-------" );
  for( int s = 0; s < NUM_SYNAPSES_PER_NEURON; ++s ) {
    printf( "----------" );
  }
  printf( " \n" );
  for( int n = 0; n < NUM_TOTAL_NEURONS; ++n ) {
    printf( "%3d :   ", n );
    for( int s = 0; s < NUM_SYNAPSES_PER_NEURON; ++s ) {
      printf( " %03d ", GetDelayOfConnection( n, s ));
    }
    printf( "\n" );
  }
  printf( "Press enter ...\n" );
  getchar();
}

#endif   // __RUN_REFACTORED_VERSION__
