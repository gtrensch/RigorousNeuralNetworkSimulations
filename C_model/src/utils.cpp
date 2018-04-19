/*
 *  utils.cpp
 *
 *  This file is part of the refactored Izhikevich polychronization model application.
 *
 *  Copyright (C) 2018, Author: G. Trensch
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
  while( firings[idx][TIME] < 0 ) {
    idx++;
  }

  for( ; idx < numFirings; ++idx ) {
    fprintf( pFile, "%06d %03d %03d\n", simulationSecond, firings[idx][TIME], firings[idx][NEURON] );
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
