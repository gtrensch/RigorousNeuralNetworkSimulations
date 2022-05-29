#ifndef __GLOBALS_H__
#define __GLOBALS_H__
/*
 *  globals.h
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

#ifdef __IMPORT_GLOBAL_VARIABLES__
  #define EXTERN extern
#else
  #define EXTERN
#endif

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   N E T W O R K   P A R A M E T E R S
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#define NUM_EXCITATORY_NEURONS      (int)(800)
#define NUM_INHIBITORY_NEURONS      (int)(200)
#define NUM_SYNAPSES_PER_NEURON     (int)(100)
#define MAX_SYNAPSE_DELAY           (int)(20)
#define MAX_SYNAPTIC_STRENGTH       (double)(10.0)

// down-scaled versions 20 neurons

// #define NUM_EXCITATORY_NEURONS      (int)(16)
// #define NUM_INHIBITORY_NEURONS      (int)(4)
// #define NUM_SYNAPSES_PER_NEURON     (int)(5)
// #define MAX_SYNAPSE_DELAY           (int)(5)
// #define MAX_SYNAPTIC_STRENGTH       (double)(10.0)

// down-scaled version 200 neurons

// #define NUM_EXCITATORY_NEURONS      (int)(160)
// #define NUM_INHIBITORY_NEURONS      (int)(40)
// #define NUM_SYNAPSES_PER_NEURON     (int)(20)
// #define MAX_SYNAPSE_DELAY           (int)(10)
// #define MAX_SYNAPTIC_STRENGTH       (double)(10.0)

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   M A C R O   D E F I N I T I O N S
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#define GET_RANDOM_INT(max)         (int(floor(max * ( 0.9999999 * double(rand()) / RAND_MAX))))
#define GET_RANDOM_FLOAT(max)       (max * ( 0.9999999 * double(rand()) / RAND_MAX))

#define NUM_TOTAL_NEURONS           (int)(NUM_EXCITATORY_NEURONS + NUM_INHIBITORY_NEURONS)
#define MAX_NUM_FIRINGS             (long int)(1000 * NUM_TOTAL_NEURONS)
#define INIT_EXC_SYNAPTIC_WEIGHT    (double)(6.0)
#define INIT_INH_SYNAPTIC_WEIGHT    (double)(-5.0)

#define RS_A                        (double)(0.02)
#define RS_D                        (double)(8.0)
#define RS_V_INIT                   (double)(-65.0)
#define RS_U_INIT                   (double)(0.2 * RS_V_INIT)

#define FS_A                        (double)(0.1)
#define FS_D                        (double)(2.0)
#define FS_V_INIT                   (double)(-65.0)
#define FS_U_INIT                   (double)(0.2 * FS_V_INIT)

#define RS_FS_B                     (double)(0.2)
#define RS_FS_C                     (double)(-65.0)

#define RS_FS_THR                   (double)(30.0)

#define I_EXT                       (double)(20.0)

#define I_OFFS_RS                   (double)(0.0)
#define I_OFFS_FS                   (double)(0.0)

#define RC_NOID                     (int)(-1)
#define RC_ERROR_EXIT               (int)(-2)

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   G L O B A L   D A T A   S T R U C T U R E S
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
EXTERN int      matrixOfPostSynapticNeurons        [ NUM_TOTAL_NEURONS ][ NUM_SYNAPSES_PER_NEURON ];
EXTERN double   matrixOfSynapticWeights            [ NUM_TOTAL_NEURONS ][ NUM_SYNAPSES_PER_NEURON ];
EXTERN double   matrixOfSynapticWeights_derivatives[ NUM_TOTAL_NEURONS ][ NUM_SYNAPSES_PER_NEURON ];

EXTERN int      numEntriesPerDelay[NUM_TOTAL_NEURONS][MAX_SYNAPSE_DELAY];
EXTERN int      listOfSynsByNeuronAndDelay[NUM_TOTAL_NEURONS][MAX_SYNAPSE_DELAY][NUM_SYNAPSES_PER_NEURON];

EXTERN int      numPreSynapticNeuronsOfTarget[ NUM_TOTAL_NEURONS ];

// REFACTOR COMMENT: assumption that it fits
EXTERN int      listOfPresynapticNeurons[ NUM_TOTAL_NEURONS ][ 3 * NUM_SYNAPSES_PER_NEURON ];
EXTERN int      listOfPresynapticDelays [ NUM_TOTAL_NEURONS ][ 3 * NUM_SYNAPSES_PER_NEURON ];

// REFACTOR COMMENT: dead code removed
// double*  listOfPointersToSynapticWeights  [NUM_TOTAL_NEURONS][3 * NUM_SYNAPSES_PER_NEURON];

EXTERN double*  listOfPointersToSynapticWeights_derivatives[NUM_TOTAL_NEURONS][3 * NUM_SYNAPSES_PER_NEURON];

EXTERN double   LTP[NUM_TOTAL_NEURONS][1001 + MAX_SYNAPSE_DELAY];
EXTERN double   LTD[NUM_TOTAL_NEURONS];

EXTERN long int  numFirings;
EXTERN long int  firings[MAX_NUM_FIRINGS][2];
#define IDX_TIME    (int)(0)
#define IDX_NEURON  (int)(1)

EXTERN double   I_ext[NUM_TOTAL_NEURONS];
EXTERN double   i_offs[NUM_TOTAL_NEURONS];

EXTERN double   a[NUM_TOTAL_NEURONS];                                                // neuronal dynamics parameters
EXTERN double   d[NUM_TOTAL_NEURONS];
EXTERN double   v[NUM_TOTAL_NEURONS];                                                // activity variables
EXTERN double   u[NUM_TOTAL_NEURONS];

#if( ODE_SOLVER_REFINEMENT )
EXTERN long int spike[NUM_TOTAL_NEURONS];                                            // spike count in an interval  (precise threshold detection)
#endif

#if( LOG_MIN_MAX_V_U )
EXTERN double   v_max;
EXTERN double   v_min;
EXTERN double   u_max;
EXTERN double   u_min;
#endif

EXTERN FILE*    pFileStimulusOutput;
EXTERN FILE*    pFileStimulusInput;

#endif   // __RUN_REFACTORED_VERSION__
#endif   // __GLOBALS_H__
