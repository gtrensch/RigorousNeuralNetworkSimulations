/*
 *  izhikevich.h
 *
 *  Copyright (C) 2016, G. Trensch, Forschungszentrum Jülich, JSC, Simulation & Data Laboratory Neuroscience
 *
 *  The Izhikevich console application is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  It is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this application. If not, see <http://www.gnu.org/licenses/>.
 *
 */

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   TYPE DEFINITIONS
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
typedef struct  _CNTXT      CNTXT;
typedef struct  _MAIN_ARGS  MAIN_ARGS;

typedef float (*PFUNCD)(float, float, void*);

typedef enum {
  STANDARD_EULER                    = 1,
  SYMPLECTIC_EULER                  = 2,
  EXPLICIT_EULER_SPINNAKER          = 3,
  EXPLICIT_EULER_NEST               = 4,
  EXPLICIT_IZHIKEVICH_NUMERICS_NEST = 5,
  GSL_LIBRARY                       = 6,
  HEUN_METHOD                       = 7,
  ADAPTIVE_EULER                    = 8
} ENUM_ODESOLVER_METHOD;

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   MACRO DEFINITIONS
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#define TRUE               1
#define FALSE              0
#define S_CNTXT            sizeof( CNTXT )
#define NULLCHR            '\0'
#define S_NAME             256

#define PLOT_COLOR_RED     1
#define PLOT_COLOR_YELLOW  2
#define PLOT_COLOR_GREEN   3
#define PLOT_COLOR_GREY    7
#define PLOT_COLOR_BROWN   8
#define PLOT_COLOR_BLUE    9  
#define PLOT_COLOR_CYAN    11
#define PLOT_COLOR_MAGENTA 13

#define PLOT_LINE_FULL     1
#define PLOT_LINE_DOTTED   2
#define PLOT_LINE_DASHED1  3
#define PLOT_LINE_DASHED2  5

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   COMMAND LINE PARAMETER DEFINITIONS AND DATA STRUCTURES
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// Option -c ... current shape
#define OPT_C  'c'

typedef enum {
  CURRENT_STEP_150      = 1,
  CURRENT_PULSE         = 2,
  CURRENT_NEGATIVE_STEP = 3
} ENUM_CURRENT_SHAPE;

#define STR_CURRENT_STEP_150       "step150"
#define STR_CURRENT_PULSE          "pulse"
#define STR_CURRENT_NEGATIVE_STEP  "nstep"

// Option -n ... neuron type
#define OPT_N  'n'

typedef enum {
  NEURON_TYPE_RS  = 1,       // regular spiking
  NEURON_TYPE_IB  = 2,       // intrinsically spiking
  NEURON_TYPE_CH  = 3,       // chattering
  NEURON_TYPE_FS  = 4,       // fast spiking
  NEURON_TYPE_TC  = 5,       // thalamo cortical 
  NEURON_TYPE_RZ  = 6,       // resonance
  NEURON_TYPE_LTS = 7        // low threshold spiking
} ENUM_NEURON_TYPE;

#define STR_NEURON_TYPE_RS   "RS"
#define STR_NEURON_TYPE_IB   "IB"
#define STR_NEURON_TYPE_CH   "CH"
#define STR_NEURON_TYPE_FS   "FS"
#define STR_NEURON_TYPE_TC   "TC"
#define STR_NEURON_TYPE_RZ   "RZ"
#define STR_NEURON_TYPE_LTS  "LTS"

// Option -f ... svg output file name
#define OPT_F  'f'

// Option -m ... ODE solver method
#define OPT_M  'm'

#define STR_ODE_SOLVER_METHOD_STANDARD_EULER         (char*)"stdEuler"
#define STR_ODE_SOLVER_METHOD_SYMPLECTIC_EULER       (char*)"sympEuler"
#define STR_ODE_SOLVER_METHOD_EULER_SPINNAKER        (char*)"SpiNNaker"  
#define STR_ODE_SOLVER_METHOD_EULER_NEST             (char*)"EulerNEST"
#define STR_ODE_SOLVER_METHOD_IZHIKEVICH_NEST        (char*)"IzhikevichNEST"
#define STR_ODE_SOLVER_METHOD_GSL                    (char*)"GSL"
#define STR_ODE_SOLVER_METHOD_HEUN                   (char*)"Heun"
#define STR_ODE_SOLVER_METHOD_HIGH_RESOLUTION_EULER  (char*)"highResEuler"
#define STR_ODE_SOLVER_METHOD_ADAPTIVE_EULER         (char*)"adaptiveEuler"

// Option -h ... help
#define OPT_H  'h'

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   DATA STRUCTURES
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
struct _MAIN_ARGS {
  ENUM_CURRENT_SHAPE  currentShape;
  ENUM_NEURON_TYPE    neuronType;
  char                outputSvgFilename[S_NAME];
  bool                fOdeSolverStandardEuler;
  bool                fOdeSolverSymplecticEuler;  
  bool                fOdeSolverEulerSpiNNaker;
  bool                fOdeSolverEulerNEST;
  bool                fOdeSolverIzhikevichNEST;
  bool                fOdeSolverGSL;
  bool                fOdeSolverHeun;
  bool                fOdeSolverHighResolutionEuler;
  bool                fOdeSolverAdaptiveEuler;
};

struct _CNTXT {
  MAIN_ARGS  args;
  char      neuronTypeName[S_NAME];
  float     neuronParam_A;
  float     neuronParam_B;
  float     neuronParam_C;
  float     neuronParam_D;
  float     neuronParam_THR;
  float     neuronParam_I_EXT;
  float     neuronParam_V_INIT;
  float     neuronParam_U_INIT;
  float     legendYOffset;
};

