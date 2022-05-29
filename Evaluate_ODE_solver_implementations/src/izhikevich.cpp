/*
 *  izhikevich.cpp
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

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
//  Comparison of different ODE solver implementations for solving the dynamics of the
//  Izhikevich neuron model [1].
//
//  References:
//
//  [1] Izhikevich, E. M. (2003). Simple model of spiking neurons. Trans. Neur. Netw. 14, 1569–1572.
//      doi:10.1109/TNN.2003.820440
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +


#include <plstream.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>
#include "izhikevich.h"
#include "izhikevich_model.h"
#include "ode_solvers.h"

#define PRINT_DEBUG  printf

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   SIMULATION PARAMETERS
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#define SIMULATION_TIME                    (float)(600.0)    // in milliseconds
#define INTEGRATION_STEP_SIZE              (float)(1.000)    // in milliseconds
#define HIGHRES_INTEGRATION_STEP_SIZE      (float)(0.001)    // in milliseconds
#define PLOT_SHIFT                         (float)(0.000)    // shift plot for debug, e.g., 5.0


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   FORWARD DECLARATIONS
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void drawGraph( CNTXT* pCntxt, plstream* pPlStream, int numDataPoints, PLFLT* pTimesArray, PLFLT* pVoltagesArray, int plot_color, int plot_line_style, char* pLegendText );
void exitProgram( CNTXT* pCntxt );
void finalizeContext( CNTXT* pCntxt );
void finalizePlot( CNTXT* pCntxt, plstream* pPlStream );
void finalizeSimulation( CNTXT* pCntxt, PLFLT* pTimesArray, PLFLT* pVoltagesArray, float* pCurrentArray );
void initializeContext( CNTXT** ppCntxt );
void initializePlot( CNTXT* pCntxt, plstream** pPlStream );
void initializeSimulation( CNTXT* pCntxt, float simulationTime, float integrationStepSize, PLFLT** ppTimesArray, PLFLT** ppVoltagesArray, PLFLT** ppFrequencyArray, float** ppCurrentArray );
int  runSimulation( CNTXT* pCntxt, ENUM_ODESOLVER_METHOD method, float simulationTime, float integrationStepSize, PLFLT* pTimesArray, PLFLT* pVoltagesArray, PLFLT* pFrequencyArray, float* pCurrentArray );
void setCurrentShape( CNTXT* pCntxt, float* pCurrentArray, int numDataPoints, float integrationStepSize );

void processMainArguments( CNTXT* pCntxt, int argc, char** argv );
void processOPT_C( CNTXT* pCntxt, int argc, char** argv );
void processOPT_F( CNTXT* pCntxt, int argc, char** argv );
void processOPT_H( CNTXT* pCntxt, int argc, char** argv );
void processOPT_M( CNTXT* pCntxt, int argc, char** argv );
void processOPT_N( CNTXT* pCntxt, int argc, char** argv );


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   MAIN ENTRY
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
int main (int argc, char** argv)
{
  CNTXT*     pCntxt          = NULL;
  PLFLT*     pVoltagesArray  = NULL;
  PLFLT*     pTimesArray     = NULL;
  PLFLT*     pFrequencyArray = NULL;
  float*     pCurrentArray   = NULL;
  int        numDataPoints   = 0;
  plstream*  pPlStream       = NULL;

  initializeContext( &pCntxt );
  processMainArguments( pCntxt, argc, argv );
  initializePlot( pCntxt, &pPlStream );
  
  // Explicit Forward Euler
  if( pCntxt->args.fOdeSolverStandardEuler ) { 
    PRINT_DEBUG( "    Debug: + + + perform standard Euler + + +\n" );
    initializeSimulation( pCntxt, SIMULATION_TIME, INTEGRATION_STEP_SIZE, &pTimesArray, &pVoltagesArray, &pFrequencyArray, &pCurrentArray );
    numDataPoints = runSimulation( pCntxt, STANDARD_EULER, SIMULATION_TIME, INTEGRATION_STEP_SIZE, pTimesArray, pVoltagesArray, pFrequencyArray, pCurrentArray );
    drawGraph( pCntxt, pPlStream, numDataPoints, pTimesArray, pVoltagesArray, PLOT_COLOR_RED, PLOT_LINE_FULL, STR_ODE_SOLVER_METHOD_STANDARD_EULER );
    finalizeSimulation( pCntxt, pTimesArray, pVoltagesArray, pCurrentArray );
  }

  // Symplectic Forward Euler
  if( pCntxt->args.fOdeSolverSymplecticEuler ) {
    PRINT_DEBUG( "    Debug: + + + symplectic Euler + + +\n" );
    initializeSimulation( pCntxt, SIMULATION_TIME, INTEGRATION_STEP_SIZE, &pTimesArray, &pVoltagesArray, &pFrequencyArray, &pCurrentArray );
    numDataPoints = runSimulation( pCntxt, SYMPLECTIC_EULER, SIMULATION_TIME, INTEGRATION_STEP_SIZE, pTimesArray, pVoltagesArray, pFrequencyArray, pCurrentArray );
    drawGraph( pCntxt, pPlStream, numDataPoints, pTimesArray, pVoltagesArray, PLOT_COLOR_GREEN, PLOT_LINE_FULL, STR_ODE_SOLVER_METHOD_SYMPLECTIC_EULER );
    finalizeSimulation( pCntxt, pTimesArray, pVoltagesArray, pCurrentArray );
  }    

  // SpiNNaker implementation
  if( pCntxt->args.fOdeSolverEulerSpiNNaker ) {
    PRINT_DEBUG( "    Debug: + + + perform explicit Euler - SpiNNaker implementation + + +\n" );
    initializeSimulation( pCntxt, SIMULATION_TIME, INTEGRATION_STEP_SIZE, &pTimesArray, &pVoltagesArray, &pFrequencyArray, &pCurrentArray );
    numDataPoints = runSimulation( pCntxt, EXPLICIT_EULER_SPINNAKER, SIMULATION_TIME, INTEGRATION_STEP_SIZE, pTimesArray, pVoltagesArray, pFrequencyArray, pCurrentArray );
    drawGraph( pCntxt, pPlStream, numDataPoints, pTimesArray, pVoltagesArray, PLOT_COLOR_MAGENTA, PLOT_LINE_FULL, STR_ODE_SOLVER_METHOD_EULER_SPINNAKER );
    finalizeSimulation( pCntxt, pTimesArray, pVoltagesArray, pCurrentArray );
  }

  // NEST implementation
  if( pCntxt->args.fOdeSolverEulerNEST ) {
    PRINT_DEBUG( "    Debug: + + + explicit Euler - NEST implementation + + +\n" );  
    initializeSimulation( pCntxt, SIMULATION_TIME, INTEGRATION_STEP_SIZE, &pTimesArray, &pVoltagesArray, &pFrequencyArray, &pCurrentArray );
    numDataPoints = runSimulation( pCntxt, EXPLICIT_EULER_NEST, SIMULATION_TIME, INTEGRATION_STEP_SIZE, pTimesArray, pVoltagesArray, pFrequencyArray, pCurrentArray );
    drawGraph( pCntxt, pPlStream, numDataPoints, pTimesArray, pVoltagesArray, PLOT_COLOR_YELLOW, PLOT_LINE_FULL, STR_ODE_SOLVER_METHOD_EULER_NEST );
    finalizeSimulation( pCntxt, pTimesArray, pVoltagesArray, pCurrentArray );
  }

  // Original Izhikevich implementation
  if( pCntxt->args.fOdeSolverIzhikevichNEST ) {
    PRINT_DEBUG( "    Debug: + + + explicit Izhikevich numerics - NEST implementation + + +\n" );
    initializeSimulation( pCntxt, SIMULATION_TIME, INTEGRATION_STEP_SIZE, &pTimesArray, &pVoltagesArray, &pFrequencyArray, &pCurrentArray );
    numDataPoints = runSimulation( pCntxt, EXPLICIT_IZHIKEVICH_NUMERICS_NEST, SIMULATION_TIME, INTEGRATION_STEP_SIZE, pTimesArray, pVoltagesArray, pFrequencyArray, pCurrentArray );
    drawGraph( pCntxt, pPlStream, numDataPoints, pTimesArray, pVoltagesArray, PLOT_COLOR_CYAN, PLOT_LINE_FULL, STR_ODE_SOLVER_METHOD_IZHIKEVICH_NEST );
    finalizeSimulation( pCntxt, pTimesArray, pVoltagesArray, pCurrentArray );
  }
 
  // GSL library, rkf45
  if( pCntxt->args.fOdeSolverGSL ) {
    PRINT_DEBUG( "    Debug: + + + GSL library + + +\n" );
    initializeSimulation( pCntxt, SIMULATION_TIME, INTEGRATION_STEP_SIZE, &pTimesArray, &pVoltagesArray, &pFrequencyArray, &pCurrentArray );
    numDataPoints = runSimulation( pCntxt, GSL_LIBRARY, SIMULATION_TIME, INTEGRATION_STEP_SIZE, pTimesArray, pVoltagesArray, pFrequencyArray, pCurrentArray );
    drawGraph( pCntxt, pPlStream, numDataPoints, pTimesArray, pVoltagesArray, PLOT_COLOR_GREY, PLOT_LINE_FULL, STR_ODE_SOLVER_METHOD_GSL );
    finalizeSimulation( pCntxt, pTimesArray, pVoltagesArray, pCurrentArray );
  }

  // Heun
  if( pCntxt->args.fOdeSolverHeun ) {
    PRINT_DEBUG( "    Debug: + + + Heun's method + + +\n" );
    initializeSimulation( pCntxt, SIMULATION_TIME, INTEGRATION_STEP_SIZE, &pTimesArray, &pVoltagesArray, &pFrequencyArray, &pCurrentArray );
    numDataPoints = runSimulation( pCntxt, HEUN_METHOD, SIMULATION_TIME, INTEGRATION_STEP_SIZE, pTimesArray, pVoltagesArray, pFrequencyArray, pCurrentArray );
    drawGraph( pCntxt, pPlStream, numDataPoints, pTimesArray, pVoltagesArray, PLOT_COLOR_BROWN, PLOT_LINE_FULL, STR_ODE_SOLVER_METHOD_HEUN );
    finalizeSimulation( pCntxt, pTimesArray, pVoltagesArray, pCurrentArray );
  }

  // Adaptive Euler
  if( pCntxt->args.fOdeSolverAdaptiveEuler ) { 
    PRINT_DEBUG( "    Debug: + + + perform adaptive Euler + + +\n" );
    initializeSimulation( pCntxt, SIMULATION_TIME, INTEGRATION_STEP_SIZE, &pTimesArray, &pVoltagesArray, &pFrequencyArray, &pCurrentArray );
    numDataPoints = runSimulation( pCntxt, ADAPTIVE_EULER, SIMULATION_TIME, INTEGRATION_STEP_SIZE, pTimesArray, pVoltagesArray, pFrequencyArray, pCurrentArray );
    drawGraph( pCntxt, pPlStream, numDataPoints, pTimesArray, pVoltagesArray, PLOT_COLOR_BROWN, PLOT_LINE_FULL, STR_ODE_SOLVER_METHOD_ADAPTIVE_EULER );
    finalizeSimulation( pCntxt, pTimesArray, pVoltagesArray, pCurrentArray );
  }

  // High resolution Euler
  if( pCntxt->args.fOdeSolverHighResolutionEuler ) {
    PRINT_DEBUG( "    Debug: + + + reference high resolution standard Euler + + +\n" );
    initializeSimulation( pCntxt, SIMULATION_TIME, HIGHRES_INTEGRATION_STEP_SIZE, &pTimesArray, &pVoltagesArray, &pFrequencyArray, &pCurrentArray );
    numDataPoints = runSimulation( pCntxt, STANDARD_EULER, SIMULATION_TIME, HIGHRES_INTEGRATION_STEP_SIZE, pTimesArray, pVoltagesArray, pFrequencyArray, pCurrentArray );
    drawGraph( pCntxt, pPlStream, numDataPoints, pTimesArray, pVoltagesArray, PLOT_COLOR_BLUE, PLOT_LINE_DOTTED, STR_ODE_SOLVER_METHOD_HIGH_RESOLUTION_EULER );
    finalizeSimulation( pCntxt, pTimesArray, pVoltagesArray, pCurrentArray );
  }

  finalizePlot( pCntxt, pPlStream ); 
  finalizeContext( pCntxt ); 
  printf( "    Terminated normally\n" );
  return 0;
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
//   HELPER FUNCTIONS
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
//   Plot graph, write to svg file.
//   Multiple graphs can be plotted and overlaid.
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void drawGraph( CNTXT* pCntxt, plstream* pPlStream, int numDataPoints, PLFLT* pTimesArray, PLFLT* pVoltagesArray, int plot_color, int plot_line_style, char* pLegendText )
{
  PRINT_DEBUG( "    Debug: plot graph, i.e. write to svg file: %s\n", pCntxt->args.outputSvgFilename );
  pllsty( plot_line_style );
  plcol0( plot_color );
  pPlStream->line( numDataPoints, pTimesArray, pVoltagesArray );
  plptex( 80.0, pCntxt->legendYOffset, 0.0, 0.0, 0.5, pLegendText );
  pCntxt->legendYOffset += 5;
  return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
//   Free program context area and exit
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void exitProgram( CNTXT* pCntxt )
{
  finalizeContext( pCntxt );
  exit( 0 );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
//   Free program context area.
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void finalizeContext( CNTXT* pCntxt )
{
  if( pCntxt != NULL ) {
    free( pCntxt );
  }
  return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
//   Delete Plplot stream
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void finalizePlot( CNTXT* pCntxt, plstream* pPlStream )
{
  delete pPlStream;
  return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
//   Free data arrays
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void finalizeSimulation( CNTXT* pCntxt, PLFLT* pTimesArray, PLFLT* pVoltagesArray, float* pCurrentArray )
{
  PRINT_DEBUG( "    Debug: finalize\n" );
  free( pTimesArray );
  free( pVoltagesArray );
  free( pCurrentArray );
  return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
//   Allocate and initialize the program's context area
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void initializeContext( CNTXT** ppCntxt )
{
  PRINT_DEBUG( "    Debug: initialize context\n" );
  *ppCntxt = (CNTXT*)malloc( S_CNTXT );
  if( *ppCntxt == NULL ) {
    printf( "ERRROR: Failed to allocate program context area\n");
    exitProgram( NULL );
  }
  memset( *ppCntxt, NULLCHR, S_CNTXT );
  return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
//   Set up Plplot.
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void initializePlot( CNTXT* pCntxt, plstream** pPlStream )
{
  *pPlStream = new plstream();
  plsdev("svg");
  plsfnam( pCntxt->args.outputSvgFilename );
  plscolbg( 255, 255, 255 );
  plscmap0n(16);
  (*pPlStream)->init();

  (*pPlStream)->env( 0, SIMULATION_TIME, -80, 50, 0, 0 );
  //                 |  |                 |   |   |  |
  //                 |  |                 |   |   |  +---  axis - 0 draw axis box and numeric labels
  //                 |  |                 |   |   +------  just - 0 axis scales independently
  //                 |  |                 |   +----------  ymax
  //                 |  |                 +--------------  ymin
  //                 |  +--------------------------------  xmax
  //                 +-----------------------------------  xmin

  (*pPlStream)->lab( "Time [ms]", "Voltage [mv]", pCntxt->neuronTypeName );

  return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
//   Allocate and set up data structures
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
void initializeSimulation( CNTXT* pCntxt, float simulationTime, float integrationStepSize, PLFLT** ppTimesArray, PLFLT** ppVoltagesArray, PLFLT** ppFrequencyArray, float** ppCurrentArray )
{

  PRINT_DEBUG( "    Debug: initialize simulation and prepare data arrays\n" );

  PLFLT*  pVoltagesArray  = NULL;
  PLFLT*  pTimesArray     = NULL;
  PLFLT*  pFrequencyArray = NULL;
  float*  pCurrentArray   = NULL;
  int     numDataPoints   = floor(simulationTime / integrationStepSize) + 10;   // due to rounding errors there might be a few more datapoints returned from the simulation run; just add 10 to be on the save side
  
  PRINT_DEBUG( "    Debug: number of data points calculated: %d\n", numDataPoints );

  pVoltagesArray = (PLFLT*)malloc( numDataPoints * sizeof(PLFLT) );
  pTimesArray = (PLFLT*)malloc( numDataPoints * sizeof(PLFLT) );
  pFrequencyArray = (PLFLT*)malloc( numDataPoints * sizeof(PLFLT) );
  pCurrentArray = (float*)malloc( numDataPoints * sizeof(float) );

  if( NULL == pVoltagesArray || NULL == pTimesArray || NULL == pFrequencyArray || NULL == pCurrentArray ) {
    printf( "ERROR: failed to allocate memory\n" );
    exitProgram( pCntxt );
  }
  memset( pVoltagesArray, 0x00, numDataPoints * sizeof(PLFLT) );
  memset( pTimesArray, 0x00, numDataPoints * sizeof(PLFLT) );
  memset( pFrequencyArray, 0x00, numDataPoints * sizeof(PLFLT) );
  memset( pCurrentArray, 0x00, numDataPoints * sizeof(float) );

  setCurrentShape( pCntxt, pCurrentArray, numDataPoints, integrationStepSize );

  *ppTimesArray = pTimesArray;
  *ppVoltagesArray = pVoltagesArray;
  *ppFrequencyArray = pFrequencyArray;
  *ppCurrentArray = pCurrentArray;

  // set Izhikevich model constants
  IZHIKEVICH_A = pCntxt->neuronParam_A;
  IZHIKEVICH_B = pCntxt->neuronParam_B;
  IZHIKEVICH_C = pCntxt->neuronParam_C;
  IZHIKEVICH_D = pCntxt->neuronParam_D;
  IZHIKEVICH_THR = pCntxt->neuronParam_THR;        

  return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
//   Depending on parameter "method", call the corresponding ODE solver and
//   perform the simulation.
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
int runSimulation( CNTXT* pCntxt, ENUM_ODESOLVER_METHOD method, float simulationTime, float integrationStepSize
                 , PLFLT* pTimesArray, PLFLT* pVoltagesArray, PLFLT* pFrequencyArray, float* pCurrentArray )
{

  PFUNCD pODE1dxdt = (PFUNCD)dvdt;      // Izhikevich ODE1 
  PFUNCD pODE2dydt = (PFUNCD)dudt;      // Izhikevich ODE2 
  float  v_t = pCntxt->neuronParam_V_INIT;
  float  u_t = pCntxt->neuronParam_U_INIT;
  float  v_tplus1 = 0;
  float  u_tplus1 = 0;
  float  t = 0;
  int    dataPointIndex = 0;

  PRINT_DEBUG( "    Debug: run simulation\n" );

  switch( method ) {
    case STANDARD_EULER:
      t = PLOT_SHIFT * 1;
      while( t <= simulationTime ) {
        pTimesArray[dataPointIndex] = t;     
        pVoltagesArray[dataPointIndex] = v_t;
        ode2dSolver_StandardEuler( pODE1dxdt, pODE2dydt, integrationStepSize, v_t, u_t, &v_tplus1, &u_tplus1, &pCurrentArray[dataPointIndex] ); 
        if( threshold( &v_tplus1, &u_tplus1 )) {
          pFrequencyArray[dataPointIndex] = 30;
        };
        v_t = v_tplus1;
        u_t = u_tplus1;
        t += integrationStepSize;
        dataPointIndex++;
      }
      break;

    case SYMPLECTIC_EULER:
      t = PLOT_SHIFT * 2;
      while( t <= simulationTime ) {
        pTimesArray[dataPointIndex] = t;     
        pVoltagesArray[dataPointIndex] = v_t;
        ode2dSolver_SymplecticEuler( pODE1dxdt, pODE2dydt, integrationStepSize, v_t, u_t, &v_tplus1, &u_tplus1, &pCurrentArray[dataPointIndex] );
        threshold( &v_tplus1, &u_tplus1 );
        v_t = v_tplus1;
        u_t = u_tplus1;
        t += integrationStepSize;
        dataPointIndex++;
      }
      break;

    case EXPLICIT_EULER_SPINNAKER:
      t = PLOT_SHIFT * 3;
      while( t <= simulationTime ) {
        pTimesArray[dataPointIndex] = t;     
        pVoltagesArray[dataPointIndex] = v_t;
        explizitSpiNNaker( integrationStepSize, &v_t, &u_t, pCurrentArray[dataPointIndex] );
        t += integrationStepSize;
        dataPointIndex++;
      }
      break;

    case EXPLICIT_EULER_NEST:
      t = PLOT_SHIFT * 4;
      while( t <= simulationTime ) {
        pTimesArray[dataPointIndex] = t;     
        pVoltagesArray[dataPointIndex] = v_t;
        explizitNEST1( integrationStepSize, &v_t, &u_t, pCurrentArray[dataPointIndex] );
        t += integrationStepSize;
        dataPointIndex++;
      }    
      break;

    case EXPLICIT_IZHIKEVICH_NUMERICS_NEST:
      t = PLOT_SHIFT * 5;
      while( t <= simulationTime ) {
        pTimesArray[dataPointIndex] = t;     
        pVoltagesArray[dataPointIndex] = v_t;
        explizitNEST2( integrationStepSize, &v_t, &u_t, pCurrentArray[dataPointIndex] );
        t += integrationStepSize;
        dataPointIndex++;
      }
      break;

    case GSL_LIBRARY:
      t = PLOT_SHIFT * 6;
      gslODESolver_Init();
      while( t <= simulationTime ) {
        pTimesArray[dataPointIndex] = t;     
        pVoltagesArray[dataPointIndex] = v_t;
        gslODESolver_Integrate( integrationStepSize, &v_t, &u_t, pCurrentArray[dataPointIndex] );
        t += integrationStepSize;
        dataPointIndex++;
      }
      break;

    case HEUN_METHOD:
      t = PLOT_SHIFT * 7;
      while( t <= simulationTime ) {
        pTimesArray[dataPointIndex] = t;     
        pVoltagesArray[dataPointIndex] = v_t;
        ode2dSolver_Heun( pODE1dxdt, pODE2dydt, integrationStepSize, v_t, u_t, &v_tplus1, &u_tplus1, &pCurrentArray[dataPointIndex] ); 
        threshold( &v_tplus1, &u_tplus1 );
        v_t = v_tplus1;
        u_t = u_tplus1;
        t += integrationStepSize;
        dataPointIndex++;
      }
      break;

    case ADAPTIVE_EULER:
      t = PLOT_SHIFT * 8;
      while( t <= simulationTime ) {
        pTimesArray[dataPointIndex] = t;     
        pVoltagesArray[dataPointIndex] = v_t;
        ode2dSolver_AdaptiveEuler( pODE1dxdt, pODE2dydt, integrationStepSize, v_t, u_t, &v_tplus1, &u_tplus1, &pCurrentArray[dataPointIndex] ); 
        threshold( &v_tplus1, &u_tplus1 );
        v_t = v_tplus1;
        u_t = u_tplus1;
        t += integrationStepSize;
        dataPointIndex++;
      }
      break;


    default:
      printf( "ERROR: ODE solver method not defined\n");
      exitProgram( pCntxt );
  }
  
  PRINT_DEBUG( "    Debug: number of data points: %d\n", dataPointIndex );
  return( dataPointIndex );  // return number of data points, hence the dataPointIndex
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
//   Fill the current-data-array with the desired current-shape.
//   The size of the array, i.e. the number od datapoints depend on the 
//   resolution, hence the integration step size.
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void setCurrentShape( CNTXT* pCntxt, float* pCurrentArray, int numDataPoints, float integrationStepSize )
{

  PRINT_DEBUG( "    Debug: initialize current-shape array\n" );

  switch( pCntxt->args.currentShape ) {

    case CURRENT_STEP_150: 
      for( int i = 0; i < numDataPoints; ++i ) {
        pCurrentArray[i] = pCntxt->neuronParam_I_EXT;
      }
      for( int i = 0; i < floor( 150 / integrationStepSize ); ++i ) {
        pCurrentArray[i] = 0;
      }
      // printf( "             i     +------------------------      \n" );
      // printf( "             0 ----+                              \n" );
      // printf( "                  150                    t [ms]   \n" );  
      break;

    case CURRENT_PULSE: 
      for( int i = 0; i < numDataPoints; ++i ) {
        pCurrentArray[i] = pCntxt->neuronParam_I_EXT;
      }
      for( int i = 0; i < floor( 150 / integrationStepSize ); ++i ) {
        pCurrentArray[i] = 0;
      }
      for( int i = floor( 200 / integrationStepSize ) ; i < floor( (200 + 2) / integrationStepSize ); ++i ) {
        pCurrentArray[i] = pCntxt->neuronParam_I_EXT + 5 ;
      }
      // printf( "                              +--+                \n" );
      // printf( "                              |  | i + 5 pA       \n" );
      // printf( "             i     +----------+  +----------      \n" );
      // printf( "             0 ----+                              \n" );
      // printf( "                  150        200 202     t [ms]   \n" );  
      break;

    case CURRENT_NEGATIVE_STEP: 
      for( int i = 0; i < numDataPoints; ++i ) {
        pCurrentArray[i] = - pCntxt->neuronParam_I_EXT;
      }
      for( int i = floor( 150 / integrationStepSize ); i < numDataPoints; ++i ) {
        pCurrentArray[i] = 0;
      }
      // printf( "             0     +------------------------      \n" );
      // printf( "            -i ----+                              \n" );
      // printf( "                  150                    t [ms]   \n" );  
      break;

    default:
      printf( "ERROR: Failed to set up current shape\n" );
      exitProgram( pCntxt ) ;
  }  

  return;
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
//   HELPER FUNCTIONS TO PROCESS PROGRAM MAIN ARGUMENTS
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
//   Check program options
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
void processMainArguments( CNTXT* pCntxt, int argc, char** argv )
{
  while( --argc ) {

    if( argv[argc][0] != '-' ) {
      printf( "ERROR: Parameter error\n");
      processOPT_H( pCntxt, argc, argv );
      exitProgram( pCntxt );
    }

    switch( argv[argc][1] ) {

      case OPT_C:
        processOPT_C( pCntxt, argc, argv );
        break;

      case OPT_F:
        processOPT_F( pCntxt, argc, argv );
        break;

      case OPT_H:
        processOPT_H( pCntxt, argc, argv );
        exitProgram( pCntxt );
        break;

      case OPT_M:
        processOPT_M( pCntxt, argc, argv );
        break;

      case OPT_N:
        processOPT_N( pCntxt, argc, argv );
        break;

      default:
        printf( "ERROR: Parameter error\n");
        processOPT_H( pCntxt, argc, argv );
        exitProgram( pCntxt );
    }
  }

  // parameter checks

  if( pCntxt->args.currentShape == 0 ) {
    printf( "ERROR: option -c not set\n" );
    processOPT_H( pCntxt, argc, argv );
    exitProgram( pCntxt );
  }
  if( pCntxt->args.outputSvgFilename[0] == NULLCHR ) {
    printf( "ERROR: option -f not set\n" );
    processOPT_H( pCntxt, argc, argv );
    exitProgram( pCntxt );
  }
  if( pCntxt->neuronParam_A == 0 ) {
    printf( "ERROR: option -n not set\n" );
    processOPT_H( pCntxt, argc, argv );
    exitProgram( pCntxt );
  }
  if(( pCntxt->args.fOdeSolverStandardEuler || pCntxt->args.fOdeSolverSymplecticEuler ||  pCntxt->args.fOdeSolverEulerSpiNNaker || 
       pCntxt->args.fOdeSolverEulerNEST || pCntxt->args.fOdeSolverIzhikevichNEST ||  pCntxt->args.fOdeSolverHighResolutionEuler ||
       pCntxt->args.fOdeSolverGSL || pCntxt->args.fOdeSolverHeun || pCntxt->args.fOdeSolverAdaptiveEuler ) == 0 ) {
    printf( "ERROR: option -m not set\n" );
    processOPT_H( pCntxt, argc, argv );
    exitProgram( pCntxt );
  }

  printf( "+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +\n" );
  printf( "+           INZHIKEVICH MODEL SIMULATION                      +\n" );
  printf( "+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +\n" );
  printf( "Following parameters are in use for this run:\n" );
  printf( "    SVG Output File Name: %s\n", pCntxt->args.outputSvgFilename );
  printf( "    Neuron Type         : %s\n", pCntxt->neuronTypeName );
  
  switch( pCntxt->args.currentShape ) {

    case CURRENT_STEP_150: 
      printf( "    Current Shape       : step current at 150ms\n" );
      printf( "                          i     +------------------------      \n" );
      printf( "                          0 ----+                              \n" );
      printf( "                               150                    t [ms]   \n" );  
      break;

    case CURRENT_PULSE: 
      printf( "    Current Shape       : step current at 150ms with pulse at 200ms for 2ms\n" );
      printf( "                                           +--+                \n" );
      printf( "                                           |  | i + 5 pA       \n" );
      printf( "                          i     +----------+  +----------      \n" );
      printf( "                          0 ----+                              \n" );
      printf( "                               150        200 202     t [ms]   \n" );  
      break;

    case CURRENT_NEGATIVE_STEP: 
      printf( "    Current Shape       : negative step current at 150ms\n" );
      printf( "                           0     +------------------------      \n" );
      printf( "                          -i ----+                              \n" );
      printf( "                                150                    t [ms]   \n" );  
      break;
  }

  printf( "    ODE solver selected : \n" );
  if( pCntxt->args.fOdeSolverStandardEuler ) {
      printf( "                          Standard Euler\n" );
  }
  if( pCntxt->args.fOdeSolverSymplecticEuler ) {
      printf( "                          Symplextic Euler\n" );
  }
  if( pCntxt->args.fOdeSolverEulerSpiNNaker ) {
      printf( "                          Explicit Euler SpiNNaker implementation\n" );
  }
  if( pCntxt->args.fOdeSolverEulerNEST ) {
      printf( "                          Explicit Euler NEST implementation\n" );
  }
  if( pCntxt->args.fOdeSolverIzhikevichNEST ) {
      printf( "                          Izhikevich numerics NEST implementation\n" );
  }
  if( pCntxt->args.fOdeSolverGSL ) {
      printf( "                          GSL library\n" );
  }
  if( pCntxt->args.fOdeSolverHeun ) {
      printf( "                          Heun's method\n" );
  }
  if( pCntxt->args.fOdeSolverHighResolutionEuler ) {
      printf( "                          Standard Euler with a high integration step resolution, i.e. 1 micro second\n" );
  }
  if( pCntxt->args.fOdeSolverAdaptiveEuler ) {
      printf( "                          Standard Euler with an adaptive integration step size\n" );
  }


  return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
//   Set the current shape
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
void processOPT_C( CNTXT* pCntxt, int argc, char** argv ) 
{
  if( strncmp( &argv[argc][2], STR_CURRENT_STEP_150, sizeof( STR_CURRENT_STEP_150 )) == 0 ) {
    pCntxt->args.currentShape = CURRENT_STEP_150;
  }
  else if ( strncmp( &argv[argc][2], STR_CURRENT_PULSE, sizeof( STR_CURRENT_PULSE )) == 0 ) {
    pCntxt->args.currentShape = CURRENT_PULSE;
  }
  else if ( strncmp( &argv[argc][2], STR_CURRENT_NEGATIVE_STEP, sizeof( STR_CURRENT_NEGATIVE_STEP )) == 0 ) {
    pCntxt->args.currentShape = CURRENT_NEGATIVE_STEP;
  }
  else {
    printf( "ERROR: Parameter error, failed to set current shape\n" );
    exitProgram( pCntxt );
  }

  return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
//   Set svg output file name
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
void processOPT_F( CNTXT* pCntxt, int argc, char** argv )
{
  strncpy( pCntxt->args.outputSvgFilename, &argv[argc][2], S_NAME - 1 );
  pCntxt->args.outputSvgFilename[S_NAME - 1] = NULLCHR;
  return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
//   Show help text
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void processOPT_H( CNTXT* pCntxt, int argc, char** argv )
{
  printf( "USAGE  : izhikevich -n<neuron type> -c<current shape> -f<svg output path/file.svg> -m<ODE solver method> ... -m<ODE solver method>\n" );
  printf( "         -f     filename         svg output file name, e.g. svg\\RS_step150.svg\n" );
  printf( "         -n     RS               regular spiking\n" );
  printf( "                IB               intrinsically bursting\n" );
  printf( "                CH               chattering\n" );
  printf( "                FS               fast spiking\n" );  
  printf( "                TC               thalamo cortical\n" );
  printf( "                RZ               resonance\n" );
  printf( "                LTS              low threshold spiking\n" );
  printf( "\n" );      
  printf( "         -c     step150          step current at 150ms\n" );
  printf( "                                 i     +------------------------      \n" );
  printf( "                                 0 ----+                              \n" );
  printf( "                                      150                    t [ms]   \n" );  
  printf( "                pulse            step current at 150ms with pulse at 200ms for 2ms\n" );
  printf( "                                                  +--+                \n" );
  printf( "                                                  |  | i + 5 pA       \n" );
  printf( "                                 i     +----------+  +----------      \n" );
  printf( "                                 0 ----+                              \n" );
  printf( "                                      150        200 202     t [ms]   \n" );    
  printf( "                nstep            negative step current at 150ms\n" );
  printf( "                                  0     +------------------------      \n" );
  printf( "                                 -i ----+                              \n" );
  printf( "                                      150                    t [ms]   \n" );  
  printf( "\n" );
  printf( "         -m     stdEuler         standard Euler\n" );
  printf( "                sympEuler        symplectic Euler\n" );
  printf( "                SpiNNaker        explicit Euler SpiNNaker implementation\n" );
  printf( "                EulerNEST        explicit Euler NEST implementation\n" );  
  printf( "                IzhikevichNEST   Izhikevich numerics implementation in NEST\n" );  
  printf( "                GSL              use GSL library\n" );  
  printf( "                Heun             Heun method\n" );    
  printf( "                highResEuler     standard Euler with a high resolution integration step size, i.e. 1 micro second\n" );  
  printf( "                adaptiveEuler    standard Euler with an adaptive integration step size\n" );  
  printf( "\n" );
  printf( "EXAMPLE: compare different ODE solver methods for a regular spiking neuron with a positive external step current at 150 milliseconds\n" );
  printf( "         izhikevich -nRS -cstep150 -fsvg/RS_step150.svg -mstdEuler -msympEuler -mSpiNNaker -mEulerNEST -mIzhikevichNEST -mhighResEuler\n" );
  printf( "\n" );

  return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
//   Set ODE solver method(s)
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
void processOPT_M( CNTXT* pCntxt, int argc, char** argv )
{
  if( strncmp( &argv[argc][2], STR_ODE_SOLVER_METHOD_STANDARD_EULER, sizeof( STR_ODE_SOLVER_METHOD_STANDARD_EULER )) == 0 ) {
    pCntxt->args.fOdeSolverStandardEuler = TRUE;
  }
  else if( strncmp( &argv[argc][2], STR_ODE_SOLVER_METHOD_SYMPLECTIC_EULER, sizeof( STR_ODE_SOLVER_METHOD_SYMPLECTIC_EULER )) == 0 ) {
    pCntxt->args.fOdeSolverSymplecticEuler = TRUE;
  }
  else if( strncmp( &argv[argc][2], STR_ODE_SOLVER_METHOD_EULER_SPINNAKER, sizeof( STR_ODE_SOLVER_METHOD_EULER_SPINNAKER )) == 0 ) {
    pCntxt->args.fOdeSolverEulerSpiNNaker = TRUE;
  }
  else if( strncmp( &argv[argc][2], STR_ODE_SOLVER_METHOD_EULER_NEST, sizeof( STR_ODE_SOLVER_METHOD_EULER_NEST )) == 0 ) {
    pCntxt->args.fOdeSolverEulerNEST = TRUE;
  }
  else if( strncmp( &argv[argc][2], STR_ODE_SOLVER_METHOD_IZHIKEVICH_NEST, sizeof( STR_ODE_SOLVER_METHOD_IZHIKEVICH_NEST )) == 0 ) {
    pCntxt->args.fOdeSolverIzhikevichNEST = TRUE;
  }
  else if( strncmp( &argv[argc][2], STR_ODE_SOLVER_METHOD_GSL, sizeof( STR_ODE_SOLVER_METHOD_GSL )) == 0 ) {
    pCntxt->args.fOdeSolverGSL = TRUE;
  }
  else if( strncmp( &argv[argc][2], STR_ODE_SOLVER_METHOD_HEUN, sizeof( STR_ODE_SOLVER_METHOD_HEUN )) == 0 ) {
    pCntxt->args.fOdeSolverHeun = TRUE;
  }
  else if( strncmp( &argv[argc][2], STR_ODE_SOLVER_METHOD_HIGH_RESOLUTION_EULER, sizeof( STR_ODE_SOLVER_METHOD_HIGH_RESOLUTION_EULER )) == 0 ) {
    pCntxt->args.fOdeSolverHighResolutionEuler = TRUE;
  }
  else if( strncmp( &argv[argc][2], STR_ODE_SOLVER_METHOD_ADAPTIVE_EULER, sizeof( STR_ODE_SOLVER_METHOD_ADAPTIVE_EULER )) == 0 ) {
    pCntxt->args.fOdeSolverAdaptiveEuler = TRUE;
  }
  else {
    printf( "ERROR: Parameter error, invalid ODE solver method\n" );
    exitProgram( pCntxt );
  }    

  return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
//   Set model parameters depending on the neuron type
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
void processOPT_N( CNTXT* pCntxt, int argc, char** argv )
{
  if( strncmp( &argv[argc][2], STR_NEURON_TYPE_RS, sizeof( STR_NEURON_TYPE_RS )) == 0 ) {
    pCntxt->args.neuronType = NEURON_TYPE_RS;
    strncpy( pCntxt->neuronTypeName, RS_TYPENAME, sizeof( RS_TYPENAME ));
    pCntxt->neuronParam_A = RS_A;
    pCntxt->neuronParam_B = RS_B;
    pCntxt->neuronParam_C = RS_C;
    pCntxt->neuronParam_D = RS_D;
    pCntxt->neuronParam_THR = RS_THR;
    pCntxt->neuronParam_I_EXT = RS_I_EXT;
    pCntxt->neuronParam_V_INIT = RS_V_INIT;
    pCntxt->neuronParam_U_INIT = RS_U_INIT;
  }
  else if( strncmp( &argv[argc][2], STR_NEURON_TYPE_IB, sizeof( STR_NEURON_TYPE_IB )) == 0 ) {
    pCntxt->args.neuronType = NEURON_TYPE_IB;
    strncpy( pCntxt->neuronTypeName, IB_TYPENAME, sizeof( IB_TYPENAME ));
    pCntxt->neuronParam_A = IB_A;
    pCntxt->neuronParam_B = IB_B;
    pCntxt->neuronParam_C = IB_C;
    pCntxt->neuronParam_D = IB_D;
    pCntxt->neuronParam_THR = IB_THR;
    pCntxt->neuronParam_I_EXT = IB_I_EXT;
    pCntxt->neuronParam_V_INIT = IB_V_INIT;
    pCntxt->neuronParam_U_INIT = IB_U_INIT;
  }
  else if( strncmp( &argv[argc][2], STR_NEURON_TYPE_CH, sizeof( STR_NEURON_TYPE_CH )) == 0 ) {
    pCntxt->args.neuronType = NEURON_TYPE_CH;
    strncpy( pCntxt->neuronTypeName, CH_TYPENAME, sizeof( CH_TYPENAME ));
    pCntxt->neuronParam_A = CH_A;
    pCntxt->neuronParam_B = CH_B;
    pCntxt->neuronParam_C = CH_C;
    pCntxt->neuronParam_D = CH_D;
    pCntxt->neuronParam_THR = CH_THR;
    pCntxt->neuronParam_I_EXT = CH_I_EXT;
    pCntxt->neuronParam_V_INIT = CH_V_INIT;
    pCntxt->neuronParam_U_INIT = CH_U_INIT;
  }
  else if( strncmp( &argv[argc][2], STR_NEURON_TYPE_FS, sizeof( STR_NEURON_TYPE_FS )) == 0 ) {
    pCntxt->args.neuronType = NEURON_TYPE_FS;
    strncpy( pCntxt->neuronTypeName, FS_TYPENAME, sizeof( FS_TYPENAME ));
    pCntxt->neuronParam_A = FS_A;
    pCntxt->neuronParam_B = FS_B;
    pCntxt->neuronParam_C = FS_C;
    pCntxt->neuronParam_D = FS_D;
    pCntxt->neuronParam_THR = FS_THR;
    pCntxt->neuronParam_I_EXT = FS_I_EXT;
    pCntxt->neuronParam_V_INIT = FS_V_INIT;
    pCntxt->neuronParam_U_INIT = FS_U_INIT;
  }
  else if( strncmp( &argv[argc][2], STR_NEURON_TYPE_TC, sizeof( STR_NEURON_TYPE_TC )) == 0 ) {
    pCntxt->args.neuronType = NEURON_TYPE_TC;
    strncpy( pCntxt->neuronTypeName, TC_TYPENAME, sizeof( TC_TYPENAME ));
    pCntxt->neuronParam_A = TC_A;
    pCntxt->neuronParam_B = TC_B;
    pCntxt->neuronParam_C = TC_C;
    pCntxt->neuronParam_D = TC_D;
    pCntxt->neuronParam_THR = TC_THR;
    pCntxt->neuronParam_I_EXT = TC_I_EXT;
    pCntxt->neuronParam_V_INIT = TC_V_INIT;
    pCntxt->neuronParam_U_INIT = TC_U_INIT;
  }
  else if( strncmp( &argv[argc][2], STR_NEURON_TYPE_RZ, sizeof( STR_NEURON_TYPE_RZ )) == 0 ) {
    pCntxt->args.neuronType = NEURON_TYPE_RZ;
    strncpy( pCntxt->neuronTypeName, RZ_TYPENAME, sizeof( RZ_TYPENAME ));
    pCntxt->neuronParam_A = RZ_A;
    pCntxt->neuronParam_B = RZ_B;
    pCntxt->neuronParam_C = RZ_C;
    pCntxt->neuronParam_D = RZ_D;
    pCntxt->neuronParam_THR = RZ_THR;
    pCntxt->neuronParam_I_EXT = RZ_I_EXT;
    pCntxt->neuronParam_V_INIT = RZ_V_INIT;
    pCntxt->neuronParam_U_INIT = RZ_U_INIT;
  }
  else if( strncmp( &argv[argc][2], STR_NEURON_TYPE_LTS, sizeof( STR_NEURON_TYPE_LTS )) == 0 ) {
    pCntxt->args.neuronType = NEURON_TYPE_LTS;
    strncpy( pCntxt->neuronTypeName, LTS_TYPENAME, sizeof( LTS_TYPENAME ));
    pCntxt->neuronParam_A = LTS_A;
    pCntxt->neuronParam_B = LTS_B;
    pCntxt->neuronParam_C = LTS_C;
    pCntxt->neuronParam_D = LTS_D;
    pCntxt->neuronParam_THR = LTS_THR;
    pCntxt->neuronParam_I_EXT = LTS_I_EXT;
    pCntxt->neuronParam_V_INIT = LTS_V_INIT;
    pCntxt->neuronParam_U_INIT = LTS_U_INIT;
  }  
  else {
    printf( "ERROR: Failed to set neuron type and neuron type properties\n" );
    exitProgram( pCntxt );
  }

  return;
}
