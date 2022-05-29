/*
 *  ode_solver.h
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
//  ODE solver implementations
//
//  References:
//
//  [1] Izhikevich, E. M. (2003). Simple model of spiking neurons. Trans. Neur. Netw. 14, 1569–1572.
//      doi:10.1109/TNN.2003.820440
//
//  [2] Hopkins, M., and Furber, S. (2015). Accuracy and efficiency in fixed-point neural ode solvers.
//      Neural Comput. 27, 2148–2182. doi: 10.1162/NECO_a_00772
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +



#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   Explicit Forward Euler
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void ode2dSolver_StandardEuler( PFUNCD  pODE1dxdt          // IN        function pointers to the 2d coupled ODE system
                              , PFUNCD  pODE2dydt          // IN 
                              , float   h                  // IN        step width
                              , float   x_t                // IN        x(t)
                              , float   y_t                // IN        y(t)
                              , float*  x_tplus1           // OUT       x(t+1)
                              , float*  y_tplus1           // OUT       y(t+1)
                              , void*   optParm = NULL )   // IN / OUT  Optional parameters
{
  *x_tplus1 = x_t + h * pODE1dxdt( x_t, y_t, optParm );
  *y_tplus1 = y_t + h * pODE2dydt( x_t, y_t, optParm );

  return;
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   Symplectic Euler
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void ode2dSolver_SymplecticEuler( PFUNCD  pODE1dxdt          // IN        function pointers to the 2d coupled ODE system
                                , PFUNCD  pODE2dydt          // IN 
                                , float   h                  // IN        step width
                                , float   x_t                // IN        x(t)
                                , float   y_t                // IN        y(t)
                                , float*  x_tplus1           // OUT       x(t+1)
                                , float*  y_tplus1           // OUT       y(t+1)
                                , void*   optParm = NULL )   // IN / OUT  Optional parameters
{
  *x_tplus1 = x_t + h * pODE1dxdt( x_t, y_t, optParm );
  *y_tplus1 = y_t + h * pODE2dydt( *x_tplus1, *y_tplus1, optParm );    // symplectic

  return;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   Explicite solver reduction, SpiNNaker implementations, see [2]
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void explizitSpiNNaker( float   h
                      , float*  v
                      , float*  u
                      , float   i )
{  
  float lastV1 = *v;
  float lastU1 = *u;
  float a = IZHIKEVICH_A;
  float b = IZHIKEVICH_B;
  float input_this_timestep = i;

  float pre_alph = (140.0) + input_this_timestep - lastU1;
  float alpha    = pre_alph + ( (5.0) + (0.0400) * lastV1) * lastV1;
  float eta = lastV1 + 0.5*(h * alpha);

  float beta = 0.5*(h * (b * lastV1 - lastU1) * a); 

  *v += h * (pre_alph - beta + ( (5.0) + (0.0400) * eta) * eta);

  *u += a * h * (-lastU1 - beta + b * eta);         

  threshold( v, u );  
  return;
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   Forward Euler, NEST implemenattion
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void explizitNEST1( float   h
                  , float*  v
                  , float*  u
                  , float   i )
{  
  float v_old = *v;
  float u_old = *u;
    
  float a = IZHIKEVICH_A;
  float b = IZHIKEVICH_B;

  *v += h * ( 0.04 * v_old * v_old + 5.0 * v_old + 140.0 - u_old + i );  
  *u += h * a * ( b * v_old - u_old );

  threshold( v, u );  
  return;
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   Forward Euler, original Izhikevich implementation, see [1]
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void explizitNEST2( float   h
                  , float*  v
                  , float*  u
                  , float   i )
{
    float v_old = *v;
    float u_old = *u;
    
    float a = IZHIKEVICH_A;
    float b = IZHIKEVICH_B;

    *v += h * 0.5 * ( 0.04 * v_old * v_old + 5.0 * v_old + 140.0 - u_old + i );  // for numetrical stablility
    *v += h * 0.5 * ( 0.04 * v_old * v_old + 5.0 * v_old + 140.0 - u_old + i );  // see [1]
    *u += h * a * ( b * v_old - u_old );

    threshold( v, u );
    return;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   GNU Scientific Library, Runge-Kutta-Fehlberg (rkf45)
// =   (based on the code generated by NESTML)
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
int gslODESolver_OdeSystem( double t, const double *y, double *f, void* params );

const int V_m_INDEX = 0;
const int U_m_INDEX = 1;

float param_I_e = 0;

static int countGSL = 0;

// gsl_odeiv.h
// typedef struct  {
//   int (* function) (float t, const float y[], float dydt[], void * params);
//   int (* jacobian) (float t, const float y[], float * dfdy, float dfdt[], void * params);
//   size_t dimension;
//   void * params;
// }
// gsl_odeiv_system;
gsl_odeiv_system    systemOfODEs = { 0 };
gsl_odeiv_step*     pSteppingFunction  = NULL;
gsl_odeiv_control*  pControlObject     = NULL;
gsl_odeiv_evolve*   pEvolutionFunction = NULL;


void gslODESolver_Integrate( float   h
                           , float*  v
                           , float*  u
                           , float   i )
{
  int    gslRc = GSL_SUCCESS;
  double t = 0;
  double t1 = h;
  double stateVector[2];
  double stepSize = h;          // the step size is adjusted by the evolve-function

  param_I_e = i;

  stateVector[V_m_INDEX] = *v;
  stateVector[U_m_INDEX] = *u;
  


  while( t < t1 ) {                                        // There are 3 loops !!!
                                                           // - The outer loop of the simulation.
                                                           // - This loop: t < stepSize < t1, where the step size is adjusted in each iteration.
                                                           // - The loop within 'evolve-function' calling the ODE's 'systemOfODEs.function = &gslODESolver_OdeSystem'
    gslRc = gsl_odeiv_evolve_apply( pEvolutionFunction
                                  , pControlObject         // -> accuracy
                                  , pSteppingFunction      // -> algorithm
                                  , &systemOfODEs          // -> ODEs
                                  , &t
                                  , t1
                                  , &stepSize
                                  , stateVector );         // the state vectors need to be updated in the equations as they are reset between steps

    if( gslRc != GSL_SUCCESS ) {
      printf( "ERROR: gsl error\n" );
      exit( -1 );
    }

    // threshold detection
    if( stateVector[V_m_INDEX] >= IZHIKEVICH_THR ) {
      stateVector[V_m_INDEX] = IZHIKEVICH_C;
      stateVector[U_m_INDEX] += IZHIKEVICH_D;
    }

  }

  // new state
  *v = stateVector[V_m_INDEX];
  *u = stateVector[U_m_INDEX];

  return;
}

void gslODESolver_Init()
{
  pSteppingFunction = gsl_odeiv_step_alloc( gsl_odeiv_step_rkf45, 2 );      // instance of a stepping function for a 2D system of ODEs
  pControlObject = gsl_odeiv_control_y_new( 1e-6, 0.0 );                    // create a new control object keeping the local error in each step
                                                                            // within an absolute error of 1e-6 and a relative error of 0.0 with respect to
                                                                            // the solution y_i(t)
  pEvolutionFunction = gsl_odeiv_evolve_alloc ( 2 );                        // instance of an evolution function for a 2D system of ODEs
  
  systemOfODEs.function = &gslODESolver_OdeSystem;
  systemOfODEs.jacobian = NULL;
  systemOfODEs.dimension = 2;
  systemOfODEs.params = &param_I_e;
 
  return;
}

int gslODESolver_OdeSystem( double t, const double *y, double *f, void* params )
{
  float I_e = *(float*)params;

  f[ V_m_INDEX ] = 0.04 * y[V_m_INDEX] * y[V_m_INDEX] + 5.0 * y[V_m_INDEX] + 140 - y[U_m_INDEX] + I_e;
  f[ U_m_INDEX ] = IZHIKEVICH_A * ( IZHIKEVICH_B * y[V_m_INDEX] - y[U_m_INDEX] );

  countGSL++;
  printf( "countGSL: %d\n", countGSL );

  return GSL_SUCCESS;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   Heun
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void ode2dSolver_Heun( PFUNCD  pODE1dxdt          // IN        function pointers to the 2d coupled ODE system
                     , PFUNCD  pODE2dydt          // IN 
                     , float   h                  // IN        step width
                     , float   x_t                // IN        x(t)
                     , float   y_t                // IN        y(t)
                     , float*  x_tplus1           // OUT       x(t+1)
                     , float*  y_tplus1           // OUT       y(t+1)
                     , void*   optParm = NULL )   // IN / OUT  Optional parameters
{
  float x_Euler = x_t + h * pODE1dxdt( x_t, y_t, optParm );
  float y_Euler = y_t + h * pODE2dydt( x_t, y_t, optParm );
  *x_tplus1 = x_t + h * pODE1dxdt( 0.5 * (x_t + x_Euler), 0.5 * (y_t + y_Euler), optParm );
  *y_tplus1 = y_t + h * pODE2dydt( 0.5 * (x_t + x_Euler), 0.5 * (y_t + y_Euler), optParm );

  return;
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   Standard Euler with an adaptive step size
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void ode2dSolver_AdaptiveEuler( PFUNCD  pODE1dxdt          // IN        function pointers to the 2d coupled ODE system
                              , PFUNCD  pODE2dydt          // IN 
                              , float   h                  // IN        step width
                              , float   x_t                // IN        x(t)
                              , float   y_t                // IN        y(t)
                              , float*  x_tplus1           // OUT       x(t+1)
                              , float*  y_tplus1           // OUT       y(t+1)
                              , void*   optParm = NULL )   // IN / OUT  Optional parameters
{
  float t        = 0;
  float stepSize = h;   // initial step size, will be adjusted in every step
  float x_interm = x_t;
  float y_interm = y_t;

  // determine integration step size
  float error_accepted   = 0.2;
  float error_calculated = 0;
  int   k = 1;
  float x_hi;

  static float simulationTime_ms = 0;
  int count1 = 0;
  int count2 = 0;
  static int count3 = 0;
  static int count4 = 0;

  float tempError[32] = { 0.0 };
  int index = 0;
  int i = 0;
  int fSpike = FALSE;
  
  if( count4 == 0 ) {
    printf( "[INITIAL VALUES]      h: %11.6f    x_t: %11.6f    y_t: %11.6f   \n", h, x_t, y_t );
  }
  
  do {
    stepSize = h/k;
   
    x_hi = x_t + stepSize * dvdt( x_t, y_t, optParm ); 

    error_calculated =  fabs(fabs(x_hi) - fabs(x_t));

    if( index == 31 ) {
      tempError[ index ] = 99.99;
    }
    else {
      tempError[ index ] = error_calculated;
      index++;
    }


    k *= 2;
    count1++;
  }  while(!(( error_calculated < error_accepted ) || ( stepSize <= 0.001 )));

  while( t < h ) {

    x_interm = x_interm + stepSize * dvdt( x_interm, y_interm, optParm );
    y_interm = y_interm + stepSize * dudt( x_interm, y_interm );
    fSpike = threshold( &x_interm, &y_interm );

    t += stepSize;
    count2++;
    count3++;

    if( fSpike ) {
      printf( "[%8.3f]   > > > S P I K E < < <   +  %11.6f ms \n", simulationTime_ms, t );
    }

  };

  #if( 0 )
    printf( "[%8.3f]   V: %11.6f   U: %11.6f   stepSize: %11.6f   CNT(find stepSize): %6d   CNT(calc): %6d   CNT(total): %6d "
          , simulationTime_ms, x_interm, y_interm, stepSize, count1, count2, count3 );
    printf( "      ERRORs: " );
    for( i = 0; i < index; ++i ) {
      printf( "%11.6f, ", tempError[ i ] );
    }
    printf( "\n" );
  #endif

  *x_tplus1 = x_interm;
  *y_tplus1 = y_interm;

  count4++;
  simulationTime_ms += h;
  return;
}
