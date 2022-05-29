/*
 *  izhikevich_model.h
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
// =   IZHIKEVICH MODEL DEFINITION
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// IZHIKEVICH MODEL CONSTANTS SET UP AS GLOBAL VARIABLES
float IZHIKEVICH_A   = 0;
float IZHIKEVICH_B   = 0;
float IZHIKEVICH_C   = 0;
float IZHIKEVICH_D   = 0;
float IZHIKEVICH_THR = 0;

#define TRUE  1
#define FALSE 0

// ODE 1
float dvdt( float  v         // IN
          , float  u         // IN
          , void*  pI )      // IN  external current (optional for the ODE solver)
{
  if( pI == NULL ) {
    printf( "ERROR: external current expected\n" );
    exit( -1 );
  } 
  // SpiNNaker FP error simulation 
  // return( 0.03997802 * v * v + 5.0 * v + 140.0 - u + *(float*)pI );
  return( 0.04 * v * v + 5.0 * v + 140.0 - u + *(float*)pI );
}

// ODE 2
float dudt( float  v         // IN
          , float  u )       // IN
{
  return( IZHIKEVICH_A * ( IZHIKEVICH_B * v - u ));
}

// Threshold
int threshold( float*  v     // IN / OUT
              , float*  u )   // IN / OUT             
{
  if( *v >= IZHIKEVICH_THR ) {
    *v = IZHIKEVICH_C;
    *u += IZHIKEVICH_D;
    return( TRUE );
  }
  else {
    return( FALSE );    
  }
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   IZHIKEVICH MODEL PARAMETERS
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// parameter set for regular spiking
#define RS_A        (float)(0.02)   
#define RS_B        (float)(0.2)
#define RS_C        (float)(-65.0)
#define RS_D        (float)(8.0)
#define RS_THR      (float)(30.0)
#define RS_I_EXT    (float)(5.0)
#define RS_V_INIT   (float)(-75.0)   
#define RS_U_INIT   (float)(0.0)
#define RS_TYPENAME "IZHIKEVICH RS NEURON TYPE"

// parameter set for intrinsically bursting
#define IB_A        (float)(0.02)   
#define IB_B        (float)(0.2)
#define IB_C        (float)(-55.0)
#define IB_D        (float)(4.0)
#define IB_THR      (float)(30.0)
#define IB_I_EXT    (float)(5.0)
#define IB_V_INIT   (float)(-75.0)   
#define IB_U_INIT   (float)(0.0)   
#define IB_TYPENAME "IZHIKEVICH IB NEURON TYPE"

// parameter set for chattering
#define CH_A        (float)(0.02)   
#define CH_B        (float)(0.2)
#define CH_C        (float)(-50.0)
#define CH_D        (float)(2.0)
#define CH_THR      (float)(30.0)
#define CH_I_EXT    (float)(5.0)
#define CH_V_INIT   (float)(-75.0)   
#define CH_U_INIT   (float)(0.0)
#define CH_TYPENAME "IZHIKEVICH CH NEURON TYPE"

// parameter set for fast spiking
#define FS_A        (float)(0.1)   
#define FS_B        (float)(0.2)
#define FS_C        (float)(-65.0)
#define FS_D        (float)(2.0)
#define FS_THR      (float)(30.0)
#define FS_I_EXT    (float)(5.0)
#define FS_V_INIT   (float)(-75.0)   
#define FS_U_INIT   (float)(0.0)
#define FS_TYPENAME "IZHIKEVICH FS NEURON TYPE"

// parameter set for thalamo cortical
#define TC_A        (float)(0.02)   
#define TC_B        (float)(0.25)
#define TC_C        (float)(-65.0)
#define TC_D        (float)(0.05)
#define TC_THR      (float)(30.0)
#define TC_I_EXT    (float)(5.0)
#define TC_V_INIT   (float)(-75.0)   
#define TC_U_INIT   (float)(0.0)
#define TC_TYPENAME "IZHIKEVICH TC NEURON TYPE"

// parameter set for resonator
// requires extra current step at e.g. 150 ms
#define RZ_A        (float)(0.1)   
#define RZ_B        (float)(0.25)
#define RZ_C        (float)(-65.0)
#define RZ_D        (float)(2.0)
#define RZ_THR      (float)(30.0)
#define RZ_I_EXT    (float)(0.7)
#define RZ_V_INIT   (float)(-75.0)   
#define RZ_U_INIT   (float)(0.0)
#define RZ_TYPENAME "IZHIKEVICH RZ NEURON TYPE"

// parameter set for low-threshold spiking
#define LTS_A        (float)(0.02)   
#define LTS_B        (float)(0.25)
#define LTS_C        (float)(-65.0)
#define LTS_D        (float)(2.0)
#define LTS_THR      (float)(30.0)
#define LTS_I_EXT    (float)(5.0)
#define LTS_V_INIT   (float)(-75.0)   
#define LTS_U_INIT   (float)(0.0)
#define LTS_TYPENAME "IZHIKEVICH LTS NEURON TYPE"