// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = Set up desired ODE solver implementation
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#define IMPL01  false    // original SpiNNaker ESR implementation (corresponds to Iteration I)
#define IMPL02  false    // original Izhikevich implementation
#define IMPL03  false    // fixed-step size (h=1/16) symplectic forward Euler, precise threshold detection
#define IMPL04  true     // fixed-step size (h=1/16) symplectic forward Euler, precise threshold detection, FP conversion (corresponds to Iteration III)

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



#include "neuron_model_izh_impl.h"
#include <debug.h>

static global_neuron_params_pointer_t global_params;

/*! \brief For linear membrane voltages, 1.5 is the correct value. However
 * with actual membrane voltage behaviour and tested over an wide range of
 * use cases 1.85 gives slightly better spike timings.
 */
static const REAL SIMPLE_TQ_OFFSET = REAL_CONST(1.85);

#if( IMPL01 )
static inline void _rk2_kernel_midpoint(REAL h, neuron_pointer_t neuron,
                                        REAL input_this_timestep) {

    // to match Mathematica names
    REAL lastV1 = neuron->V;
    REAL lastU1 = neuron->U;
    REAL a = neuron->A;
    REAL b = neuron->B;

    REAL pre_alph = REAL_CONST(140.0) + input_this_timestep - lastU1;
    REAL alpha = pre_alph
                 + ( REAL_CONST(5.0) + REAL_CONST(0.0400) * lastV1) * lastV1;
    REAL eta = lastV1 + REAL_HALF(h * alpha);

    // could be represented as a long fract?
    REAL beta = REAL_HALF(h * (b * lastV1 - lastU1) * a);

    neuron->V += h * (pre_alph - beta
                      + ( REAL_CONST(5.0) + REAL_CONST(0.0400) * eta) * eta);

    neuron->U += a * h * (-lastU1 - beta + b * eta);
}
#endif

#if( IMPL02 )
// h = 1.0, i.e., simulation time-step must be set to 1.0 !
static inline void _originalIzhikevich( REAL h, neuron_pointer_t neuron,
                                        REAL input_this_timestep) {
   REAL v = neuron->V;
   REAL u = neuron->U;
   REAL param_a = neuron->A;
   REAL param_b = neuron->B;

   v += REAL_HALF(( 0.04k * v + 5.0k) * v + 140.0k - u + input_this_timestep );   // for numerical stability
   v += REAL_HALF(( 0.04k * v + 5.0k) * v + 140.0k - u + input_this_timestep );   // time step is 0.5 ms
   u += param_a * ( param_b * v - u );

   neuron->V = v;
   neuron->U = u;
}
#endif

#if( IMPL03 )
// h = 1.0, i.e., simulation time-step must be set to 1.0 !
static inline void _fixedStepSizeEuler( REAL h, neuron_pointer_t neuron,
                                        REAL input_this_timestep) {
   #define THRESHOLD_DETECTED  1024.0k      // some high value above threshold to trigger spike event
                                            // Hiding and encoding information in the data to control
                                            // the program flow is not good practice! 
   REAL v = neuron->V;
   REAL u = neuron->U;
   REAL param_a = neuron->A;
   REAL param_b = neuron->B;
   REAL a;
   REAL b;
   bool spikeEvent = false;
   int steps = 16;
 
   for( int cycles = 0; cycles; ++cycles ) {
     a = 0.04k * v;
     a *= v;

     b = 5.0k * v + 140.0k - u + input_this_timestep;
     v += 0.0625k * (a + b);                           // 0.0625 = 1/16
  
     u += 0.0625k * (param_a * ( param_b * v - u ));   // 0.0625 = 1/16
 
     if( v >= 30.0k ) {
       spikeEvent = true;
       v = neuron->C;
       u += neuron->D;
     }
   }
 
   if( spikeEvent ) {
     neuron->V = v + THRESHOLD_DETECTED;
   }
   else { 
     neuron->V = v;
   }
   neuron->U = u;    
}
#endif

#if( IMPL04 )
// h = 1.0, i.e., simulation time-step must be set to 1.0 !
static inline void _fixedStepSizeEulerWithFPConv( REAL h, neuron_pointer_t neuron,
                                                  REAL input_this_timestep) {
   #define THRESHOLD_DETECTED  1024.0k      // some high value above threshold to trigger spike event

   REAL v = neuron->V;
   REAL u = neuron->U;
   REAL param_a = neuron->A;
   REAL param_b = neuron->B;

   REAL param_a_RS = 5.12k;                // 0.02 RS param_a s8.23
   REAL param_a_FS = 25.6k;                // 0.1  FS param_a s8.23
   REAL param_b_FSRS = 51.2k;              // 0.2  FS param_a s8.23
   REAL param_a_s8_23;
   REAL a;
   REAL b;
   REAL c;
   REAL d;
   bool spikeEvent = false;

   int steps = 16;
 
   for( int cycles = 0; cycles < (steps + 2); ++cycles ) {  // (steps + 2) to supress unrolling the loop ..
     if(cycles < steps ) {                                  // .. to prevent ITCM overflow (2 empty cyles)
                                                            // (Corresponding compiler option did not work.)

       a = 10.24k * v;                                      // byte left-shift
       a *= 0.00390625k;                                    // byte right-shift
       a *= v;

       b = 5.0k * v + 140.0k - u + input_this_timestep;
       v += 0.0625k * (a + b);                              // 0.0625 = 1/16
  
       if( param_a > 0.05k ) {                              // identify neuron type
         param_a_s8_23 = param_a_FS;
       }
       else {
         param_a_s8_23 = param_a_RS;
       }

       d = param_b_FSRS * v;
       d *= 0.00390625k;                                   // byte right-shift 
       d -= u;

       c = 0.0625k * (param_a_s8_23 * ( d ));              // 0.0625 = 1/16
 
       c *= 0.00390625k;                                   // byte right-shift 
       u += c;

       if( v >= 30.0k ) {
         spikeEvent = true;
         v = neuron->C;
         u += neuron->D;
       }
     }
   }
 
   if( spikeEvent ) {
     neuron->V = v + THRESHOLD_DETECTED;
   }
   else { 
     neuron->V = v;
   }
   neuron->U = u;
}
#endif

void neuron_model_set_global_neuron_params(
        global_neuron_params_pointer_t params) {
    global_params = params;
}


state_t neuron_model_state_update( input_t exc_input, input_t inh_input, input_t external_bias,
                                   neuron_pointer_t neuron) {

    input_t input_this_timestep = exc_input - inh_input + external_bias + neuron->I_offset;

#if( IMPL01 )
    // the best AR update so far
    _rk2_kernel_midpoint(neuron->this_h, neuron, input_this_timestep);
#endif

#if( IMPL02 )
    _originalIzhikevich(neuron->this_h, neuron, input_this_timestep);
#endif

#if( IMPL03 )
    _fixedStepSizeEuler(neuron->this_h, neuron, input_this_timestep);
#endif

#if( IMPL04 )
    _fixedStepSizeEulerWithFPConv(neuron->this_h, neuron, input_this_timestep);
#endif

    neuron->this_h = global_params->machine_timestep_ms;

    return neuron->V;
}


void neuron_model_has_spiked(neuron_pointer_t neuron) {
#if( IMPL01 )
    // reset membrane voltage
    neuron->V = neuron->C;

    // offset 2nd state variable
    neuron->U += neuron->D;

    // simple threshold correction - next timestep (only) gets a bump
    neuron->this_h = global_params->machine_timestep_ms * SIMPLE_TQ_OFFSET;
#endif

#if( IMPL02 )
    neuron->V = neuron->C;
    neuron->U += neuron->D;    
#endif

#if( IMPL03 || IMPL04 )
    neuron->V -= THRESHOLD_DETECTED;
#endif
}

state_t neuron_model_get_membrane_voltage(neuron_pointer_t neuron) {
    return neuron->V;
}

void neuron_model_print_state_variables(restrict neuron_pointer_t neuron) {
    log_debug("V = %11.4k ", neuron->V);
    log_debug("U = %11.4k ", neuron->U);
}

void neuron_model_print_parameters(restrict neuron_pointer_t neuron) {
    log_debug("A = %11.4k ", neuron->A);
    log_debug("B = %11.4k ", neuron->B);
    log_debug("C = %11.4k ", neuron->C);
    log_debug("D = %11.4k ", neuron->D);

    log_debug("I = %11.4k \n", neuron->I_offset);
}

