#ifndef __PARAMS_H__
#define __PARAMS_H__
/*
 *  params.h
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

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   S I M U L A T I O N   S E T U P S
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#include "params_create_5_selected_states.h"
// #include "params_simulate_selected_state_1.h"
// #include "params_simulate_selected_state_2.h"
// #include "params_simulate_selected_state_3.h"
// #include "params_simulate_selected_state_4.h"
// #include "params_simulate_selected_state_5.h"


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   S O U R C E   C O D E   V E R S I O N   A N D   D E B U G
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#define __DEBUG__                  false        // print connection, weight and delay matrix
#define __RUN_REFACTORED_VERSION__ true         // run the refactored (i.e., this) version of poly_spnet.cpp
                                                // set to false to run the original version from https://www.izhikevich.org/publications/spnet.htm
#define LOG_MIN_MAX_V_U            false
#define REPORT_ON_MULTIPLE_SPIKES  false
#define RECORD_SPIKE_DATA          true

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   O D E   S O L V E R   P A R A M E T E R S
// =
// =   The parameters also apply to the poly_spnet.cpp source.
// =   It has been modified to be able to run with the refined ODE solver and
// =   test for the number of polychronous groups.
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#define ODE_SOLVER_REFINEMENT      true
#define ODE_SOLVER_NOT_SYMPLECTIC  true
#define ODE_SOLVER_STEPS           (int)(10)

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   N E T W O R K   G E N E R A T I O N   P A R A M E T E R S
// =
// =   These parameters can be set to create a network with a better balanced number of
// =   outgoing connections.
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#define BALANCED_OUTGOING_CONNS    false
#define LIMIT_OUTGOING_CONNS       (int)(107)

#endif   // __PARAMS_H__
