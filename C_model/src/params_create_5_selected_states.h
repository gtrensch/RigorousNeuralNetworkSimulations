#ifndef __PARAMS_CREATE_5_SELECTED_STATES_H__
#define __PARAMS_CREATE_5_SELECTED_STATES_H__
/*
 *  params_create_5_selected_states.h
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

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = P R O G R A M   P A R A M E T E R S
// =
// = Select five network states, after 1, 2, 3, 4, and 5 hours.
// = Run the simulation for 5+ hours with STDP on, and save the five network states, that is, the connection and
// = delay matrix (same for all states), and the weight matrices: w(t1), ... w(t5).
// = Additionally, the random input data is recorded as well as the network activity data.
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#define PARAM_INFO_STRING                     "CREATE FIVE NETWORK STATES WITH STDP ON"

#define SIM_TIME                              (int)(60*60*5+120)      // seconds

#define GENERATE_NETWORK_FROM_EXTERNAL_DATA   false
#define USE_EXTERNAL_STIMULUS                 false
#define USE_STDP                              true
#define INIT_WITH_RANDOM_WEIGHTS              false

#define SELECTED_STATE_1_AFTER_N_SECONDS      (int)(60*60*1)          // seconds
#define SELECTED_STATE_2_AFTER_N_SECONDS      (int)(60*60*2)
#define SELECTED_STATE_3_AFTER_N_SECONDS      (int)(60*60*3)
#define SELECTED_STATE_4_AFTER_N_SECONDS      (int)(60*60*4)
#define SELECTED_STATE_5_AFTER_N_SECONDS      (int)(60*60*5)

#define OUTFILE_CONNECTIONS                   "../../data/conMatrix.dat"
#define OUTFILE_DELAYS                        "../../data/delayMatrix.dat"

#define OUTFILE_STATE_1_WEIGHT_MATRIX         "../../data/weightMatrix_after1h.dat"
#define OUTFILE_STATE_2_WEIGHT_MATRIX         "../../data/weightMatrix_after2h.dat"
#define OUTFILE_STATE_3_WEIGHT_MATRIX         "../../data/weightMatrix_after3h.dat"
#define OUTFILE_STATE_4_WEIGHT_MATRIX         "../../data/weightMatrix_after4h.dat"
#define OUTFILE_STATE_5_WEIGHT_MATRIX         "../../data/weightMatrix_after5h.dat"

#define OUTFILE_STIMULUS                      "../../data/randomNetworkInput.dat"
#define OUTFILE_FIRINGS                       "../../data/firings.dat"

#endif   // __PARAMS_CREATE_5_SELECTED_STATES_H__
