#ifndef __PARAMS_SIMULATE_SELECTED_STATE_1_H__
#define __PARAMS_SIMULATE_SELECTED_STATE_1_H__
/*
 *  params_simulate_selected_state_1.h
 *
 *  Copyright (C) 2018*  This file is part of the refactored Izhikevich polychronization model application.
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
// = Reload the 1st selected network state, run the simulation for 60 seconds, and record the network activity data.
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#define PARAM_INFO_STRING                     "RELOAD 1st NETWORK STATE AND SIMULATE W/O STDP"

#define SIM_TIME                              (int)(60)      // in seconds

#define GENERATE_NETWORK_FROM_EXTERNAL_DATA   true
#define USE_EXTERNAL_STIMULUS                 true
#define USE_STDP                              false
#define INIT_WITH_RANDOM_WEIGHTS              false

#define INFILE_CONNECTION_MATRIX              "../../data/conMatrix.dat"
#define INFILE_DELAY_MATRIX                   "../../data/delayMatrix.dat"
#define INFILE_STIMULUS                       "../../data/randomNetworkInput.dat"
#define INFILE_WEIGHT_MATRIX                  "../../data/weightMatrix_after1h.dat"

#define OUTFILE_FIRINGS                       "../../data/firings_after1h.dat"

// can be specified for verification purposes
// #define OUTFILE_CONNECTIONS                "../../data/conMatrix_initial.dat"
// #define OUTFILE_WEIGHTS_INITIAL            "../../data/weightMatrix_initial.dat"
// #define OUTFILE_DELAYS                     "../../data/delayMatrix_initial.dat"

// generate connection lists for import into PyNN for SpiNNaker
#define OUTFILE_CON_WEIGHT_DELAY_INITIAL      "../../data/pythonConWeightDelay_after1h.py"

// generate connect function calls for import into FPGA NC node APC program
#define OUTFILE_HNC_NODE_CONNECT_CALLS        "../../data/hnc_node_gen_network_after1h.h"

// generate dictionary with connection information for import into NEST
#define OUTFILE_NEST_CONNECTIVITY_DICTS       "../../data/dictsConWeightDelayNEST_after1h.py"

#endif   // __PARAMS_SIMULATE_SELECTED_STATE_1_H__
