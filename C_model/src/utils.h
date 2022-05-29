#ifndef __UTILS_H__
#define __UTILS_H__
/*
 *  utils.h
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

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   F O R W A R D   D E C L A R A T I O N S
/// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
void ExportConnectionMatrixToFile( const char* pFileName );
void ExportWeightMatrixToFile( const char* pFileName );
void ExportDelayMatrixToFile( const char* pFileName );
void ExportConnectionMatrixWeightAndDelay( const char* pFileName );
void ExportHNCNodeConnectCalls( const char* pFileName );
void ExportNESTConnectionDicts( const char* pFileName );

void ImportConnectionMatrixFromFile( const char* pFileName );
void ImportWeightMatrixFromFile( const char* pFileName );
void ImportDelayMatrixFromFile( const char* pFileName );

void RecordNetworkActivityToFile( const char* pFileName, int simulationSecond, int numFirings );
void RecordRandomStimulusToFile( const char* pFileName, int simTimeSecond, int simTimeMillisecond, int inputNeuron );
int  GetNextExternalStimulusFromFile( const char* pFileName, int simTimeSecond, int t );

int  GetDelayOfConnection( int preSynNeuron, int synapse );
void DeleteFile( const char* pFileName );

void PrintMatrixOfPostSynapticNeurons();
void PrintMatrixOfSynapticWeights();
void PrintMatrixOfSynapticDelays();

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// =   M A C R O   D E F I N I T I O N S
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#define SIZE_OF_ENTRY_NEURON       (int)(5)
#define SIZE_OF_ENTRY_WEIGHT       (int)(11)
#define SIZE_OF_ENTRY_DELAY        (int)(11)

#endif   // __RUN_REFACTORED_VERSION__
#endif   // __UTILS_H__
