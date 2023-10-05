/*
################################################################################
# File:        example.cc
# Description: LAGO example on how to include and use LAGO ANNA main classes
# Author:      Hernán Asorey
# Email:       asoreyh@gmail.com
# Date:        2012
# 
# Copyright:   [2012] [The LAGO Collaboration]
# License:     BSD-3-Clause
# See the LICENSE file in the project root for full license information.
################################################################################
*/
#define _FILE_OFFSET_BITS 64

#include "lago_defs.h"
#include "lago_data.h"
#include "lago_file.h"

using namespace std;

#define MAXPULSEPERSEC 1000000

// This is just an example to use LAGO libraries

int main (int argc, char **argv) {
	if (argc < 2) {
		cerr << "Error: data file is needed. Please check." << endl;
		cerr << endl;
		cerr << "Usage: " << argv[0] << " <file>" << endl;
		return 1;
	}
	LagoFile Input;
	int NbPulses=0;
	LagoGeneric Data;
	LagoEvent *Pulse;
	Pulse=(LagoEvent*)malloc(MAXPULSEPERSEC*sizeof(LagoEvent));
	for (int i=0;i<MAXPULSEPERSEC;i++) 
		Pulse[i].Init();
	Input.Open(argv[1]);
  
  //
  // File reading and processing
  //
  cerr << "Reading file" << endl;
  while(NbPulses!=-1) {
    NbPulses=Input.ReadOneSecond(&Data,Pulse,MAXPULSEPERSEC);
    cout << Data.second << " " << Data.pressure << " " << Data.temperature << " " << NbPulses << endl;
	/* 
	 * So, Data have one second of pulses, here we should our analysis
	 * say, looking for excesses for grb analysis, something like this:
	 * if (NbPulses>0)
	 *		TreatSecond(&Data,Pulse,NbPulses);
	 */
  }
  return 0;
}