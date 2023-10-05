/*
################################################################################
# File:        dump.cc
# Description: LAGO dump of pulses, mainly intended to debug and develop classes
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
	Input.SetForce(1);
  
  //
  // File reading and processing
  //
  double rt=0., ft=0., tt=0.;
  cerr << "Reading file" << endl;
  while(NbPulses!=-1) {
    NbPulses=Input.ReadOneSecond(&Data,Pulse,MAXPULSEPERSEC);
	for (int i=0; i<NbPulses; i++) {
		rt = Pulse[i].GetRiseTime(2);
		ft = Pulse[i].GetFallTime(2);
		tt = Pulse[i].GetFullTime(2);
		if (rt>0 && ft>0)
			cout << rt << " " << ft << " " << tt << endl;
	}
  }
  return 0;
}