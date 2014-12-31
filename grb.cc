/*  grb.cc  --  LAGO grb analysis program
 *  
 *  Copyright (C) 2012-TODAY The LAGO Project, http://lagoproject.org, lago-pi@lagoproject.org
 *
 *  Original authors: Hernán Asorey, Xavier Bertou
 *  e-mail: asoreyh@cab.cnea.gov.ar  (asoreyh@gmail.com)
 *  Laboratorio de Detección de Partículas y Radiación (LabDPR)
 *  Centro Atómico Bariloche - San Carlos de Bariloche, Argentina */

/*  LICENSE GPLv3
 *  This program is free software: you can redistribute it and/or modify 
 *  it under the terms of the GNU General Public License as published by 
 *  the Free Software Foundation, either version 3 of the License, or 
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>. */

/*  -*- coding: utf8 -*-
 *  try to preserve encoding  */
/************************************************************************/
#define _FILE_OFFSET_BITS 64

#include "lago_data.h"
#include "lago_file.h"

using namespace std;

#define MAXPULSEPERSEC 1000000

// This is just an example to use LAGO libraries

int main (int argc, char **argv) {
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
	/* TODO
    if (NbPulses>0) 
		TreatSecond(&Data,Pulse,NbPulses);
	*/
  }
  return 0;
}
