/*
################################################################################
# File:        lago_file.h
# Description: handles different lago data structure files
# Author:      Hernán Asorey
# Email:       asoreyh@gmail.com
# Date:        2012
# 
# Copyright:   [2012] [The LAGO Collaboration]
# License:     BSD-3-Clause
# See the LICENSE file in the project root for full license information.
################################################################################
*/
#ifndef LAGO_FILE_H
#define LAGO_FILE_H

#ifndef LAGO_DEFS_H
#include "lago_defs.h"
#endif

#ifndef LAGO_DATA_H
#include "lago_data.h"
#endif

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
using namespace std;

// FIXME Since it is raw data, printing should be done at lago_file::ReadOneSecond file
int iraw = 0, raw_seconds = 0, raw_limit = 10;
ofstream raw;

class LagoFile {
  public:
    LagoFile() {
      Init();
    }
    FILE * fi;
    int verbose;
    int force;
    int NbPulses;
    int _nextsecond;

    void Init() {
      fi=NULL;
      verbose=1;
      force=0;
      NbPulses=0;
      _nextsecond=0;
    }

    void SetVerbose(int f) {
      verbose=f;
    }

    void SetForce(int f) {
      force=f;
    }

    int Open(char *nfi) {
      char tmpc[256];
      if (strncmp(nfi+strlen(nfi)-4,".bz2",4)==0) {
        if (verbose) 
          fprintf(stderr,"File '%s' ending by bz2, reading it via bzip2\n",nfi);
        snprintf(tmpc,256,"bzip2 -d -c %s",nfi);
        fi=popen(tmpc,"r");
      }
      else if((fi=fopen(nfi,"r")) == NULL) {
        if (verbose)
          fprintf(stderr,"Can't open file '%s' - Switching to SDTIN\n",nfi);
        fi=stdin;
      }
      return 1;
    }

    int ReadOneSecond(LagoGeneric *Data,LagoEvent *Pulse,int MAXPULSEPERSEC) {
      char line[256];
      char tmpbuf[256];
      for (int i=0;i<NbPulses;i++) 
        Pulse[i].Init();
      NbPulses=0;
      Data->second=_nextsecond;
      while(fgets(line,250,fi)) {
        if (iraw && raw_seconds <= raw_limit)
          raw << line;
        if (line[0] == '#') {  // generic data
          switch(line[2]) {  // global switch
            case '#':
              break;  // comment
            case 'v':  // version check
              double ver;
              sscanf(line, "# v %lf\n", &ver);
              if (ver > 1 && ver != DATAVERSION) {
                std::cerr << "This code is prepared to work with LAGO raw data version " << DATAVERSION << std::endl;
                std::cerr << "Raw data was obtained using LAGO version " << ver << std::endl;
                if (!force) {
                  std::cerr << "See you later." << std::endl;
                  exit(1);
                }
              }
              break;
            case 'x':         // extra data
              switch (line[4]) {
                case 's':     // sensors
                  sscanf(line, "# x s %lf C %lf hPa", &(Data->temperature), &(Data->pressure)); // not reading altitude
                  break;
                case 'f':     // clock frequency 
                  sscanf(line, "# x f %d\n", &(Data->clockfrequency));
                  break;
                case 'h':     // time
                  sscanf(line, "# x h %s %s %d\n", tmpbuf, tmpbuf, &_nextsecond);
                  raw_seconds++;
                  if (NbPulses)  
                    return NbPulses;
                  Data->second=_nextsecond;
                  break;
              }
              break;
            case 't':         // time and trigger line
              sscanf(line, "# t %d %d\n", &(Pulse[NbPulses].trigger), &(Pulse[NbPulses].clockcount));
              break;
            case 'c':         // trigger counter, pulse is over
              sscanf(line, "# c %d\n", &(Pulse[NbPulses].counter));
              NbPulses++;
              if (NbPulses>MAXPULSEPERSEC-1) {
                std::cerr << "Error, too many pulses in one second (or missing second flag?)" << std::endl;
                std::cerr << "Overwritting pulses..." << std::endl;
                NbPulses=MAXPULSEPERSEC-1;
              }
              break;
          }
        }
        else { // reading pulse
          int c0=0, c1=0, c2=0;
          sscanf(line, "%d %d %d\n", &c0, &c1, &c2);
          Pulse[NbPulses].fill(c0,c1,c2);
        }
      }
      if (NbPulses) 
        return NbPulses;
      return -1; // end of file
    }
};
#endif
