/*
################################################################################
# File:        lago_data.h
# Description: Classes, methods and attributes for handling LAGO data defs.
# Author:      Hernán Asorey
# Email:       asoreyh@gmail.com
# Date:        2012
# 
# Copyright:   [2012] [The LAGO Collaboration]
# License:     BSD-3-Clause
# See the LICENSE file in the project root for full license information.
################################################################################
*/
#ifndef LAGO_DATA_H
#define LAGO_DATA_H

#ifndef LAGO_DEFS_H
#include "lago_defs.h"
#endif

#include <iostream>
#include <math.h>

class LagoGeneric {
  public:
    int second;
    double temperature;
    double pressure;
    int clockfrequency;
    int current;

    void Init() {
      second = 0;
      temperature = 0.;
      pressure = 0.;
      clockfrequency = 0;
      current = 0;
    }
};

class LagoEvent {
  public:
    LagoEvent() {
      Init();
    };
    int trigger;
    int counter;
    int clockcount;
    int trace[CHANNELS][TRACELEN];
    int currentbinfilled;

    int fill(int v1, int v2, int v3) {
      if (currentbinfilled<TRACELEN) {
        trace[0][currentbinfilled]=v1;
        trace[1][currentbinfilled]=v2;
        trace[2][currentbinfilled]=v3;
        currentbinfilled++;
        return 1;
      }
      return 0;
    }

    int nanotime(int clockfrequency) {
      return (int)(floor(((double)(1.e9*clockcount))/clockfrequency));
    }

    void Init() {
      currentbinfilled=0;
    }

    int IsTriggered(int channel, int level=-1) { // 0, 1, or 2
      if (level==-1) // use electronic trigger
        return ((trigger & (1<<channel))==(1<<channel));
      else { // use specific level
        return 0; // FIXME implement code
      }
    }

    int GetPeak(int channel) {
      int peak=0, aux=0;
      for (int i=0;i<currentbinfilled;i++) {
        aux = trace[channel][i]-BASELINE;
        if (aux > peak)
          peak=aux;
      }
      return peak;
    }

	double GetPulseBase(int channel) {
		/* get first and the last two bins and return the average as 
		 * the bl for this pulse
		 * */
		//double bl=trace[channel][0]+trace[channel][TRACELEN-2]+trace[channel][TRACELEN-1];
		//return (bl/3.);
		return trace[channel][0];
	}

    int GetBase(int channel) {
      return GetPulseBase(channel);
    }
	
	int GetValAtPos(int channel, int pos) {
      return (trace[channel][pos - 1] - BASELINE); 
    }

    int GetValAtTrigger(int channel) {
      return (GetValAtPos(channel, TRIGGERBIN)); // Pulse triggered at 3rd bin
    }
    
    int GetCharge(int channel, int negativepulse=0, int max=4095) {
      int charge=0;
      if (negativepulse) {
        // negative pulse for Boyita and other detectors where pulse
        // can become negative due to the shaper response
        for (int i=0;i<currentbinfilled;i++) 
          if (trace[channel][i]>BASELINE) 
			  charge+=trace[channel][i]-BASELINE;
          else 
			  charge+=BASELINE-trace[channel][i];
      } else {
        for (int i=0;i<currentbinfilled;i++) 
          charge+=trace[channel][i];
        charge -= currentbinfilled * BASELINE;
      }
      if (charge>max) 
        charge=max;
      if (charge<0) 
        charge=0;
      return charge;
    }

	double GetTfraction(int channel, int fraction) {
		if (fraction < 10 || fraction > 100)
			return -1;
		double charge = (double) GetCharge(channel);
		if (!charge)
			return -2;
		int i=0;
		double s = 0.;
		double dy = 0.;
		double x = 0.;
		double limit = fraction / 100.;
		for (i=0; i<currentbinfilled; i++) {
			dy = 1.0 * (trace[channel][i] - BASELINE) / charge;
			s += dy;
			if (s >= limit)
				break;
		}
		if (!i)
			return -3;
		if (!dy)// no changes from previous bin? should not happen thanks to the equal sign in the comparisson, but just in case
			return ((i-1.)*BIN); // return previous bin
		x = (1.0 * i - ((s - limit) / dy)) * BIN;
		return x;
	}

	double GetCharTime(int channel, int up, int low) {
		if (up < low)
			return -1; // stupid user
		double Tu=GetTfraction(channel,up);
		double Tl=GetTfraction(channel,low);
		if (Tu >= 0 && Tl >= 0) 
			return (Tu-Tl);
		else
			return -1;
	}
	
	double GetRiseTime (int channel) {
		return GetCharTime(channel,50,10);
	}
	
	double GetFallTime (int channel) {
		return GetCharTime(channel,90,50);
	}
	
	double GetFullTime (int channel) {
		return GetCharTime(channel,90,10);
	}
    
	void dump() {
      std::cout 
		  << "# trg: " << trigger << " cnt: " << counter << " clk: " << clockcount << " " 
		  << std::endl;
	  for (int i=0; i<currentbinfilled; i++) {
		  std::cout << "# " << i+1 << ": ";
		  for (int j=0; j<CHANNELS; j++) {
			  std::cout << trace[j][i] << " ";
		  }
		  std::cout << std::endl;
	  }
	  std::cout << "# bl: ";
	  for (int j=0; j<CHANNELS; j++)
			  std::cout << GetPulseBase(j) << " "; 
	  std::cout << std::endl;
	  std::cout << "# pk: ";
	  for (int j=0; j<CHANNELS; j++)
			  std::cout << GetPeak(j) << " "; 
	  std::cout << std::endl;
	  std::cout << "# cg: ";
	  for (int j=0; j<CHANNELS; j++)
			  std::cout << GetCharge(j) << " "; 
	  std::cout << std::endl;
	  std::cout << std::endl;
    }
};

#endif
