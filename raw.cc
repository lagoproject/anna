/*
 * raw.cc  --  analyze "raw" (L0) data and produce the preprocessed "preliminary" (L1) data set.
 *
 *  Copyright (C) 2012-TODAY The LAGO Project, http://lagoproject.org, lago-pi@lagoproject.org
 *
 *  Original authors: Hernán Asorey, Xavier Bertou
 *  e-mail: asoreyh@cab.cnea.gov.ar  (asoreyh@gmail.com)
 *  Laboratorio de Detección de Partículas y Radiación (LabDPR)
 *  Centro Atómico Bariloche - San Carlos de Bariloche, Argentina
 */

/* LICENSE BSD-3-Clause
Copyright (c) 2012, The LAGO Project
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*  -*- coding: utf8 -*-
 *  try to preserve encoding  */
/************************************************************************/
#define _FILE_OFFSET_BITS 64

#include <string.h>
#include <string.h>
#include <vector>
#include <fstream>
#include "lago_defs.h"
#include "lago_data.h"
#include "lago_file.h"

using namespace std;

// Local defines
#define MAXPULSEPERSEC 1000000
#define MAXTIMEINVECTOR 40000
#define SCL_LEVELS 4
#define SCL_BINS 200
#define SCL_TIME (1./SCL_BINS)
#define SCL_CLOCK (SCL_TIME/25e-9)

// action flags. Should be improved with a class
int ical=0, itim=0, isol=0, iall=0, imon=0, force=0, ivalid=0;
int izip=0, iscl=0, isclg=0, itrg=0, irte=0, iflx=0, ineg[3]= {0,0,0};

// OLD SCALER ANALYSIS
// levels should be given above baseline or threshold (-lt modifier).
// With the new electronics there is no undershoot analysis
int scl_level[CHANNELS][SCL_LEVELS];
int scl_scalers[CHANNELS][SCL_LEVELS][SCL_BINS];
int scl_default[SCL_LEVELS] = { 30, 50, 100, 300};
int scl_prv_trigger = 0;
int scl_flux[CHANNELS];

// FLUX ANALYSIS
int flx_default = 60, flx_time = 0;
int flx_tots_sum = 0, flx_tots_sum_2 = 0;
int flx_sum[CHANNELS], flx_sum_2[CHANNELS];
double flx_aux = 0.;
double flx_temp_sum = 0., flx_temp_sum_2 = 0.;
double flx_pres_sum = 0., flx_pres_sum_2 = 0.;
double flx_temp_avg = 0., flx_temp_dev = 0.;
double flx_pres_avg = 0., flx_pres_dev = 0.;
double flx_tots_avg = 0., flx_tots_dev = 0.;
double flx_avg[CHANNELS], flx_dev[CHANNELS];

// OFFLINE TRIGGER
int trg_level[CHANNELS];
int trg_default = 85;

// TIME DIFFERENCE (AND ALSO ALL PULSE DATA) ANALYSIS
int tim_pc = 0, tim_pt = 0, tim_dc = 0, tim_dt = 0;
double tim_ap = 0., tim_map = 0.;

// ALL PULSE DATA ANALYSIS
int all_trg = 0;

// USE ONLY PULSES TRIGGERED AT CHANNEL FOR CALIBRATION HISTOGRAMS
int icaltrg = 0; // No! by default
int chargetmp = 0, peaktmp = 0;
// DATA STREAMS
ofstream cal, tim;
FILE *scl, *sol, *all, *rte, *flx, *mon;

// BASELINE ANALYSIS
double mon_avg_bl[CHANNELS];
double mon_av2_bl[CHANNELS];
double mon_pulse_bl;
double mon_avg_bl_tmp;
double mon_dev_bl_tmp;
int mon_bl_counts;

// AUXILIARY ARRAYS
int Peak[CHANNELS][1024];
int Base[CHANNELS][1024];
int Charge[CHANNELS][4096];
int ECharge[CHANNELS][4096];
int Time[CHANNELS][MAXTIMEINVECTOR];

int Peak_minute[CHANNELS][1024];
int Charge_minute[CHANNELS][4096];

void TreatSecond(LagoGeneric *Data, LagoEvent*Pulse, int NbPulses) {
	if (((Data->second)%36)==0)
		cerr << "Done " << ((Data->second)%3600)/36 << "%\r";
	if (!NbPulses) {
		if (iscl)
			fprintf(scl, "# s %d %d %.2f %.2f\n", Data->second, Data->clockfrequency, Data->temperature, Data->pressure);
		return;
	}
	if (iall)
		fprintf(all, "# %d %d %.2f %.2f\n", Data->second, Data->clockfrequency, Data->temperature, Data->pressure);
	if (imon) 
		fprintf(mon, "%d %d %.2f %.2f", Data->second, Data->clockfrequency, Data->temperature, Data->pressure);

	for (int j=0; j<CHANNELS; j++)
		mon_avg_bl[j] = mon_av2_bl[j] = 0.;
    mon_pulse_bl=0.;
	mon_avg_bl_tmp = 0.;
	mon_dev_bl_tmp = 0.;
    mon_bl_counts=0;

	for (int j=0; j<CHANNELS; j++)
		scl_flux[j] = 0;
	int scl_idx = 0;
	int scl_peak = 0;

	// processing pulses
	for (int i=0; i<NbPulses; i++) {
		// impossing external trigger
		int trg_drop = 0;
		if (itrg)
			for (int j=0; j<CHANNELS; j++)
				if (Pulse[i].IsTriggered(j))
					if (Pulse[i].GetValAtTrigger(j) < trg_level[j])
						trg_drop++;
		if (itrg && trg_drop)
			continue;
		//we can use this pulse 
		for (int j=0; j<CHANNELS; j++) {
			peaktmp = Pulse[i].GetPeak(j);
			chargetmp = Pulse[i].GetCharge(j,ineg[j]);
			if (icaltrg) {
				if (Pulse[i].IsTriggered(j)) {
					Peak[j][peaktmp]++;
					Charge[j][chargetmp]++;
				}
			} else {
				Peak[j][peaktmp]++;
				Charge[j][chargetmp]++;
			}
			Base[j][Pulse[i].GetBase(j)]++;
			scl_flux[j]++;
			if (iscl) {
				scl_idx = (int)(Pulse[i].clockcount / SCL_CLOCK);
				scl_peak = peaktmp;
				if (isclg)
					scl_peak += BASELINE;
				for (int k=0; k<SCL_LEVELS; k++) {
					if (scl_peak >= scl_level[j][k]) {
						scl_scalers[j][k][scl_idx]++;
					}
				}
			}
		}

		if (iall || itim) {			
			// for both cases we need to calculate time difference between consecutive pulses...
			tim_dc = Pulse[i].counter - tim_pc;
			if (tim_dc < 0)
				tim_dc += 40000000;
			tim_dt = Pulse[i].clockcount - tim_pt;
			if (tim_dt < 0)
				tim_dt += 40000000;

			if (itim) {// Time difference analysis
				if (tim_dc == 1) { // use this pulse for time analysis
					if (tim_dt < MAXTIMEINVECTOR) { // this pulse could fit in the times array
						for (int j=0; j<CHANNELS; j++) {
							if (Pulse[i].IsTriggered(j)) { // only use triggered channel
								if (tim_map) { //aop filter is enabled
									tim_ap = Pulse[i].GetCharge(j,ineg[j]) / Pulse[i].GetPeak(j);
									if (tim_ap > tim_map) // that's ok
										Time[j][tim_dt]++; //Finally
								}
								else {
									Time[j][tim_dt]++; //Finally
								}
							}
						}
					}
				}
			}
			if (iall) {
				fprintf(all, "%d %d %d", tim_dt * 25, tim_dc, Pulse[i].trigger);
				for (int j=0; j<CHANNELS; j++)
					if (Pulse[i].IsTriggered(j)) 
						fprintf(all, " %d %d", Pulse[i].GetCharge(j,ineg[j]), Pulse[i].GetPeak(j));
				for (int j=0; j<CHANNELS; j++)
					if (Pulse[i].IsTriggered(j)) 
						fprintf(all, " %.2f %.2f", Pulse[i].GetRiseTime(j), Pulse[i].GetFallTime(j));
				fprintf(all, "\n");
			}
			tim_pc = Pulse[i].counter;
			tim_pt = Pulse[i].clockcount;
		} // end of tim and all sector

		if (imon) { // baseline analysis for each pulse 
			for (int j=0; j<CHANNELS; j++) {
				mon_pulse_bl = Pulse[i].GetPulseBase(j); 
				mon_avg_bl[j] += mon_pulse_bl;
				mon_av2_bl[j] += mon_pulse_bl * mon_pulse_bl;
			}
			mon_bl_counts++;
		}
	} // close loop for all pulses

	if (imon) { // print baselines
		if (mon_bl_counts) {
			for (int j=0; j<CHANNELS; j++) {
				mon_avg_bl_tmp = mon_avg_bl[j] / mon_bl_counts;
				mon_dev_bl_tmp = sqrt(mon_av2_bl[j] / mon_bl_counts - mon_avg_bl_tmp * mon_avg_bl_tmp);
				fprintf(mon, " %.3f %.3f", mon_avg_bl_tmp, mon_dev_bl_tmp);
			}
			fprintf(mon, "\n");
		}
		else {
			for (int j=0; j<CHANNELS; j++)
				fprintf(mon, " 0.000 0.000");
			fprintf(mon, "\n");
		}
	}
	if (iflx) {
		if (flx_time == flx_default) {
			// averages
			flx_temp_avg = flx_temp_sum / flx_time;
			flx_pres_avg = flx_pres_sum / flx_time;
			flx_tots_avg = flx_tots_sum / flx_time;
			for (int j=0; j<CHANNELS; j++)
				flx_avg[j] = flx_sum[j] / flx_time;

			// devs
			if (flx_time > 1)
				flx_aux = flx_time / (flx_time - 1);
			else
				flx_aux = 1.;
			flx_temp_dev = sqrt(flx_aux * ( flx_temp_sum_2 / flx_time - flx_temp_avg * flx_temp_avg));
			flx_pres_dev = sqrt(flx_aux * ( flx_pres_sum_2 / flx_time - flx_pres_avg * flx_pres_avg));
			flx_tots_dev = sqrt(flx_aux * ( flx_tots_sum_2 / flx_time - flx_tots_avg * flx_tots_avg));
			for (int j=0; j<CHANNELS; j++)
				flx_dev[j] = sqrt(flx_aux * ( flx_sum_2[j] / flx_time -  flx_avg[j] *  flx_avg[j]));

			// print
			fprintf (flx, "%d %.2f %.2f %.2f", Data->second - int(flx_time / 2.), flx_temp_avg, flx_pres_avg, flx_tots_avg);
			for (int j=0; j<CHANNELS; j++)
				fprintf (flx, " %.2f", flx_avg[j]);
			fprintf (flx, " %.2f %.2f %.2f", flx_temp_dev, flx_pres_dev, flx_tots_dev);
			for (int j=0; j<CHANNELS; j++)
				fprintf (flx, " %.2f", flx_dev[j]);
			fprintf (flx, "\n");
			// clean
			flx_time = 0;
			flx_temp_sum = flx_pres_sum = flx_temp_sum_2 = flx_pres_sum_2 = 0.;
			flx_tots_sum = flx_tots_sum_2 = 0;
			for (int j=0; j<CHANNELS; j++)
				flx_sum[j] = flx_sum_2[j] = 0;
		}
		flx_time++;
		flx_temp_sum += Data->temperature;
		flx_pres_sum += Data->pressure;
		flx_tots_sum  += NbPulses;
		flx_temp_sum_2 += Data->temperature * Data->temperature;
		flx_pres_sum_2 += Data->pressure * Data->pressure;
		flx_tots_sum_2  += NbPulses * NbPulses;
		for (int j=0; j<CHANNELS; j++) {
			flx_sum[j] += scl_flux[j];
			flx_sum_2[j] += scl_flux[j] * scl_flux[j];
		}
	}

	if (isol) {
		if (((Data->second)%60)==0) { // every minute, histograms
			fprintf(sol, "# q %d %d %.2f %.2f\n", Data->second, Data->clockfrequency, Data->temperature, Data->pressure);
			for (int i=0; i<CHANNELS; i++) {
				fprintf(sol, "0 %d %d %d %.2f %.2f", i, Data->second, Data->clockfrequency, Data->temperature, Data->pressure);
				for (int j=0; j<4096; j++) {
					fprintf(sol, " %d", Charge[i][j]-Charge_minute[i][j]);
					Charge_minute[i][j]=Charge[i][j];
				}
				fprintf(sol, "\n");
			}
			fprintf(sol, "# p %d %d %.2f %.2f\n", Data->second, Data->clockfrequency, Data->temperature, Data->pressure);
			for (int i=0; i<CHANNELS; i++) {
				fprintf(sol, "1 %d %d %d %.2f %.2f", i, Data->second, Data->clockfrequency, Data->temperature, Data->pressure);
				for (int j=0; j<1024; j++) {
					fprintf(sol, " %d", Peak[i][j]-Peak_minute[i][j]);
					Peak_minute[i][j]=Peak[i][j];
				}
				fprintf(sol, "\n");
			}
		}
	}

	if (irte) {
		fprintf(rte, "%d %.2f %.2f", Data->second, Data->temperature, Data->pressure);
		fprintf(rte, " %d",NbPulses);
		for (int j=0; j<CHANNELS; j++)
			fprintf(rte, " %d", scl_flux[j]);
		fprintf(rte, "\n");
	}

	if (iscl) {
		fprintf(scl, "# s %d %d %.2f %.2f", Data->second, Data->clockfrequency, Data->temperature, Data->pressure);
		fprintf(scl, " %d",NbPulses);
		for (int j=0; j<CHANNELS; j++) {
			fprintf(scl, " %d", scl_flux[j]);
		}
		fprintf(scl, "\n");
		for (int l=0; l<SCL_BINS; l++) {
			for (int j=0; j<CHANNELS; j++) {
				for (int k=0; k<SCL_LEVELS; k++) {
					fprintf(scl,"%d ", scl_scalers[j][k][l]);
					scl_scalers[j][k][l] = 0;
				}
			}
			fprintf(scl,"\n");
		}
	}
}

void Usage(char *prog, int verbose=0)
{
	cout << "\t" << PROJECT << " " << CODEVERSION << endl;
	cout << endl;
	cout << "\t(c) 2012-Today, The LAGO Project, http://lagoproject.org" << endl;
	cout << "\t(c) 2012, LabDPR, http://labdpr.cab.cnea.gov.ar" << endl;
	cout << endl;
	cout << "\tThe LAGO Project, lago@lagoproject.org" << endl;
	cout << endl;
	cout << "\tDPR Lab. 2012" << endl;
	cout << "\tX. Bertou and H. Asorey, asoreyh@cab.cnea.gov.ar" << endl;
	cout << endl;
	cout << "\tRaw data analysis suite for the LAGO Project (L0 -> L1)" << endl; 
	cout << "Usage: " << prog << " [modifiers] <options> <values> <raw_file>[.dat[.bz2]]" << endl;
	cout << endl;
	cout << "\tIf 'raw_file' does not exist, switch to STDIN and use <raw_file> as" << endl;
    cout << "\tbasename for the output files." << endl;
	cout << endl;
	cout << "\tOptions and values:"<< endl;
	cout << "\t-t <a/p> \tproduces .tim time difference histogram file,"  << endl;
	cout << "\t         \tfiltering by a/p if a value is given (default = no)" << endl;
	cout << "\t-n <time>\tproduces .rte scaler file with total number of pulses" << endl;
	cout << "\t         \tper second. It is enabled by default when -l option" << endl;
	cout << "\t         \tis used." << endl;
	cout << "\t         \tIf time (in seconds) is given, it also produces a .flx" << endl;
	cout << "\t         \tflux file with averaged rates over <time> seconds." << endl;
	cout << "\t-v <chn> \tchecks if channel (1-3) is working (bits and baseline)" << endl;
	cout << "\t-N <mask>\tindicate what channels have negative undershoots" << endl;
	cout << "\t         \tUse channel mask 1-7." << endl;
	cout << "\t-u <tr i>\tImpose an offline trigger level for each channel" << endl;
	cout << "\t         \tDefault value: " << trg_default << " ADC)" << endl;
	cout << "\t-l <ch> <t_i>\tdefines the " << SCL_LEVELS << " thresholds t_i for the old" << endl;
	cout << "\t         \tlago-like scalers analysis on channel <ch>."<< endl;
	cout << "\t         \tFor example: -l 1 5 15 30 50, defines subchannels" << endl;
	cout << "\t         \tthresholds for ch 1" << endl;
	cout << "\t         \tDefault thresholds (using -l <ch>): ";
	for (int i=0; i<SCL_LEVELS; i++)
		cout << scl_default[i] << " ";
	cout << endl;
	cout << endl;
	cout << "\tModifiers (note they are case sensitive!):" << endl;
	cout << "\t-h\tprints help and exits" << endl;
	cout << "\t-c\tproduces the .cal calibration file" << endl;
	cout << "\t-s\tproduces the .sol solar physics file" << endl;
	cout << "\t-r\tproduces the .raw " << raw_limit << " second raw file copy" << endl;
	cout << "\t-m\tproduces the .mon monitoring file" << endl;
	cout << "\t-f\tforce analysis for older data versions than " << CODEVERSION << endl;
	cout << "\t-z\tproduces bzip2 compressed files (solar and scaler data)"<<endl;
	cout << "\t-g\tproduces .scl scaler data file (old lago analysis), but" << endl;
	cout << "\t  \tabsolute thresholds (should be used with -l)" << endl;
	cout << "\t-a\tproduces the .all.bz2 compressed file containing charge," << endl;
	cout << "\t  \tpeak, dc and dt, rise time and fall time for all pulses" << endl;
	cout << "\t-i\tOnly include pulses that trigger each channel to fill" << endl; 
	cout << "\t  \tthe calibration histograms of each channel (by default," << endl; 
	cout << "\t  \tall pulses are included." << endl;
    cout << "\t-?\tprints this help and exits" << endl << endl;
	cout << endl;
	exit(1);
}

int main (int argc, char *argv[])
{
	LagoFile Input;

	char nfi[256];
	char *ifiname=NULL;
	int NbPulses=0;
	int scl_ch_set=0;

	// setting scaler levels (same for all channels)
	for (int i = 0; i < CHANNELS; i++)
		for (int j = 0; j < SCL_LEVELS; j++)
			scl_level[i][j] =  scl_default[j];

	// setting scaler counters to 0
	for (unsigned int i = 0; i < CHANNELS; i++)
		for (unsigned int j = 0; j < SCL_LEVELS; j++)
			for (unsigned int l = 0; l < SCL_LEVELS; l++)
				scl_scalers[i][j][l]=0;
	for (int i = 0; i < CHANNELS; i++)
		scl_flux[i] = 0;
	for (int i = 0; i < CHANNELS; i++)
		trg_level[i] = trg_default;

	LagoGeneric Data;
	LagoEvent *Pulse;
	Pulse=(LagoEvent*)malloc(MAXPULSEPERSEC*sizeof(LagoEvent));

	for (int i=0; i<MAXPULSEPERSEC; i++)
		Pulse[i].Init();

	for (int i=1; i<argc; i++) {
		char *tmparg=argv[i];
		if (*tmparg=='-')
			switch (*(tmparg+1)) {
			case 'f':
				force++;
				Input.SetForce(force);
				break;
			case 'c':
				ical=1;
				break;
			case 's':
				isol=1;
				break;
			case 'z':
				izip=1;
				break;
			case 't':
				itim=1;
				if (atof(argv[i+1])) { // false if not a/p value was given
					i++;
					tim_map = atof(argv[i]);
				}
				break;
			case 'r':
				iraw=1;
				break;
			case 'a':
				iall=1;
				break;
			case 'u':
				itrg=1;
				if (atoi(argv[i+1])) { // false if not trigger value was given. Use default levels
					for (int j = 0; j < CHANNELS; j++) {
						i++;
						if (atoi(argv[i])) {
							trg_level[j] = atoi(argv[i]);
						}
						else { // user didn't set the four levels (mandatory): try again.
							cerr << " Error: you must set one trigger level per channel" << endl;
							Usage(argv[0]);
							break;
						}
					}
				}
				break;
			case 'm':
				imon=1;
				break;
			case 'n':
				irte=1;
				if (atof(argv[i+1])) {
					i++;
					flx_default = atof(argv[i]);
					iflx = 1;
				}
				break;
			case 'g':
				isclg=1;
				break;
			case 'i':
				icaltrg=1;
				break;
			case 'l':
				iscl=1;
				irte=1;
				if (atoi(argv[i+1])) { // false if not levels present. Use default levels
					i++;
					scl_ch_set = atoi(argv[i]);
					if (scl_ch_set > 0 && scl_ch_set <= CHANNELS) {
						scl_ch_set -= 1;
						for (int j = 0; j < SCL_LEVELS; j++) {
							i++;
							if (atoi(argv[i])) {
								scl_level[scl_ch_set][j] = atoi(argv[i]);
							}
							else { // user didn't set the four levels (mandatory): try again.
								cerr << " Error: you must set the four levels for channel " << scl_ch_set + 1 << endl << endl;
								Usage(argv[0]);
								break;
							}
						}
					}
					else {
						cerr << " Error: channel number not valid " << endl << endl;
						Usage(argv[0]);
						break;
					}
				}
				break;
			case 'N':
				i++;
				for (int j=0; j<CHANNELS; j++)
					ineg[j] = (atoi(argv[i]) & (1<<j));
				if (atoi(argv[i])<1 || atoi(argv[i])>7) {
					cerr << " Error: channel mask not valid " << endl << endl;
					Usage(argv[0]);
				}
				break;
			case 'v':
				i++;
				ivalid=atoi(argv[i]);
				if (ivalid<1 || ivalid>CHANNELS) {
					cerr << "Error: channel should be within 1-3" << endl;
					return 1;
				}
				break;
			case 'h':
			default:
				Usage(argv[0]);
				break;
			}
		else ifiname=argv[i];
	}
	if (!ifiname) {
		cerr << endl << "Error: Missing filename" << endl << endl;
		Usage(argv[0]);
	}
	snprintf(nfi,256,"%s",ifiname);
	if (isclg && !iscl) {
		cerr << endl << "Error: Try using -g -l" << endl << endl;
		Usage(argv[0]);
	}

	Input.Open(nfi);

	for (int j = 0; j < CHANNELS; j++)
		trg_level[j] -= BASELINE;

	char *ifile2;
	char ifile[256];
	ifile2=ifiname;
	if (strrchr(ifile2,'/')!=NULL) {
		ifile2=strrchr(ifile2,'/')+1;
	}
	snprintf(ifile, 256,"l1_%s_%s",CODEVERSION,ifile2);
	if (strrchr(ifile,'.')!=NULL) {
		if (strncmp(strrchr(ifile,'.'),".bz2",4)==0) { // remove .bz2 if present
			*(strrchr(ifile,'.'))='\0';
		}
	}
	if (strrchr(ifile,'.')!=NULL) {
		if (strncmp(strrchr(ifile,'.'),".dat",4)==0) { // remove .dat if present
			*(strrchr(ifile,'.'))='\0';
		}
	}

	for (int i=0; i<CHANNELS; i++) {
		for (int j=0; j<1024; j++)
			Peak[i][j]=Peak_minute[i][j]=Base[i][j]=0;
		for (int j=0; j<4096; j++)
			ECharge[i][j]=Charge[i][j]=Charge_minute[i][j]=0;
		for (int j=0; j<MAXTIMEINVECTOR; j++)
			Time[i][j]=0;
	}

	if (ivalid)
		fprintf(stderr,"Validating channel %d\n",ivalid);

	if (ical) {
		snprintf(nfi,256,"%s.cal",ifile);
		cal.open(nfi);
	}

	if (itim) {
		snprintf(nfi,256,"%s.tim",ifile);
		tim.open(nfi);
	}

	if (iraw) {
		snprintf(nfi,256,"%s.raw",ifile);
		raw.open(nfi);
	}

	if (iall) {
		snprintf(nfi,256,"bzip2 -9z > %s.all.bz2",ifile);
		if ((all = popen(nfi,"w"))==NULL) {
			fprintf(stderr,"Failed to open compressed all particles file. Abort.\n");
			exit(1);
		}
	}

	if (imon) {
		snprintf(nfi,256,"%s.mon",ifile);
		if ((mon = fopen(nfi,"w"))==NULL) {
			fprintf(stderr,"Failed to open monitoring file. Abort.\n");
			exit(1);
		}
	}
	if (irte) {
		snprintf(nfi,256,"%s.rte",ifile);
		if ((rte = fopen(nfi,"w"))==NULL) {
			fprintf(stderr,"Failed to open rate file. Abort.\n");
			exit(1);
		}
	}

	if (iflx) {
		snprintf(nfi,256,"%s.flx",ifile);
		if ((flx = fopen(nfi,"w"))==NULL) {
			fprintf(stderr,"Failed to open flux file. Abort.\n");
			exit(1);
		}
	}

	if (isol) {
		if (izip) {
			snprintf(nfi,256,"bzip2 -9z > %s.sol.bz2",ifile);
			if ((sol = popen(nfi,"w"))==NULL) {
				fprintf(stderr,"Failed to open compressed solar file. Abort.\n");
				exit(1);
			}
		}
		else {
			snprintf(nfi,256,"%s.sol",ifile);
			if ((sol = fopen(nfi,"w"))==NULL) {
				fprintf(stderr,"Failed to open scaler file. Abort.\n");
				exit(1);
			}
		}
	}

	if (iscl) {
		if (izip) {
			snprintf(nfi,256,"bzip2 -9z > %s.scl.bz2",ifile);
			if ((scl = popen(nfi,"w"))==NULL) {
				fprintf(stderr,"Failed to open compressed scaler file. Abort.\n");
				exit(1);
			}
		}
		else {
			snprintf(nfi,256,"%s.scl",ifile);
			if ((scl = fopen(nfi,"w"))==NULL) {
				fprintf(stderr,"Failed to open scaler file. Abort.\n");
				exit(1);
			}
		}
	}

	if (ical) {
		cal << "# # # p 1 cal " << PROJECT << " " << CODEVERSION << endl;
		cal << "# # L1 level file (processed raw data, use at your own risk or contact lago@lagoproject.org)" << endl;
		cal << "# # This is a file containing the charge and peak calibration histograms" << endl;
		cal << "# # Format is ch1 ch2 ch3 pk1 pk2 pk3" << endl;
		if (itrg)
			cal << "# # An offline trigger of " << trg_level[0] << " " << trg_level[1] << " " << trg_level[2] << " ADC above baseline has been used for each channel respectively." << endl;
		if (icaltrg)
			cal << "# # For each channel we only used triggered pulses at this particular channel (-i option)" << endl;
	}

	if (itim) {
		tim << "# # # p 1 tim " << PROJECT << " " << CODEVERSION << endl;
		tim << "# # L1 level file (processed raw data, use at your own risk or contact lago@lagoproject.org)" << endl;
		tim << "# # This is a file containing the time difference histogram" << endl;
		tim << "# # Format is #time_difference(ns) #count for each channel" << endl;
		if (tim_map)
			tim << "# # Pulses were discarded if (a/p < " << tim_map << ")"<< endl;
		if (itrg)
			tim << "# # An offline trigger of " << trg_level[0] << " " << trg_level[1] << " " << trg_level[2] << " ADC above baseline has been used for each channel respectively." << endl;
	}

	if (iraw) {
		raw << "# # # p 1 raw " << PROJECT << " " << CODEVERSION << endl;
		raw << "# # L1 level file (processed raw data, use at your own risk or contact lago@lagoproject.org)" << endl;
		raw << "# # This is a RAW file containing the first 10 seconds of data" << endl;
	}

	if (isol) {
		fprintf(sol, "# # # p 1 sol %s %s\n", PROJECT, CODEVERSION);
		fprintf(sol, "# # L1 level file (processed raw data, use at your own risk or contact lago@lagoproject.org)\n");
		fprintf(sol, "# # This is a Solar data file.\n");
		fprintf(sol, "# # These are one minute charge and peak histograms, with monitoring information\n");
		fprintf(sol, "# # Format is # q/p second frequency temperature pressure\n");
		fprintf(sol, "# # (q for charge and p for peak)\n");
		fprintf(sol, "# # followed by 0/1 0/1/2 second frequency temperature pressure and 1024 or 4096 values\n");
		fprintf(sol, "# # where 0/1 stands for charge (0) or peak (1) and 0/1/2 is the channel\n");
		if (itrg)
			fprintf(sol, "# # An offline trigger of %d %d %d ADC above baseline has been used for each channel respectively.\n", trg_level[0], trg_level[1], trg_level[2]);
	}

	if (iall) {
		fprintf(all, "# # # p 1 all %s %s\n", PROJECT, CODEVERSION);
		fprintf(all, "# # L1 level file (processed raw data, use at your own risk or contact lago@lagoproject.org)\n");
		fprintf(all, "# # This is a file containing all processed data\n");
		fprintf(all, "# # Format is # second frequency temperature pressure\n");
		fprintf(all, "# #   time_to_prev_pulse counts_to_prev_pulse channels_triggered_mask ch_ch1 pk_ch1 ch_ch2 pk_ch2 ch_ch3 pk_ch3 risetime_ch1 falltime_ch1 risetime_ch2 falltime_ch2 risetime_ch3 falltime_ch3 (channel <i> is printed only if it triggered)\n");
		fprintf(all, "# # note: time_to_prev_pulse is in clock cycle, and reset to 0 every second\n");
		if (itrg)
			fprintf(all, "# # An offline trigger of %d %d %d ADC above baseline has been used for each channel respectively.\n", trg_level[0], trg_level[1], trg_level[2]);
	}

	if (imon) {
		fprintf(mon, "# # # p 1 all %s %s\n", PROJECT, CODEVERSION);
		fprintf(mon, "# # L1 level file (processed raw data, use at your own risk or contact lago@lagoproject.org)\n");
		fprintf(mon, "# # This is a monitoring file.\n");
		fprintf(mon, "# # Format is second frequency temperature pressure average_baseline_chN dev_baseline_chN\n");  
	}

	if (iscl) {
		fprintf(scl, "# # # p 1 scl %s %s\n", PROJECT, CODEVERSION);
		fprintf(scl, "# # L1 level file (processed raw data, use at your own risk or contact lago@lagoproject.org)\n");
		fprintf(scl, "# # This is the old-lago-like scaler file.\n");
		fprintf(scl, "# # Format is:\n");
		fprintf(scl, "# # # s second frequency temperature pressure Total_Rate Rate_per_channel:_(%d)_columns\n",CHANNELS);
		fprintf(scl, "# # followed by %d lines (%.1e s) containing %d counting levels for each channel\n", SCL_BINS, SCL_TIME, SCL_LEVELS);
		fprintf(scl, "# # These are the thresholds used:\n");
		for (int i = 0; i < CHANNELS; i++) {
			fprintf(scl, "# # Channel %d : ", i+1);
			for (int j = 0; j < SCL_LEVELS; j++)
				fprintf(scl, "l%d = %d ", j+1, scl_level[i][j]);
			fprintf(scl,"\n");
		}
		fprintf(scl,"# # Scalers thresholds are ");
		if (isclg)
			fprintf(scl, "considered as absolute levels.\n");
		else
			fprintf(scl, "compared with the baseline level for each channel (%d adc)\n", BASELINE);
		if (itrg)
			fprintf(scl, "# # An offline trigger of %d %d %d ADC above baseline has been used for each channel respectively.\n", trg_level[0], trg_level[1], trg_level[2]);
	}
	if (irte) {
		fprintf(rte, "# # # p 1 rte %s %s\n", PROJECT, CODEVERSION);
		fprintf(rte, "# # L1 level file (processed raw data, use at your own risk or contact lago@lagoproject.org)\n");
		fprintf(rte, "# # This is the rate (pulse per second) file.\n");
		fprintf(rte, "# # Format is:\n");
		fprintf(rte, "# # second temperature pressure Total_Rate Rate_per_channel:_(%d)_columns\n",CHANNELS);
	}
	if (iflx) {
		fprintf(flx, "# # # p 1 flx %s %s\n", PROJECT, CODEVERSION);
		fprintf(flx, "# # L1 level file (processed raw data, use at your own risk or contact lago@lagoproject.org)\n");
		fprintf(flx, "# # This is the flux file.\n");
		fprintf(flx, "# # Format is:\n");
		fprintf(flx, "# # second_at_half_interval avg_temp avg_press avg_rate avg_rate_per_channel sigma_temp sigma_press sigma_rate sigma_rate_per_channel\n");
		fprintf(flx, "# # Average interval used: %d seconds\n", flx_default);
	}

	/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	//
	// File reading and processing
	//
	cerr << "Reading file" << endl;
	while(NbPulses!=-1) {
		NbPulses=Input.ReadOneSecond(&Data,Pulse,MAXPULSEPERSEC);
		if (NbPulses>0) TreatSecond(&Data,Pulse,NbPulses);
	}


	// Filling calibration histograms
	for (int j=0; j<4096; j++) {
		for (int i=0; i<CHANNELS; i++)  {
			cal << Charge[i][j] << " ";
		}
		if (j<1024) for (int i=0; i<CHANNELS; i++)
				cal << Peak[i][j] << " ";
		cal << endl;
	}

	//
	// Channel validation code
	//
	if (ivalid) {
		printf("Validation of channel %d:\n",ivalid);
		// validity checks...
		ivalid=ivalid-1; // 0 to 2
		// baseline
		double b=0,b2=0,n=0;
		for (int j=0; j<1024; j++) {
			b+=Base[ivalid][j]*1.*j;
			b2+=Base[ivalid][j]*1.*j*j;
			n+=Base[ivalid][j];
		}
		if (n<1000) {
			printf("@ Problem: not enough events on this channel (%d)\n",int(n));
		}
		else {
			//printf("%lf %lf %lf\n",b,b2,n);
			float m=b*1./n;
			float v=sqrt(b2*1./n-m*m);
			int bb=0;
			for (int j=0; j<1024; j++)
				if (j<m-5*v || j>m+5*v) bb+=Base[ivalid][j];
			printf("  Baseline: %lf +- %lf, %lf%% at more than 5 sigmas\n",m,v,bb*100./n);
			if (m-v>BASELINE || m+v<BASELINE) printf("@ Problem: Baseline is not within 1 sigma of %d\n",BASELINE);
			if (v>2) printf("@ Problem: Baseline fluctuation is high, more than 2 ADC\n");
			// peak histogram
			int BitHist[10][10];
			for (int i=0; i<10; i++)
				for (int j=0; j<10; j++) BitHist[i][j]=0;
			int tv=0;
			for (int j=BASELINE; j<1024; j++) {
				tv+=Peak[ivalid][j-BASELINE];
				for (int i=0; i<10; i++)
					for (int k=0; k<10; k++)
						if (int(pow(2,i))==(j&(int(pow(2,i)))))
							if (int(pow(2,k))==(j&(int(pow(2,k)))))
								BitHist[i][k]+=Peak[ivalid][j-BASELINE];
			}
			for (int i=0; i<10; i++)
				for (int j=i; j<10; j++)
					if (BitHist[i][j]==0) printf("@ Problem: bit combination %d %d never happenning\n",i,j);
			//for (int j=0; j<1024; j++)
			//printf("%d %d %d\n",j,Peak[ivalid][j],Base[ivalid][j]);
		}
	}
	//
	// Time difference vector
	//
	if (itim) {
		for (int i=0; i<MAXTIMEINVECTOR; i++) {
			tim << i * 25 << " ";
			for (int j=0; j<CHANNELS; j++)
				tim << Time[j][i] << " ";
			tim << endl;
		}
	}
	//
	// Closing files
	//
	if (ical)
		cal.close();
	if (itim)
		tim.close();
	if (iall)
		pclose(all);
	if (iraw)
		raw.close();
	if (imon)
		fclose(mon);
	if (irte)
		fclose(rte);
	if (iflx)
		fclose(flx);
	if (izip) {
		if (isol)
			pclose(sol);
		if (iscl)
			pclose(scl);
	}
	else {
		if (isol)
			fclose(sol);
		if (iscl)
			fclose(scl);
	}
}
