/*  sol.cc  -- code oriented to perform Space Weather related analysis from 1-minute histograms
 *
 *  Copyright (C) 2012-TODAY The LAGO Project, http://lagoproject.org, lago-pi@lagoproject.org
 *
 *  Original authors: Hernán Asorey
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
#include "sol.h"

/* Functions */

int Open(char *nfi)
{
	char tmpc[256];
	if (strncmp(nfi+strlen(nfi)-4,".bz2",4)==0) {
		fprintf(stderr,"## WARNING ## file '%s' ending by bz2, reading it via bzip2\n",nfi);
		snprintf(tmpc,256,"bzip2 -d -c %s",nfi);
		fi=popen(tmpc,"r");
	}
	else if((fi=fopen(nfi,"r")) == NULL) {
		fprintf(stderr,"## WARNING ## Can't open file '%s' - Switching to SDTIN\n",nfi);
		fi=stdin;
	}
	return 1;
}

inline double sign(double x)
{
	return (x < 0.0 ? -1.0 : 1.0);
}

inline double log10(double x)
{
	return log(x)/log(10.);
}

double itpl_pos(double *pos_avg, double *chr_der, int i)
{
	double y0=chr_der[i];
	double x0=pos_avg[i];
	double y1=chr_der[i+1];
	double x1=pos_avg[i+1];
	return (x0 - y0 * (x1 - x0)/(y1-y0));
}

double itpl_chr(double *pos_avg, double *chr_avg, int i, double x)
{
	double y0=chr_avg[i];
	double x0=pos_avg[i];
	double y1=chr_avg[i+1];
	double x1=pos_avg[i+1];
	return (y0 + (y1 - y0) * (x - x0) / (x1 - x0));
}

double fitHisto(int xi, int xf, int* h, double* r)
{
	double y1, x1, x2, x3, x4, xy, x2y, x, y, v1, v2, v;
	v = v1 = v2 = x = y = y1 = x1 = x2 = x3 = x4 = xy = x2y = 0.0;
	int n = 0;
	for (int i = xi; i <= xf; i++) {
		x = i;
		y = (double) h[i];
		x1 += x;
		y1 += y;
		xy += x * y;
		x *= i;
		x2 += x;
		x2y += x * y;
		x *= i;
		x3 += x;
		x *= i;
		x4 += x;
		n++;
	}
	v1 = x1 * x2 * x2y - n * x2y * x3 - x2 * x2 * xy + n * x4 * xy + x2 * x3 * y1 - x1 * x4 * y1;
	v2 = x1 * x1 * x2y - n * x2 * x2y + n * x3 * xy + x2 * x2 * y1 - x1 * (x2 * xy + x3 * y1);
	v = v1 / (2.0 * v2); // note v1 = -a1 and v2 = a2, so we have v = -b / 2a !
	*r = 0.;
	return v;
}

void Usage(char *prog)
{
	cout << endl << PROJECT << " " << VERSION << endl;
	cout << endl << "This is the Space Weather dedicated analysis (L1 -> L2SW)" << endl<< endl;
	cout << "  Analyze .sol output files to identify solar modulation phenomena." << endl;
	cout << "  Usage: " << prog << " [flags] <filename>.sol.bz2" << endl;
	cout << "    If 'file.sol.bz2' does not exits, switch to STDIN and use <filename> as" << endl;
	cout << "    basename for output file." << endl;
	cout << "    Options:"<<endl;
	cout << "      -c <ch>            : channel to be analyzed [1,2,3]" << endl;
	cout << "      -r <Qs> <Qe>       : calculates the 1-min integrated flux in the signal range Qs<S<Qe, units in ADCq" << endl;
	cout << "                         : if not Qe is given, use default (Qe=MaxCharge="<< maxcharge << ")" << endl;
	cout << "      -b <width> <start> : calculates the 1-min integrated flux for the histogram in bands of <width> ADCq" << endl;
	cout << "                           first band starts at <start> (default start=1 ADCq)" << endl;
	cout << "    Flags and Modifiers:"<<endl;
	cout << "      -p                 : do the analysis of peak histograms (by default, it use only charge histograms)" << endl;
	cout << "      -s                 : calculates the flux in counts/s instead of counts/min" << endl;
	cout << "      -n                 : disable automatic detection of features and flux calculations (default: enabled)" << endl;
	cout << "      -g                 : if auto is enabled, produced the .pos file to analyze goodness of automatic fit. If not, enables automatic detection" << endl;
	cout << "      -z                 : files output are compressed" << endl;
	cout << "      -v                 : verbose mode" << endl;
	cout << "      -h                 : prints help and exits" << endl << endl;
	exit(1);
}

int main (int argc, char *argv[])
{
	char nfi[256]; // filename container
	char *ifiname=NULL; // input

	/* reading arguments from command line */
	for (int i=1; i<argc; i++) {
		char *tmparg=argv[i];
		if (*tmparg=='-') {
			switch (*(tmparg+1)) {
			case 'p':
				ipeak=1;
				break;
			case 'c':
				i++;
				channel=atoi(argv[i]); // channel to be analyzed
				break;
			case 'b':
				iband=1;
				i++;
				width=atoi(argv[i]); // band width
				if (atoi(argv[i+1])) {// false if not start given
					i++;
					start=atoi(argv[i]); // starting point
					if (!start)
						start=1;
				}
				break;
			case 'r':
				irange=1;
				i++;
				qs = atoi(argv[i]) - 1; // band width
				if (atoi(argv[i+1])) {// false if not Qend given
					i++;
					qe = atof(argv[i]) - 1; // starting point
				}
				break;
			case 'z':
				icompress=1;
				break;
			case 's':
				iflux=1;
				fluxtime=60.;
				break;
			case 'n':
				iauto=0;
				break;
			case 'g':
				ipos=1;
				break;
			case 'v':
				iverbose=1;
				break;
			case 'h':
				Usage(argv[0]);
				break;
			default:
				cerr << "\nUnkown option: -" << (*(tmparg+1)) << endl << endl;
				Usage(argv[0]);
				break;
			}
		}
		else
			ifiname=argv[i];
	}
	/* end reading arguments from command line */

	/* Some security checks */
	if (!ifiname) {
		cerr << endl << "\n## ERROR ##  Missing filename" << endl << endl;
		Usage(argv[0]);
	}
	if (iband && irange) {
		cerr << endl << "\n## ERROR ## You can't setup band and range modes at the same time" << endl << endl;
		Usage(argv[0]);
	}
	if (!iband && !irange) {
		cerr << endl << "\n## ERROR ## You have to select band or range modes to proceed" << endl << endl;
		Usage(argv[0]);
	}
	if (irange) {
		if (qs >= qe) {
			fprintf(stderr,"\n## ERROR ## Qs(%d) has to be smaller than Qe(%d). Abort.\n",qs,qe);
			Usage(argv[0]);
		}
	}
	if (!channel) {
		fprintf(stderr,"\n## ERROR ## You have to select the channel to be analyzed. Abort.\n");
		Usage(argv[0]);
	}
	if (!iauto && ipos) {
		fprintf(stderr,"\n## WARNING ## You selected both no automatic detection and goodness of fit. Enabling automatic features detection\n");
		iauto = 1;
	}


	/* input */
	snprintf(nfi,256,"%s",ifiname);
	Open(nfi);
	if (iverbose)
		fprintf(stderr, "## STATUS ## Analyzing %s\n", nfi);

	/* preparing outputs */
	char *ifile2;
	char ifile[256], ifile0[256];
	ifile2=ifiname;
	if (strrchr(ifile2,'/')!=NULL)
		ifile2=strrchr(ifile2,'/')+1;      // remove dirs if present
	snprintf(ifile, 256,"%s",ifile2);
	if (strrchr(ifile,'.')!=NULL)
		if (strncmp(strrchr(ifile,'.'),".bz2",4)==0) // remove .bz2 if present
			*(strrchr(ifile,'.'))='\0';
	if (strrchr(ifile,'.')!=NULL)
		if (strncmp(strrchr(ifile,'.'),".sol",4)==0) // remove .sol if present
			*(strrchr(ifile,'.'))='\0';
	strcpy(ifile0,ifile);
	if (strncmp(ifile,"l1_",3)==0) {
		*(strchr(ifile,'1'))='2';  // change l1 by l2
	}
	else {
		cerr << endl << "\n## ERROR ## You should use a L1 sol file (starting with l1_)" << endl << endl;
		Usage(argv[0]);
	}

	/* end of preparing outputs */

	// log file
	snprintf(nfi,256,"%s.log", ifile);
	if ((lof = fopen(nfi,"w"))==NULL) {
		fprintf(stderr,"\n## ERROR ## Failed to open log file. Abort.\n");
		exit(1);
	}
	fprintf(lof, "# # # L2 log %s %s\n", PROJECT, VERSION);
	fprintf(lof, "# # L2 level file (lago@lagoproject.org): fit log for the sol automatic procedure\n");
	fprintf(lof, "## STATUS ## Analyzing %s\n", nfi);

	/* preparing for work */
	if (ipeak) {
		type='1';
		end=1023;
		if (qe > end) {
			fprintf(stderr,"## WARNING ## for peak analysis max. value is 1023. We will integrate up to 1023\n");
			fprintf(lof,"## WARNING ## for peak analysis max. value is 1023. We will integrate up to 1023\n");
			qe=end;
		}
	}

	if (qe > maxcharge) {
		fprintf(stderr,"## WARNING ## qe can't be greater than maxcharge=%d. Changing Qe to %d\n", maxcharge, maxcharge);
		fprintf(lof,"## WARNING ## qe can't be greater than maxcharge=%d. Changing Qe to %d\n", maxcharge, maxcharge);
		qe=end;
	}

	/* determine number of bands */
	int nbands=0;
	double nb=0.;
	if (iband) {
		nb=(1.0*(end-start)/width);
		if (nb-int(nb))
			nb+=1.;
		nbands=int(nb);
	}

	/* open outputs and headers */
	if (icompress) {  // open output file (.flx) via fopen or popen for compressed output files
		snprintf(nfi,256,"bzip2 -9z > %s.flx.bz2", ifile);
		if ((out = popen(nfi,"w"))==NULL) {
			fprintf(stderr,"\n## ERROR ## Failed to open compressed output file. Abort.\n");
			exit(1);
		}
		if (iauto) {
			snprintf(nfi,256,"bzip2 -9z > %s.auto.bz2", ifile);
			if ((muo = popen(nfi,"w"))==NULL) {
				fprintf(stderr,"\n## ERROR ## Failed to open compressed automatic flux file. Abort.\n");
				exit(1);
			}
		}
	}
	else {
		snprintf(nfi,256,"%s.flx",ifile);
		if ((out = fopen(nfi,"w"))==NULL) {
			fprintf(stderr,"\n## ERROR ## Failed to open output file. Abort.\n");
			exit(1);
		}
		if (iauto) {
			snprintf(nfi,256,"%s.auto", ifile);
			if ((muo = fopen(nfi,"w"))==NULL) {
				fprintf(stderr,"\n## ERROR ## Failed to open automatic flux file. Abort.\n");
				exit(1);
			}
		}
	}

	// headers
	fprintf(out, "# # # L2 flx %s %s\n", PROJECT, VERSION);
	fprintf(out, "# # L2 level file (lago@lagoproject.org): flux of signals as function of time\n");
	fprintf(out, "# # Analysis is done by integrating the whole channel %d charge histogram in\n", channel);
	if (iband) {
		fprintf(out, "# # %d bands of %d ADCq starting from %d ADCq and up to %d ADCq\n", nbands, width, start+1, end+1);
		fprintf(out, "# # %d columns format is:\n",nbands+4);
		if (iflux)
			fprintf(out, "# # utc pressure temperature total_flux(counts/s) band001 band002 ... band%03d\n",nbands);
		else
			fprintf(out, "# # utc pressure temperature total_flux(counts/min) band001 band002 ... band%03d\n",nbands);
	}
	if (irange) {
		fprintf(out, "# # the range (%d <= S <= %d) ADCq\n", qs+1, qe+1);
		fprintf(out, "# # %d columns Format is:\n",4);
		if (iflux)
			fprintf(out, "# # utc pressure temperature total_flux(count/s)\n");
		else
			fprintf(out, "# # utc pressure temperature total_flux(count/min)\n");
	}
	if (iauto) {
		fprintf(muo, "# # # L2 flx %s %s\n", PROJECT, VERSION);
		fprintf(muo, "# # L2 level file (lago@lagoproject.org): flux of signals as function of time\n");
		fprintf(muo, "# # Analysis is done by automaticaly detection of calibration histogram features\n");
		fprintf(muo, "# # on channel %d, and integrating the flux in four different bands of deposited energy.\n", channel);
		fprintf(muo, "# # 8 columns format is:\n");
		if (iflux)
			fprintf(muo, "# # utc pressure temperature total_flux(counts/s) lowEd_band muon_band highEd_band HAXB_band\n");
		else
			fprintf(muo, "# # utc pressure temperature total_flux(counts/min) lowEd_band muon_band highEd_band HAXB_band\n");
	}

	/* end open outputs and headers */
	if (iverbose)
		fprintf(stderr, "## STATUS ## Output files created\n");
	fprintf(lof, "## STATUS ## Output files created\n");

	/* and here is where fun begins... */
	long int lines = 0;
	int typen, chn;
	double pressure, temperature;
	long int utc;
	double tmp;
	double flux=0.;
	int histo[maxcharge];
	char line[16384];
	double fluxband[nbands];
	int currentpos=0, pos=0, i=0, is=0, j=0, k=0;

	/*
	 * New analysis done by Y. Perez and H. Asorey
	 * Obtain muon charge and other features from cal histogram
	 * and use them for determining bands of interest
	 * three features:
	 * em = feat[0] = EM Peak (handle with care);
	 * transition = feat[1] = EM-Muon transition
	 * muon = feat[2] = Muon hump
	 */
	int nfeat = 3;
	double feat[nfeat];
	for (i = 0; i < nfeat; i++)
		feat[i] = 0.;
	int nfeatband=4;
	double featband[nfeatband];
	for (i = 0; i < nfeatband; i++)
		featband[i] = 0.;
	double em = 0., transition = 0., muon = 0.;
	int muon_max = 0, em_start = 0, transition_point = 0;
	bool go=true;

	/*
	 * now, we have to read calibration (1h) histogram from ifile.cal, so:
	 * trying to open histogram file
	 */
	if (iauto) {
		snprintf(nfi,256,"%s.cal",ifile0);
		if ((cal = fopen(nfi,"r"))==NULL) {
			fprintf(stderr,"\n## ERROR ## Failed to open calibration histogram file (%s). Abort.\n",nfi);
			exit(1);
		}
		/* read the histogram */
		if (iverbose)
			fprintf(stderr, "## STATUS ## Reading calibration histogram %s\n",nfi);
		fprintf(lof, "## STATUS ## Reading calibration histogram %s\n",nfi);
		int p[3], c[3], charge[maxcharge];
		int cline=0, ncharge=0;
		while(fgets(line, 16384,cal)) {
			cline++;
			if (line[0] == '#')   // comments? histogram meta-data? for now, just skipping
				continue;
			if (cline < 1028) {
				do { // read
					check=sscanf(line, "%d %d %d %d %d %d\n", &c[0], &c[1], &c[2], &p[0], &p[1], &p[2]);
				}
				while (check != 6);
			}
			else {
				do { // read
					check=sscanf(line, "%d %d %d\n", &c[0], &c[1], &c[2]);
				}
				while (check != 3);
			}
			if (ncharge>maxcharge)
				fprintf(stderr,"## WARNING ## Too long charge histogram. It has %d lines (expected lines = %d)\n", ncharge, maxcharge);
			charge[ncharge] = c[channel];
			ncharge++;
		}
		fclose(cal);
		/*
		 * now, charge[] has the charge histogram for the analyzed channel
		 * average it!, but, graded
		 */

		/* Analysis parameters */
		const int max_em = 200;  /* estimated value for EM peak influence */
		const int max_muon = 2500; /* max expected charge for muon */
		const int low_energy_step = 5; /* average step for low energy */
		const int high_energy_step = 50; /* average step for high energy */
		// bands
		// TODO need a parser for config file...
		// in the meantime
		const double diameter = 1580.; // nahuelito, in mm
		const double height = 1560.; // nahuelito, in mm
		// END TODO
		const double muon_limit = 0.9; // TODO to avoid high energy contamination, should be < 1
		const double muon_geometry = sqrt(1. + (diameter * diameter) / (height * height) );
		const double em_factor = 5.;  // TODO start 1st band at em_peak*em_factor
		const double transition_factor = 1.; // TODO end 1st band at transition_point*transition_factor
		int sum = 0, cnt_avg=0;
		double chr_avg[maxcharge], pos_avg[maxcharge];
		for (i=0; i<maxcharge; i++)
			chr_avg[i] = pos_avg[i] = 0.;
		int step=0;
		/* low energy */
		step = low_energy_step;
		for (i=0; i<max_em; i+=step) {
			sum = 0;
			for (j=0; j<step; j++)
				sum += charge[i+j];
			chr_avg[cnt_avg] = 1. * sum / step;
			pos_avg[cnt_avg] = i + 1. * step / 2.;
			cnt_avg++;
		}
		/* high energy */
		step = high_energy_step;
		for (i=max_em; i<maxcharge; i+=step) {
			sum = 0;
			for (j=0; j<step; j++)
				sum += charge[i+j];
			chr_avg[cnt_avg] = 1. * sum / step;
			pos_avg[cnt_avg] = i + 1. * step / 2.;
			cnt_avg++;
		}
		/* Histogram are averaged. Now, derive it! */
		cnt_avg-=2; //avoid last average point because it's usally too noise
		double chr_der[cnt_avg];
		for (i=0; i<cnt_avg; i++) //taking logs to magnify the effect (thanks Mauricio Suarez)
			chr_der[i] = ((log10(chr_avg[i+1])-log10(chr_avg[i]))/(pos_avg[i+1]-pos_avg[i]));
		// looking for 0s on histogram derivative...
		double zeros[cnt_avg], zeros_pos[cnt_avg];
		int cnt_zeros=0;
		for (i=0; i<cnt_avg; i++)
			zeros_pos[i]=zeros[i]=0.;
		for (i=0; i<cnt_avg-1; i++) {
			if (!sign(chr_der[i])) { // bingo... too good to be true
				zeros_pos[cnt_zeros]=(double) i;
				zeros[cnt_zeros]=chr_avg[i]; // charge
				cnt_zeros++;
			}
			else if (sign(chr_der[i+1]) != sign(chr_der[i])) { // Bolzano's theorem!!
				zeros_pos[cnt_zeros]=itpl_pos(pos_avg, chr_der,i);
				zeros[cnt_zeros]=itpl_chr(pos_avg, chr_avg, i, zeros_pos[cnt_zeros]);
				cnt_zeros++;
			}
		}
		// check for NaN
		for (i=0; i < cnt_zeros ; i++)
			if (zeros_pos[i] != zeros_pos[i]) {//check for NaN!
				zeros_pos[i] = 0;
				zeros[i]=0;
			}
		// Looking for real features
		int max_pos[cnt_zeros];
		int min_pos[cnt_zeros];
		int cmin=0, cmax=0;
		for (i=0; i < cnt_zeros ; i++) {
			max_pos[i] = 0;
			min_pos[i] = 0;
		}
		for (i=1; i < cnt_zeros ; i++) {
			// Is this a local maximum? type=2
			if (zeros[i] > zeros[i-1] && zeros[i] > zeros[i+1] && zeros_pos[i]<max_muon) {
				max_pos[cmax] = i;
				if (!cmax && zeros_pos[i] < max_em)
					cmax++;
				if (zeros_pos[i] > max_em)
					cmax++;
			}
			// Is this a local minimum? type=1
			if (zeros[i] < zeros[i-1] && zeros[i] < zeros[i+1] && zeros_pos[i]<max_muon) {
				min_pos[cmin] = i;
				if (zeros_pos[i] > max_em)
					cmin++;
			}
		}
		if (cmax == 1 && cmin == 1) { // one is missing
			if (min_pos[0] < max_pos[0]) { // looks like we found transition_point and muon
				max_pos[1] = max_pos[0];  // moving the maximum
				int max = 0; // looking for the missing maximum
				for (i=0; i <= min_pos[1]; i++) {
					if (zeros[i] > max) {
						max = zeros[i];
						max_pos[0] = i;
					}
				}
			}
		}
		/* finally, histogram seeds */
		double zeros_fit[nfeat];
		zeros_fit[0] = zeros_pos[max_pos[0]];
		zeros_fit[1] = zeros_pos[min_pos[0]];
		zeros_fit[2] = zeros_pos[max_pos[1]];
		/* So, we have the histogram seeds, let's fit them */
		double tmpk=0.;
		if (ipos) {
			snprintf(nfi,256,"%s.pos",ifile);
			if ((pof = fopen(nfi,"w"))==NULL) {
				fprintf(stderr,"\n## ERROR ## Failed to open pos file. Abort.\n");
				exit(1);
			}
			fprintf(pof, "# # # L2 pos %s %s\n", PROJECT, VERSION);
			fprintf(pof, "# # L2 level file (lago@lagoproject.org): file to check goodness of fit for automatic detection\n");
			fprintf(pof, "# # 4 columns format is:\n");
			fprintf(pof, "# # Feature_Index Seed_for_fit Feature analyzed_file\n");
		}
		for (i=0; i < nfeat ; i++) {
			feat[i] = fitHisto(int(zeros_fit[i]*0.9), int(zeros_fit[i]*1.3), charge, &tmpk);
			if (ipos) {
				do {
					check=fprintf(pof, "%d %.3f %.3f %s\n", i, zeros_fit[i], feat[i], ifile);
				}
				while (check<=0);
			}
		}
		/*
		 * and done. We found the EM peak, the transition_point to hump, and the muon hump
		 * of the calibration histogram, in charge units =D
		 */
		// features
		em = feat[0];
		transition = feat[1];
		muon = feat[2];
		// band limits
		em_start = int(em * em_factor + 0.5);
		transition_point = int(transition * transition_factor + 0.5);
		muon_max = int(muon * muon_limit * muon_geometry + 0.5);

		/* now, we need some security checks to see if everything was fine during fits, or introduce some alerts */
		/* of course, they should be positive */
		if (!(em > 0. && transition > 0. && muon > 0.)) {
			fprintf(lof, "## WARNING ## Problems during fit: Some feature(s) have negative value(s): em=%.3f, transition=%.3f, muon=%.3f\n", em, transition, muon);
			go=false;
		}
		// and they need to have reasonable values
		if (go) {
			if (!(em < max_muon && transition < max_muon && muon < max_muon)) {
				fprintf(lof, "## WARNING ## Problems during fit: Some feature(s) have very large value(s): em=%.3f, transition=%.3f, muon=%.3f\n", em, transition, muon);
				go=false;
			}
		}
		// and the order should be em < transition < muon (should be okay due to previous checks, but...)
		if (go) {
			if (em > transition) {
				fprintf(lof, "## WARNING ## Problems during fit: em=%.3f > transition=%.3f\n", em, transition);
				go=false;
			}
			if (em > muon) {
				fprintf(lof, "## WARNING ## Problems during fit: em=%.3f > muon=%.3f\n", em, muon);
				go=false;
			}
			if (transition > muon) {
				fprintf(lof, "## WARNING ## Problems during fit: transition=%.3f > muon=%.3f\n", transition, muon);
				go=false;
			}
		}
		// and the order for band limits should be reasonables
		if (go) {
			if (!(em_start < transition_point && transition_point < muon_max && em_start < muon_max)) {
				fprintf(lof, "## WARNING ## Problems with this histogram: em_start=%d, transition_point=%d, muon_max=%d\n", em_start, transition_point, muon_max);
				go=false;
			}
		}
		// if everything was fine... we can go. If not, we produce a file called L2_*.err to inform something was wrong...
		if (go) {
			if (iverbose)
				fprintf(stderr, "## STATUS ## All test passed. Fit procedure exit without warnings. We can proceed with automatic band analysis.\n");
			fprintf(lof, "## STATUS ## All test passed. Fit procedure exit without warnings. We can proceed with automatic bands analysis.\n");
			fprintf(lof, "## STATUS ## Features found at: em=%.3f, transition=%.3f, muon=%.3f\n", em, transition, muon);
			fprintf(lof, "## STATUS ## Bands limits: em_start=%d, transition_point=%d, muon_max=%d\n", em_start, transition_point, muon_max);
			fprintf(muo, "# # Features found at: em=%.3f, transition=%.3f, muon=%.3f\n", em, transition, muon);
			fprintf(muo, "# # Bands limits: em_start=%d, transition_point=%d, muon_max=%d\n", em_start, transition_point, muon_max);
		}
		else {
			if (iverbose)
				fprintf(stderr, "## STATUS ## Some test failed for this file. Please check manually.\n");
			fprintf(lof, "## STATUS ## Some test failed for this file. Please check manually.\n");
			snprintf(nfi,256,"touch %s.err",ifile);
			system(nfi);
			if (iverbose)
				fprintf(stderr, "## STATUS ## Error indicator file %s.err created\n", ifile);
			fprintf(lof, "## STATUS ## Error indicator file %s.err created\n", ifile);
		}
	}
	/*
	 * END of new automatic analysis block by Y. Perez and H. Asorey
	 * sat ene 17 16:28:46 COT 2015
	 */

	/* analyze the .sol histograms */

	/* reading .sol file */
	while(fgets(line, 16384, fi)) {
		if (!(++lines % 1000))
			cerr << lines << "\r";
		if (line[0] == '#')   // comments? meta-data? for now, just skipping
			continue;
		else if (line[0] == type) {
			if ((line[2]-48) == channel) { // 48 is the ASCII code for 0, then '2' - 48 = 2
				currentpos=pos=i=0;
				do { // read
					check=sscanf(line, "%d %d %ld %lf %lf %lf %n", &typen, &chn, &utc, &tmp, &temperature, &pressure, &currentpos);
				}
				while (check != 6);
				if (!utc)
					continue;
				for (k=0; k<maxcharge; k++)
					histo[k]=0.;
				// reading histogram...
				while (1==sscanf(line+currentpos, "%d %n", &histo[i], &pos)) {
					currentpos+=pos;
					i++;
				}
				// done
				if (irange) {
					flux=0.;
					for (i=qs; i<=qe; i++)
						flux+=histo[i];
					do {
						check=fprintf(out, "%ld %.2f %.2f %.2f\n", utc, pressure, temperature, flux/fluxtime);
					}
					while (check<=0);
				}
				if (iband) {
					flux=0.;
					for (k=0; k<nbands; k++)
						fluxband[k]=0.;
					is = start;
					for (i=0; i<nbands; i++) {
						for (j=0; j<width; j++) {
							if((is+j)>=maxcharge)
								break;
							fluxband[i] += histo[is+j];
						}
						flux+=fluxband[i];
						is+=width;
					}
					do {
						check=fprintf(out, "%ld %.2f %.2f %.2f", utc, pressure, temperature, flux/fluxtime);
					}
					while (check<=0);
					for (i=0; i<nbands; i++) {
						do {
							check=fprintf(out, " %.2f", 1. * fluxband[i] / fluxtime);
						}
						while (check<=0);
					}
					fprintf(out, "\n");
				}
				if (iauto && go) {
					flux=0.;
					for (k=0; k<nfeatband; k++)
						featband[k]=0.;
					// EM Band, index 0, from max_em to transition
					for (i=em_start; i<transition_point; i++)
						featband[0] += histo[i];
					// Muon Band, index 1, from transition to muon_max
					for (i=transition_point; i<muon_max; i++)
						featband[1] += histo[i];
					// High energy band, index 2, from muon_max to maxcharge
					for (i=muon_max; i<maxcharge; i++)
						featband[2] += histo[i];
					// New band, index 3, from muon to maxcharge
					for (i=int(muon+0.5); i<maxcharge; i++)
						featband[3] += histo[i];
					for (i=em_start; i<maxcharge; i++)
						flux += histo[i];
					do {
						check=fprintf(muo, "%ld %.2f %.2f %.2f", utc, pressure, temperature, flux/fluxtime);
					}
					while (check<=0);
					for (i=0; i<nfeatband; i++) {
						do {
							check=fprintf(muo, " %.2f", 1. * featband[i] / fluxtime);
						}
						while (check<=0);
					}
					fprintf(muo, "\n");
				}
			}
		}
	}
	// tell something on the automatic analysis file (because it will be empty!)
	if (iauto && !go) {
		fprintf(muo, "## STATUS ## For some reason, the automatic detection of histogram features failed.\n");
		fprintf(muo, "## STATUS ## Please check %s.log file to see the gory details\n", ifile);
		fprintf(muo, "## STATUS ## Automatic band analysis was not performed\n");
	}
	// say goodbye
	if (iverbose)
		fprintf(stderr, "## STATUS ## Done. Closing files\n");
	fprintf(lof, "## STATUS ## Done. Closing files\n");
	if (icompress) {
		pclose(out);
		if (iauto)
			pclose(muo);
	}
	else {
		fclose(out);
		if (iauto)
			pclose(muo);
	}
	if (ipos)
		fclose(pof);
	fclose(lof);
}
