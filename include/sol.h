/*
################################################################################
# File:        sol.h
# Description: header file for the 'sol' LAGO ANNA application
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
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

//global variables
int irange=0, ipos=0, iband=0, ipeak=0, icompress=0, iflux = 0, iauto=1, iverbose=0;
int check = 0;
const int max_charge=4096;
int avg_time = 10;
double fluxtime=1.;
double diameter = 0., height = 0.;
int channel=0, start=0, width=0, qs=0, end_hst=max_charge, qe=max_charge;
char type='0'; //charge
int Open(char *nfi);
inline double sign(double x);
inline double log10(double x);
double itpl_pos(double *pos_avg, double *chr_der, int i);
double itpl_chr(double *pos_avg, double *chr_der, int i, double x);
double fitHisto(int xi, int xf, double* h, double* r);
void Usage(char *prog);
FILE *fi, *out, *cal, *muo, *pof, *lof;
