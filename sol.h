/*  sol.h  --  code oriented to perform Space Weather related analysis from 1-minute histograms
 *
 *  Copyright (C) 2012-TODAY The LAGO Project, http://lagoproject.org, lago-pi@lagoproject.org
 *
 *  Original authors: Hernán Asorey
 *  e-mail: asoreyh@cab.cnea.gov.ar  (asoreyh@gmail.com)
 *  Laboratorio de Detección de Partículas y Radiación (LabDPR)
 *  Centro Atómico Bariloche - San Carlos de Bariloche, Argentina */

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
