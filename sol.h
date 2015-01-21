/*  sol.h  --  code oriented to perform Space Weather related analysis from 1-minute histograms
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
const int maxcharge=4096;
double fluxtime=1.;
int channel=0, start=0, width=0, qs=0, end=maxcharge, qe=maxcharge;
char type='0'; //charge
int Open(char *nfi);
inline double sign(double x);
inline double log10(double x);
double itpl_pos(double *pos_avg, double *chr_der, int i);
double itpl_chr(double *pos_avg, double *chr_der, int i, double x);
double fitHisto(int xi, int xf, double* h, double* r);
void Usage(char *prog);
FILE *fi, *out, *cal, *muo, *pof, *lof;
