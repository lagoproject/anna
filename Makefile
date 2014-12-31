#   Makefile   -- 
#   Copyright (C) 2012-TODAY The LAGO Project, http://lagoproject.org, lago-pi@lagoproject.org
#   Original authors: Hernán Asorey
#   e-mail: asoreyh@cab.cnea.gov.ar  (asoreyh@gmail.com)
#   Laboratorio de Detección de Partículas y Radiación (LabDPR)
#   Centro Atómico Bariloche - San Carlos de Bariloche, Argentina 
#
#   LICENSE GPLv3
#   This program is free software: you can redistribute it and/or modify 
#   it under the terms of the GNU General Public License as published by 
#   the Free Software Foundation, either version 3 of the License, or 
#   (at your option) any later version.
# 
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
# 
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>. 
#
#   -*- coding: utf8 -*-
#   try to preserve encoding  
CC = g++

TARGETS = raw grb sol
CFLAGS = -Wall

all: $(TARGETS)

raw: raw.cc lago_file.h lago_data.h
	$(CC) -o $@ $< $(CFLAGS)

grb: grb.cc lago_file.h lago_data.h
	$(CC) -o $@ $< $(CFLAGS)

sol: sol.cc lago_file.h lago_data.h
	$(CC) -o $@ $< $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(TARGETS)

