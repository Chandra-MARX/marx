/*
    This file is part of the MIT PFILE Parameter Library

    Copyright (C) 2002 Massachusetts Institute of Technology

    This software was developed by the MIT Center for Space Research
    under contract SV1-61010 from the Smithsonian Institution.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
#ifndef _MIT_PARAMETER_H_
#define _MIT_PARAMETER_H_
#include "pfile.h"

/* The following routines are compatability routines for working with SAO code. */
typedef Param_File_Type *paramfile;
extern paramfile c_paramopen (char *, char **, int, char *);
extern void c_paramclose (paramfile);
extern short c_pgets (paramfile, char *);
extern int c_pgeti (paramfile, char *);
extern long c_pgetl (paramfile, char *);
extern float c_pgetf (paramfile, char *);
extern double c_pgetd (paramfile, char *);
extern char *c_pgetstr (paramfile, char *, char *, int);
extern int c_paccess (paramfile, char *);
extern char *c_paramgetpath (paramfile);
extern int c_pgetb (paramfile, char *);

#define paramopen	c_paramopen
#define paramclose	c_paramclose
#define pgets		c_pgets
#define pgeti		c_pgeti
#define pgetb		c_pgetb
#define pgetl		c_pgetl
#define pgetf		c_pgetf
#define pgetd		c_pgetd
#define pgetstr		c_pgetstr
#define paccess		c_paccess
#define paramgetpath	c_paramgetpath

extern int ParamInfo (paramfile pfile, char *name, char *mode, char *type, char *value,
		      char *min, char *max, char *prompt);

#endif
