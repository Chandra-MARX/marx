/*
    Copyright (C) 2002 MIT Center For Space Research

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
#include <stdio.h>
#include <time.h>

#include "jdfits.h"

/* This algorithm was borrowed from IDL */
double jdfits_time_t_to_mjd (time_t t)
{
   struct tm *tms;
   double day;
   int year, month;
   int a, b;
   
   tms = gmtime (&t);
   
   day = (double) tms->tm_mday +
     (tms->tm_hour + (tms->tm_min + tms->tm_sec/60.0)/60.0)/24.0;
   
   if (1582 == (year = tms->tm_year + 1900))
     {
	jdfits_error ("jdfits_time_t_to_mjd: year 1582 not covered.");
	return 0.0;
     }
	
   month = tms->tm_mon + 1;
   if (month <= 2)
     {
	month += 12;
	year -= 1;
     }
   
   a = year / 100;
   
   if (year < 1582)
     b = 0;
   else
     b = 2 - a + a/4;

   return (int)(year * 0.25) + 365.0 * (year - 1860) 
     + (int) (30.6001 * (month + 1)) + b + day - 105.5;
}

