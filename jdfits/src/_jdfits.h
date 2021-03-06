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
extern char *_jdfits_allocate_bytes_of_type (int type, unsigned int nelements);
extern JDFits_Bintable_Field_Type *_jdfits_bintable_find_column (JDFits_Bintable_Type *, char *name, unsigned int *);

extern JDFits_Keyword_Type *_jdfits_find_keyword (JDFits_Header_Type *h, char *name);

extern char *_jdfits_skip_whitespace (char *);

#ifdef HAVE_FSEEKO
# define FSEEK(a,b,c) fseeko(a,b,c)
# define FTELL(a) ftello(a)
#else
# define FSEEK(a,b,c) fseek(a,b,c)
# define FTELL(a) ftell(a)
#endif
