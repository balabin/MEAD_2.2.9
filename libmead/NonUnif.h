// This is -*- C++ -*-
#ifndef _NonUnif_h
#define _NonUnif_h 1

/*
$Id: NonUnif.h,v 1.3 1995/04/17 22:15:53 bashford Exp $
*/

/*
    This source code file is part of the MEAD (Macroscopic
    Electrostatics with Atomic Detail) package of objects and
    programs, copyright (C) 1990 by Donald Bashford of the Department
    Molecular Biology of The Scripps Research Institute.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 1, or (at your option)
    any later version.

    The GNU General Public License should be in a file called COPYING.GNU.
    Some comments about copying and distribution by D. Bashford should
    be in a file called COPYING.DB

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

    Donald Bashford can be contacted by electronic mail by the address

       bashford@scripps.edu 

    or by paper mail at

       Department of Molecular Biology
       The Scripps Research Institute
       10666 North Torrey Pines Road
       La Jolla, California  92037
*/


  struct NonUnif { // For doing calculation on non-uniform epsilon points.
    float *phiptr;
    float *ptr1, fac1;    // No, I did not forget a * here
    float *ptr2, fac2;
    float *ptr3, fac3;
    float *ptr4, fac4;
    float *ptr5, fac5;
    float *ptr6, fac6;
  };

#endif
