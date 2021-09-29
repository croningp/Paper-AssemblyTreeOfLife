/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.03
 * May 9, 2010
 *
 * Originally developed at NIST
 * Modifications and additions by IUPAC and the InChI Trust
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC);
 * you can redistribute this software and/or modify it under the terms of 
 * the GNU Lesser General Public License as published by the Free Software 
 * Foundation:
 * http://www.opensource.org/licenses/lgpl-2.1.php
 */


#ifndef ___INCHISIZE_H__
#define ___INCHISIZE_H__


typedef unsigned short AT_NUMB;
typedef unsigned short AT_RANK;
#define AT_RANK_MASK   ((AT_RANK)~0)

typedef signed short NUM_H;
#define MAX_ATOMS  1024

#define CHAR_MASK  0xFF

typedef unsigned long INCHI_MODE;

#define LEN_COORD 10
#define NUM_COORD 3
typedef char MOL_COORD[LEN_COORD*NUM_COORD + NUM_COORD-1]; /*copied 30 bytes from MOLfile */

typedef enum tagInputType { INPUT_NONE=0, INPUT_MOLFILE=1, INPUT_SDFILE=2, INPUT_INCHI_XML=3, INPUT_INCHI_PLAIN=4, INPUT_CMLFILE=5, INPUT_INCHI=6, INPUT_MAX } INPUT_TYPE;

/* other types */

#define UCINT  (int)(unsigned char)
#ifndef INCHI_US_CHAR_DEF
typedef signed char   S_CHAR;
typedef unsigned char U_CHAR;
#define INCHI_US_CHAR_DEF
#endif



#endif /* ___INCHISIZE_H__ */

