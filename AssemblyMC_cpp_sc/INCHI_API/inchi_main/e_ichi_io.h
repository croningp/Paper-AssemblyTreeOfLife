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


#ifndef __E_ICHI_IO_H__
#define __E_ICHI_IO_H__


#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif


#ifndef INCHI_LIBRARY


#define inchi_ios_print e_inchi_ios_print

#define GetMaxPrintfLength e_GetMaxPrintfLength
#define inchi_fgetsTab e_inchi_fgetsTab
#define inchi_vfprintf e_inchi_vfprintf
#define inchi_ios_str_getc e_chi_ios_str_getc
#define inchi_ios_str_gets e_chi_ios_str_gets
#define inchi_ios_str_getsTab e_chi_ios_str_getsTab

#define inchi_fprintf e_inchi_fprintf

#define inchi_fgetsLf e_inchi_fgetsLf
#define inchi_fgetsLfTab e_inchi_fgetsLfTab
#define inchi_ios_print_nodisplay e_inchi_ios_print_nodisplay
#define inchi_ios_eprint e_inchi_ios_eprint
#define inchi_ios_init e_inchi_ios_init
#define inchi_ios_flush e_inchi_ios_flush
#define inchi_ios_flush2 e_inchi_ios_flush2
#define inchi_ios_close e_inchi_ios_close
#define inchi_ios_reset e_inchi_ios_reset
#define inchi_ios_gets e_inchi_ios_gets
#define inchi_ios_getsTab e_inchi_ios_getsTab
#define inchi_ios_getsTab1 e_inchi_ios_getsTab1

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


                                    INCHI_IOSTREAM OPERATIONS 


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/


void inchi_ios_init(INCHI_IOSTREAM *ios, int io_type, FILE *f);
void inchi_ios_flush(INCHI_IOSTREAM *ios);
void inchi_ios_flush2(INCHI_IOSTREAM *ios, FILE *f2);
void inchi_ios_close(INCHI_IOSTREAM *ios);
void inchi_ios_reset(INCHI_IOSTREAM *ios);

int inchi_ios_gets( char *szLine, int len, INCHI_IOSTREAM *ios, int *bTooLongLine );
int inchi_ios_getsTab( char *szLine, int len, INCHI_IOSTREAM *ios, int *bTooLongLine );
int inchi_ios_getsTab1( char *szLine, int len, INCHI_IOSTREAM *ios, int *bTooLongLine );

int inchi_ios_print( INCHI_IOSTREAM *ios, const char* lpszFormat, ... );
int inchi_ios_print_nodisplay( INCHI_IOSTREAM *ios, const char* lpszFormat, ... );

/* Print to string buffer or to file+stderr */
int inchi_ios_eprint( INCHI_IOSTREAM *ios, const char* lpszFormat, ... );


/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


                                    PLAIN FILE OPERATIONS 


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

/* Print to file, echoing to stderr */
int inchi_fprintf( FILE* f, const char* lpszFormat, ... );

char* inchi_fgetsLf( char* line, int line_len, FILE* inp );
int inchi_fgetsLfTab( char *szLine, int len, FILE *f );


#endif



void e_PrintFileName( const char *fmt, FILE *output_file, const char *szFname );
unsigned long e_ulMyGetTickCount( int bStart );
unsigned long e_ulMyTickCountDiff( unsigned long ulTickEnd, unsigned long ulTickStart );

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif /* __E_ICHI_IO_H__ */
