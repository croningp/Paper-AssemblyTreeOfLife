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


#ifndef __CTL_DATA_H__
#define __CTL_DATA_H__
#include "e_ichisize.h"
/***********************************************/
#define STR_ERR_LEN 256

typedef struct tagStructData {
    unsigned long ulStructTime;
    int           nErrorCode;
    int           nErrorType;
    int           nStructReadError;
    int           bChiralFlag;
    char          pStrErrStruct[STR_ERR_LEN];
    long          fPtrStart;
    long          fPtrEnd;
    /* debugging info */
#if( bRELEASE_VERSION == 0 )
    int           bExtract;
#endif

} STRUCT_DATA;
/***********************************************/

#define MAX_NUM_PATHS 4


/* SDF treatment */
#define MAX_SDF_VALUE        255 /* max lenght of the SDFile data value */
#define MAX_SDF_HEADER        64 /* max length of the SDFile data header */

/***********************************************/

/* bCalcInChIHash values */
typedef enum tagInChIHashCalc 
{   
    INCHIHASH_NONE=0, 
    INCHIHASH_KEY=1, 
    INCHIHASH_KEY_XTRA1=2, 
    INCHIHASH_KEY_XTRA2=3, 
    INCHIHASH_KEY_XTRA1_XTRA2=4 
} 
INCHI_HASH_CALC;

typedef struct tagInputParms {
    char            szSdfDataHeader[MAX_SDF_HEADER+1];
    char           *pSdfLabel;
    char           *pSdfValue;
    long            lSdfId;
    long            lMolfileNumber;
#ifndef INCHI_ANSI_ONLY
    DRAW_PARMS      dp;
    PER_DRAW_PARMS  pdp;
    TBL_DRAW_PARMS  tdp;
#endif
/*
  -- Files --
  ip->path[0] => Input
  ip->path[1] => Output (INChI)
  ip->path[2] => Log
  ip->path[3] => Problem structures
  ip->path[4] => Errors file (ACD)

*/
    const char     *path[MAX_NUM_PATHS];
    int             num_paths;
    long            first_struct_number;
    long            last_struct_number;
    INPUT_TYPE      nInputType;
    INCHI_MODE       nMode;
    int             bAbcNumbers;
    /*int             bXml;*/
    int             bINChIOutputOptions; /* !(ip->bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT) */
    int             bCtPredecessors;
    int             bXmlStarted;
    int             bDisplayEachComponentINChI;

    long            msec_MaxTime;   /* was ulMaxTime; max time to run ProsessOneStructure */
    long            msec_LeftTime;

    /* unsigned long   ulMaxTime; */
    unsigned long   ulDisplTime;
    int             bDisplay;
    int             bDisplayIfRestoreWarnings; /* InChI->Struct debug */
    int             bMergeAllInputStructures;
    int             bSaveWarningStructsAsProblem;
    int             bSaveAllGoodStructsAsProblem;
    int             bGetSdfileId;
    int             bGetMolfileNumber;  /* read molfile number from the name line like "Structure #22" */
    int             bCompareComponents; /* see flags CMP_COMPONENTS, etc. */
    int             bDisplayCompositeResults;
    int             bDoNotAddH;
    int             bNoStructLabels;
    int             bChiralFlag;
    int             bAllowEmptyStructure;
    /*^^^ */
    int             bCalcInChIHash;
    int             bFixNonUniformDraw; /* correct non-uniformly drawn oxoanions and amidinium cations. */
    /*^^^ */
    INCHI_MODE       bTautFlags;
    INCHI_MODE       bTautFlagsDone;

#if( READ_INCHI_STRING == 1 )
    int             bReadInChIOptions;
#endif

/* post v.1 features */
#if( UNDERIVATIZE == 1 )
    int             bUnderivatize;
#endif
#if( RING2CHAIN == 1 )
    int             bRing2Chain;
#endif
#if ( RING2CHAIN == 1 || UNDERIVATIZE == 1 )
    int             bIngnoreUnchanged;
#endif

} INPUT_PARMS;


/*************************** INChI mode *******************************/
/* ip->nMode */
#define REQ_MODE_BASIC              0x000001    /* B    */
#define REQ_MODE_TAUT               0x000002    /* T    */
#define REQ_MODE_ISO                0x000004    /* I    */
#define REQ_MODE_NON_ISO            0x000008    /* NI   */
#define REQ_MODE_STEREO             0x000010    /* S    */
#define REQ_MODE_ISO_STEREO         0x000020    /* IS   */
#define REQ_MODE_NOEQ_STEREO        0x000040    /* SS   */
#define REQ_MODE_REDNDNT_STEREO     0x000080    /* RS   */
#define REQ_MODE_NO_ALT_SBONDS      0x000100    /* NASB */
/* new 10-10-2003 */
#define REQ_MODE_RELATIVE_STEREO    0x000200    /* REL All Relative Stereo */
#define REQ_MODE_RACEMIC_STEREO     0x000400    /* RAC All Racemic Stereo */
#define REQ_MODE_SC_IGN_ALL_UU      0x000800    /* IAUSC Ignore stereocenters if All Undef/Unknown */
#define REQ_MODE_SB_IGN_ALL_UU      0x001000    /* IAUSC Ignore stereobonds if All Undef/Unknown */
#define REQ_MODE_CHIR_FLG_STEREO    0x002000    /* SUCF  If Chiral flag then Abs otherwise Rel stereo */
/* end of 10-10-2003 */
/*^^^ 2009-12-05 */
#define REQ_MODE_DIFF_UU_STEREO     0x004000    /* SLUUD Make labels for unknown and undefined stereo different */
/*^^^ 2009-12-05 */

#define REQ_MODE_MIN_SB_RING_MASK   0x0F0000    /* RSB  */
#define REQ_MODE_MIN_SB_RING_SHFT      16

#define REQ_MODE_DEFAULT  (REQ_MODE_BASIC | REQ_MODE_TAUT | REQ_MODE_ISO | REQ_MODE_NON_ISO | REQ_MODE_STEREO)

/*********** compare components flags **********************************/
/* ip->bCompareComponents */
#define CMP_COMPONENTS              0x0001     /* perform compare components */
#define CMP_COMPONENTS_NONISO       0x0002     /* ignore isotopic */
#define CMP_COMPONENTS_NONTAUT      0x0004     /* compare non-tautomeric */

/****************** chemical identifier member definitions *************/
/* ip->bINChIOutputOptions */
#define INCHI_OUT_NO_AUX_INFO           0x0001   /* do not output Aux Info */
#define INCHI_OUT_SHORT_AUX_INFO        0x0002   /* output short version of Aux Info */
#define INCHI_OUT_ONLY_AUX_INFO         0x0004   /* output only Aux Info */
#define INCHI_OUT_EMBED_REC             0x0008   /* embed reconnected INChI into disconnected INChI */
#define INCHI_OUT_SDFILE_ONLY           0x0010   /* save input data in a Molfile instead of creating INChI */
#define INCHI_OUT_XML                   0x0020   /* output xml INChI */
#define INCHI_OUT_PLAIN_TEXT            0x0040   /* output plain text INChI */
#define INCHI_OUT_PLAIN_TEXT_COMMENTS   0x0080   /* output plain text annotation */
#define INCHI_OUT_XML_TEXT_COMMENTS     0x0100   /* output xml text annotation */
#define INCHI_OUT_WINCHI_WINDOW         0x0200   /* output into wINChI text window */
#define INCHI_OUT_TABBED_OUTPUT         0x0400   /* tab-delimited (only for plain text) */
#define INCHI_OUT_SDFILE_ATOMS_DT       0x0800   /* SDfile output H isotopes as D and T */
#define INCHI_OUT_SDFILE_SPLIT          0x1000   /* Split SDfile into components */
/*^^^ */
#define INCHI_OUT_FIX_TRANSPOSITION_CHARGE_BUG 0x2000   
                                    /* used to accomodate FIX_TRANSPOSITION_CHARGE_BUG */
/*^^^ */
#define INCHI_OUT_STDINCHI 0x4000
#define INCHI_OUT_SAVEOPT  0x8000

#define FLAG_INP_AT_CHIRAL         1
#define FLAG_INP_AT_NONCHIRAL      2
#define FLAG_SET_INP_AT_CHIRAL     4
#define FLAG_SET_INP_AT_NONCHIRAL  8

/* unknown/undefined stereo - constants  */
#define AB_PARITY_UNKN   3  /* 3 => user marked as unknown parity */
#define AB_PARITY_UNDF   4  /* 4 => parity cannot be defined because of symmetry or not well defined geometry */

#endif /* __CTL_DATA_H__ */
