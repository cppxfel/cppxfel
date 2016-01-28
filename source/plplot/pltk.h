// $Id: pltk.h 11295 2010-11-01 22:19:45Z airwin $
//
// Maurice LeBrun
// 28-Apr-95
//
// Declarations for PLplot/TK utility routines.
//

#ifndef __PLTK_H__
#define __PLTK_H__

#include "pltcl.h"
#include <tk.h>

#ifdef __cplusplus
extern "C" {
#endif

// Generic initialization routine for extended wish'es.

#ifdef BUILD_Plplot
#undef TCL_STORAGE_CLASS
#define TCL_STORAGE_CLASS    DLLEXPORT
#endif // BUILD_Vfs

EXTERN PLDLLIMPEXP_TCLTK int
Pltk_Init( Tcl_Interp *interp );

// Main program for wish, with PLplot extensions.

EXTERN PLDLLIMPEXP_TCLTK int
pltkMain( int argc, const char **argv, char *RcFileName,
          int ( *AppInit )( Tcl_Interp *interp ) );

#ifdef __cplusplus
}
#endif

#endif  // __PLTK_H__
