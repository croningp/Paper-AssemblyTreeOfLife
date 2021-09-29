/* It borrows external InChI API to
   obtain the InChI string of a molecule or fragment.
 */
#pragma once
#include "MOB_Def.h"
#include "MOB.h"


strType getInchi(MOL_BOND& mol);