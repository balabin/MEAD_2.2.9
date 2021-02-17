// -*- C++ -*-
#ifndef _meadsolvate_h
#define _meadsolvate_h

#include <vector>
using namespace std;

#include "MEAD/globals.h"
#include "MEAD/PhysCond.h"
#include "MEAD/AtomSet.h"
#include "MEAD/FDGridLevel.h"
#include "MEAD/ChargeDist.h"
#include "MEAD/AtomChargeSet.h"
#include "MEAD/DielectricEnvironment.h"
#include "MEAD/DielByAtoms.h"
#include "MEAD/ElectrolyteEnvironment.h"
#include "MEAD/ElectrolyteByAtoms.h"
#include "MEAD/FinDiffMethod.h"
#include "MEAD/ElstatPot.h"
#include "MEAD/FDGridLevel.h"

void meadsolvate(
		 int     *nqueries,
		 double  solrad,
		 double  epssol,
		 double  ionic_str,
		 int     doReactionField,
		 AtomSet atom_set,
		 double  *spacing,
		 int     *dim,
		 vector<CenteringStyle> censtl,
		 string  fieldpoint_filename,
		 string  reacfield_filename
		 );

#endif
