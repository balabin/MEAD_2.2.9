// -*- C++ -*-
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

#include "meadsolvate.h"

#include "MEAD/globals.h"
#include "MEAD/AtomSet.h"
#include "MEAD/FDGridLevel.h"


extern "C" 
void domeadsolvate_(
		    int     *atcount,
		    double  *radius,
		    double  *charge,
		    double  *coordx,
		    double  *coordy,
		    double  *coordz,
		    int     *atname_len,
		    char    *atname,
		    int     *nb_specs,
		    char    *cstyle,
		    int     *style_str_len,
		    double  *center,
		    double  *spacing,
		    int     *dim,
		    char    *runid,
		    int     *runid_len,
		    int     *icycle,
		    int     *blablevel,
		    int     *nqueries,
		    double  *solrad,
		    double  *epssol,
		    double  *ionic_str,
		    int     *doReactionField
){

  string runname = runid;

  if (*blablevel >= 1){
    blab1pt = &cout;
  }
  if (*blablevel >= 2){
    blab2pt = &cout;
  }
  if (*blablevel >= 3){
    blab3pt = &cout;
  }

  string cyclestr;
  stringstream ss;

  ss << *icycle;
  ss >> cyclestr;
    
  string pqr_runname = runname + "_" + cyclestr;
  AtomSet a(pqr_runname);
  //a.read();
  //blab3 << "Reading coords from pqr file" << endl;
  int inserted = 0;
  blab3 << "AtomSet::insert: starting" << endl;
  for (int i=0; i<*atcount; ++i){
    Atom at;
    at.resnum = 1;
    at.atname = (atname+i*(*atname_len));
    at.resname = "NONE";
    at.coord.x = *(coordx+i);
    at.coord.y = *(coordy+i);
    at.coord.z = *(coordz+i);
    at.charge = *(charge+i);
    at.rad = *(radius+i);
    at.chainid = '\0';
    at.rad = *(radius+i);
    at.chainid = '\0';
    AtomID key(at);
    if(a.contains(key)){
      ::error("AtomSet::insert: attempt to insert a duplicate atom\n");
    }
    else{
      a.insert(at);
      inserted++;
    }
  }
  blab1 << "AtomSet::insert: finished with " << inserted << " atoms" << endl;

  vector<CenteringStyle> censtl(*nb_specs);
  for (int ispec = 0; ispec < censtl.size(); ++ispec){
    string cstyle_str = cstyle+(ispec*(*style_str_len));
    if (cstyle_str == "ON_ORIGIN"){
      censtl[ispec]=ON_ORIGIN;
    }
    else if (cstyle_str == "ON_GEOM_CENT"){
      censtl[ispec]=ON_GEOM_CENT;
    }
    else
      :: error("INPUT FAILURE in domeadsolvate: ",cstyle_str,
	       " unrecognized CenteringStyle");
  }

  string fieldpoint_filename;
  string reacfield_filename;

  if (*doReactionField) {
    
    stringstream ss;
    string cyclestr;
    
    ss << *icycle;
    ss >> cyclestr;

    fieldpoint_filename = runname + "_" + cyclestr + "_meadquery" + ".dat" ;
    reacfield_filename = runname + "_" + cyclestr + ".mrf";

  }

  meadsolvate(nqueries, *solrad, *epssol, *ionic_str, *doReactionField, a, 
	      spacing, dim, censtl, fieldpoint_filename, reacfield_filename);

}
