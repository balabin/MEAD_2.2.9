// -*- C++ -*-
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

#include "meadsolinprot.h"

#include "MEAD/globals.h"
#include "MEAD/AtomSet.h"
#include "MEAD/FDGridLevel.h"


extern "C" 
void domeadsolinprot_(
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
		    char    *diel2id,
		    int     *diel2id_len,
		    int     *icycle,
		    int     *blablevel,
		    int     *nqueries,
		    double  *epsin1,
		    double  *epsin2,
		    double  *epsext,
		    double  *solrad,
		    double  *ionic_str,
		    int     *doReactionField,
		    int     *doProteinField
){

  string runname = runid;
  string diel2_name=diel2id;

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
    
  /* Set up the AtomSet for the quantum part */
  string qpt_name = runname + "_" + cyclestr;
  AtomSet qpt(qpt_name);
  //qpt.read();
  cout << "Reading coords for " << *atcount << " atoms from pqr file" << endl;
  int inserted = 0;
  blab3 << "atname_len " << *atname_len << endl;
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
    at.chainid = "QPT";
    AtomID key(at);
    if(qpt.contains(key)){
      cout << "Duplicate atom in quantum portion " << key << endl;
      ::error("AtomSet::insert: attempt to insert a duplicate atom \n");
    }
    else{
      qpt.insert(at);
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
      :: error("INPUT FAILURE in domeadsolinprot: ",cstyle_str,
	       " unrecognized CenteringStyle");
  }

  string fieldpoint_filename;
  string reacfield_filename;
  string protfield_filename;

  if (*doReactionField) {
    
    stringstream ss;
    string cyclestr;
    
    ss << *icycle;
    ss >> cyclestr;

    fieldpoint_filename = runname + "_" + cyclestr + "_meadquery" + ".dat" ;
    reacfield_filename = runname + "_" + cyclestr + ".mrf";
  }

  if (*doProteinField) {

    stringstream ss;
    string cyclestr;
    
    ss << *icycle;
    ss >> cyclestr;

    fieldpoint_filename = runname + "_" + cyclestr + "_meadquery" + ".dat" ;
    protfield_filename =  runname + "_" + cyclestr + ".mpf";
    cout << "fieldpoint filename "<< fieldpoint_filename << endl;
    cout << "protfield filename "<< protfield_filename << endl;
  }
  
  meadsolinprot(nqueries, *solrad, *epsin1, *epsin2, *ionic_str, 
		*doReactionField, qpt, spacing, dim, censtl, fieldpoint_filename, 
		reacfield_filename, *doProteinField, diel2_name, 
		protfield_filename);
}
