#include "meadsolvate.h"

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
){

  float epsin = 1.0;
  float epsvac = 1.0;
  PhysCond::set_epsext(epssol);
  PhysCond::set_solrad(solrad);
  PhysCond::set_ionicstr(ionic_str);
  FDGridLevel::set_ZeroForOutOfRange(true);

  cout << "Starting MEAD solvate calculation " << endl;
  cout << "with interior dielectric constant = " << epsin << endl;
  cout << "using the following physical conditions:"<< endl;
  PhysCond::print();
  cout << "and vacuum dielectric constant = " << epsvac << endl;
  if (blab3pt == &cout) cout << "Blab level set to " << 3 << endl;
  else if (blab2pt == &cout) cout << "Blab level set to " << 2 << endl;
  else if (blab1pt == &cout) cout << "Blab level set to " << 1 << endl;
  else cout << "No blab level set (so no blabbing) " << endl;

  ChargeDist rho(new AtomChargeSet(atom_set));

  DielectricEnvironment eps(new TwoValueDielectricByAtoms(atom_set, epsin));

  ElectrolyteEnvironment ely(new ElectrolyteByAtoms(atom_set));

  FinDiffMethod fdm;
  for (int ispec = 0; ispec < censtl.size(); ++ispec){
    fdm.add_level(*(dim+ispec), *(spacing+ispec), censtl[ispec]) ;
  }

  Coord interesting (0,0,0);
  fdm.resolve(atom_set.geom_cent(), interesting);
  cout << "Using finite difference method with lattice levels:" << endl;
  cout << fdm;

  ElstatPot phi(fdm, eps, rho, ely);
  phi.solve();
  float prod_sol = phi * rho;
  cout << "prod_sol = " << prod_sol << endl;

  PhysCond::set_epsext(epsvac);
  PhysCond::set_ionicstr(0.0);

  ElectrolyteEnvironment elyvac;  // No electrolyte is the default
  DielectricEnvironment vac_eps(new TwoValueDielectricByAtoms(atom_set, epsin));
  ElstatPot vac_phi(fdm, vac_eps, rho, elyvac);
  vac_phi.solve();
  float prod_vac = vac_phi * rho;
  cout << "prod_vac = " << prod_vac << endl;
  float solvation_energy = (prod_sol - prod_vac) / 2 * PhysCond::get_econv();
  cout << "\n\nSOLVATION ENERGY = " << solvation_energy
    << "\n(probably in kcal/mole)" << endl;

  if (doReactionField) {

    ifstream fpt (fieldpoint_filename.c_str());
    if (!fpt.good())
      ::error("domeadsolvate", "main: failed to open for reading, ",
	      fieldpoint_filename.c_str());
    ofstream rf(reacfield_filename.c_str());
    if (!rf.good())
      ::error("domeadsolvate", "main: failed to open for writing, ",
	      reacfield_filename.c_str());

    *nqueries = 0;
    while (fpt.good()) {
      Coord c;
      fpt >> c;
      if (fpt.eof()) break;
      if (!fpt.good())
	::error("domeadsolvate", "main: input failure while reading ",
		fieldpoint_filename.c_str());
      rf << phi.value(c) - vac_phi.value(c) << "\n";
      if (!rf.good())
	::error("domeadsolvate", "main: output failure while writing ",
		reacfield_filename.c_str());
      ++(*nqueries);
    }
    rf.close();
    
  }
}
