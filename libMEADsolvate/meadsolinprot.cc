#include "meadsolinprot.h"

void meadsolinprot(
		   int     *nqueries,
		   double  solrad,
		   double  epsin1,
		   double  epsin2,
		   double  ionic_str,
		   int     doReactionField,
		   AtomSet qpt,
		   double  *spacing,
		   int     *dim,
		   vector<CenteringStyle> censtl,
		   string  fieldpoint_filename,
		   string  reacfield_filename,
		   int     doProteinField,
		   string  diel2_name,
		   string  protfield_filename
){
  
  float epsvac = 1.0;
  PhysCond::set_solrad(solrad);
  PhysCond::set_ionicstr(ionic_str);
  FDGridLevel::set_ZeroForOutOfRange(true);

  cout << "Starting MEAD solinprot calculation " << endl;
  cout << "using the following physical conditions:"<< endl;
  PhysCond::print();
  cout << "  Solute interior dielectric, epsin1 = " << epsin1 << endl;
  cout << "  Protein interior dielectric constant, epsin2 = " << epsin2 << endl;
  cout << "and vacuum dielectric constant = " << epsvac << endl;
  if (blab3pt == &cout) cout << "Blab level set to " << 3 << endl;
  else if (blab2pt == &cout) cout << "Blab level set to " << 2 << endl;
  else if (blab1pt == &cout) cout << "Blab level set to " << 1 << endl;
  else cout << "No blab level set (so no blabbing) " << endl;

  ChargeDist rhogen(new AtomChargeSet(qpt));

  blab3 << "Reading protein coords." << endl;
  AtomSet diel2(diel2_name);
  diel2.read();

  ChargeDist rhofeel(new AtomChargeSet(diel2));

  DielectricEnvironment eps(new ThreeValueDielectricByAtoms(qpt, epsin1, diel2, epsin2));

  // Make the union of qpt and diel2 for defining electrolyte boundary
  // (FIXME Don thinks STL provides a function to do this.)
  AtomSet combined(qpt);
  for (AtomSet::const_iterator i = diel2.begin(); i!=diel2.end(); ++i) {
    const Atom& a = i->second;
    AtomID k = i->first;
    if (qpt.contains(k)) {
      cerr << "domeadsolinprot: WARNING:  The atom, \"" << k
	<< "\",\n   occurs in both the solute and protein atom sets.\n"
	  << "   For defining the Elelctrolyte boundary, the protein version\n"
	    << "   will be used."  << endl;
    }
    else {
      combined.insert(a);
    }
  }
  ElectrolyteEnvironment ely(new ElectrolyteByAtoms(AtomSet(combined)));

  FinDiffMethod fdm;
  for (int ispec = 0; ispec < censtl.size(); ++ispec){
    fdm.add_level(*(dim+ispec), *(spacing+ispec), censtl[ispec]) ;
  }

  Coord interesting (0,0,0);
  fdm.resolve(qpt.geom_cent(), interesting);
  cout << "Using finite difference method with lattice levels:" << endl;
  cout << fdm;

  ElstatPot phi(fdm, eps, rhogen, ely);
  phi.solve();
  float prod_sol = phi * rhogen;
  cout << "prod_sol = " << prod_sol << endl;
  float protein_interaction = phi * rhofeel;
  protein_interaction *= PhysCond::get_econv();
  cout << "Interaction of solute with protein charges = " << 
    protein_interaction << endl;
  cout << "(probably in kcal/mole)" << endl;

  float safe_eps_sol = PhysCond::get_epsext();
  PhysCond::set_epsext(epsvac);
  ElectrolyteEnvironment vac_ely;  // No electrolyte is the default
  DielectricEnvironment vac_eps(new TwoValueDielectricByAtoms(qpt, epsin1));
  ElstatPot vac_phi(fdm, vac_eps, rhogen, vac_ely);
  vac_phi.solve();
  float prod_vac = vac_phi * rhogen;
  cout << "prod_vac = " << prod_vac << endl;
  float reac_energy = (prod_sol - prod_vac) / 2 * PhysCond::get_econv();
  cout << "\n\nReaction field component of solvation = " << reac_energy
    << "\n(probably in kcal/mole)" << endl;
  float solvation_energy = reac_energy + protein_interaction;
  cout << "\n\nSOLVATION ENERGY IN PROTEIN = " << solvation_energy
    << "\n(probably in kcal/mole)" << endl;

  if (doReactionField) {

    ifstream fpt (fieldpoint_filename.c_str());
    if (!fpt.good())
      ::error("domeadsolinprot", "main: failed to open for reading, ",
	      fieldpoint_filename.c_str());
    ofstream rf(reacfield_filename.c_str());
    if (!rf.good())
      ::error("domeadsolinprot", "main: failed to open for writing, ",
	      reacfield_filename.c_str());

    *nqueries = 0;
    while (fpt.good()) {
      Coord c;
      fpt >> c;
      if (fpt.eof()) break;
      if (!fpt.good())
	::error("domeadsolinprot", "main: input failure while reading ",
		fieldpoint_filename.c_str());
      rf << phi.value(c) - vac_phi.value(c) << "\n";
      if (!rf.good())
	::error("domeadsolinprot", "main: output failure while writing ",
		reacfield_filename.c_str());
      ++(*nqueries);
    }
  }

  if (doProteinField) {

    blab1 << "Protein Field requested, so calculating field due to "
	  << "protein charges" << endl;
    PhysCond:: set_epsext(safe_eps_sol);
    ElstatPot prot_phi(fdm, eps, rhofeel, ely);
    prot_phi.solve();

    ifstream fpt (fieldpoint_filename.c_str());
    if (!fpt.good())
      ::error("domeadsolinprot", "main: failed to open for reading, ",
	      fieldpoint_filename.c_str());
    else
      cout << "Opened " << fieldpoint_filename.c_str() << " for reading."<< endl;

    ofstream pf(protfield_filename.c_str());
    if (!pf.good())
      ::error("domeadsolinprot", "main: failed to open for writing, ",
	      protfield_filename.c_str());
    else
      cout << "Opened " << protfield_filename.c_str() << " for writing."  <<endl;

    ofstream annfile("AnnsFile.out");
    if (!annfile.good())
      ::error("domeadsolinprot", "main: failed to open for writing, ",
	      "AnnsFile.out");
    else{
      cout << "Opened AnnsFile.out for writing."  <<endl;
      annfile<<" Here's a line" << endl;
      annfile.close();
    }

    while (fpt.good()) {
      Coord c;
      fpt >> c;
      if (fpt.eof()) break;
      if (!fpt.good())
	::error("domeadsolinprot", "main: input failure while reading ",
		fieldpoint_filename.c_str());
      cout << prot_phi.value(c) << endl;
      pf << prot_phi.value(c) << "\n";
      if (!pf.good())
	::error("domeadsolinprot", "main: output failure while writing ",
		protfield_filename.c_str());
    }
    
  }
    
}
