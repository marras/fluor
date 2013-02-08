#include "fluorescence.h"

/*
 * - Constructors for FluorRods class
 * - Diffusion of non-rotating rigid rods
 */

/* Simulates both diffusion of rods and fluorescence of the dye labels */
FluorRods :: FluorRods () {
	simulate_diffusion = true;
	current_frame = -1; fin = -1; natoms = 0; ROD_LENGTH = NULL;

	ReadGROMACSParameters ("grompp.mdp");

	if (natoms <= 0) LOG ("!Incorrect number of molecules (natoms). Check your grompp.mdp file.");
	if (num_rods == NULL) LOG ("!Insufficient data in grompp.mdp: rods number not set.");
//	if (ROD_LENGTH < 0) LOG ("!Insufficient/incorrect data in grompp.mdp: ROD LENGTH not set.");

	//Allocate mem for molecules & set random starting positions	
	int m = 0; // index czasteczki
	
	for (int t=0; t<types; t++) {	
	   for (int r=0; r<num_rods[t]; r++) {
		  double rod_initial_pos [3], orientation [3];
		  for (int d=0; d<3; d++) {
			  rod_initial_pos[d] = rng.rand (SIZE[d]); 	//NOTE this is the end of the rod, NOT the center of mass (doesn't matter so far)
			  orientation[d] = rng.rand (2.0) - 1.0;		 //orient: number in [-1.0,1.0]
		  }

		  Normalize (orientation);

		  //ustaw atoms_per_rod fluoroforow wzdluz calej dlugosci paleczki, jednorodnie
		  for (int a=0; a<atoms_per_rod[t]; a++) {
			  for (int d=0; d<3; d++) {
				  mol[m].x[d] = rod_initial_pos[d] + ((double) a)/((double) atoms_per_rod[t])*orientation[d]*ROD_LENGTH[t]; 
				  mol[m].init_x[d] = mol[m].x[d];
			  }
		   m++;
		  }
	   }
	}

	LOG ("Done setting up %d molecules.",m);

	//Init (); //Ta funkcja już została wywołana z konstruktora klasy bazowej Fluorescence () !! NOTE
}

/***************************************
 * Simulate diffusion for labeled RODS *
 ***************************************/

void FluorRods :: SimulateDiffusion () {
	static bool warningowanie = true;
		if (warningowanie) LOG ("*Calculating diffusion for rods.");
	warningowanie = false;

//	const double DIFFUSION_STEP = sqrt (2*diff_coeff * dT/steps_per_frame*1e-6); // dT = dt * steps_per_frame
//teraz jest zdefiniowany dla kazdej czasteczki (ale w obrebie 1 grupy beda mialy taki sam, wiec wezmy I czasteczke z danej paleczki)

	int m = 0; // index czasteczki

	for (int krok = 0; krok < steps_per_frame; krok++) {    //steps_per_frame = dT / dt (ile x gestsza rozdzielczosc sym. dyfuzji od sym. fluorescencji)
		for (int t = 0; t < types; t++) {
			for (int r=0; r < num_rods[t]; r++) {
				const double DIFFUSION_STEP = mol[m].DIFFUSION_STEP(); //aktualnie m - nr I barwnika z danej paleczki
				for (int d=0; d<3; d++) { 	//w kazdym wymiarze
					double displ = DIFFUSION_STEP * rng.randNorm(0.0,1.0);	//przesun wszystkie barwniki z tej samej paleczki o tyle samo (dlatego najpierw petla d, potem a)
					for (int a=0; a<atoms_per_rod[t]; a++) {
						mol[m+a].x[d] += displ;
						//Periodic boundary conditions
						mol[m+a].CheckSimulationBoxLimits(d);
					}
				}
				m += atoms_per_rod[t]; //przesuwamy index czasteczki do pierwszego barwnika z kolejnej paleczki
			}
		}
//	 LOG ("mol[0].x[0] = %lf",mol[0].x[0]);
	}

	//We can add more sophisticated interaction algorithms here.
	//TODO rotational diffusion! (read up!) - chbya najlepiej przeksztalcic jakimis macierzami obrotu?
}


/* Wczytywanie konfiguracji - wersja dla symulacji paleczek */

void FluorRods :: ReadGROMACSParameters (char * param_file) {
	ReadBasicGROMACSParameters (param_file);

	LOG ("Fluorescence: Reading further input parameters concerning %d types of rods...",types);
	
	std::ifstream fpar (param_file);
	char buf [200];
	
	num_rods = new int [types]; //number of molecules of in each group
	ROD_LENGTH = new double [types]; //number of molecules of in each group
	atoms_per_rod = new int [types];
	
	while (fpar.getline (buf, 200)) {
		std::string line (buf);
		
		int comment = line.find (";");
		line = line.substr(0,comment); //jesli nie znajdzie srednika, to comment=npos, czyli najwieksza liczba typu size_t - calosc sie kopiuje.

		if (line.length()==0) continue; //skip empty lines or lines with just comments.
		printf ("%d\n",line.length());
		  
		std::string zmienna (line.substr(0,line.find(" "))); //do pierwszej spacji mamy nazwe zmiennej wczytywanej z pliku
		zmienna  = zmienna.substr(0,zmienna.find("\t")); //(lub do pierwszego tabu)
		std::string wartosc = line.substr(line.find("=")+1,NPOS);
		while (wartosc[0]==' ' || wartosc[0] =='\t') wartosc = wartosc.substr(1,NPOS); //obetnij spacje i taby z lewej
		wartosc = wartosc.substr(0,wartosc.find(" ")); //obetnij spacje z prawej
		wartosc = wartosc.substr(0,wartosc.find("\t")); //obetnij taby z prawej
		printf ("%s = %s.\n", zmienna.c_str(),wartosc.c_str());
		fflush(stdout);

		if (zmienna == "rods" || zmienna == "ROD_LENGTH") {LOG ("!%s: this convention is now deprecated. Please specify the molecule type id (starting from 0), such as: NATOMS_0, NATOMS_1 etc.", zmienna.c_str());break;}

		else if (zmienna.find(std::string("rods_"))!=NPOS) { //to podejscie daje warning, ale chyba dziala HACK
		  int t = atoi (zmienna.substr(5).c_str());
		  num_rods[t] = atoi (wartosc.c_str());
		  LOG ("*Saving number of rods of [%d] type: %d",t,num_rods[t]);
		}
		else if (zmienna.find(std::string("ROD_LENGTH_"))!=NPOS) {
		  int t = atoi (zmienna.substr(11).c_str());
		  ROD_LENGTH[t] = atof (wartosc.c_str());
		  LOG ("*Saving length of rods of [%d] type: %lf",t,ROD_LENGTH[t]);
		}
	}

	fpar.close();

	for (int t=0; t<types; t++) {
	    if (num_rods[t] <= 0) LOG ("!Incorrect number of rods for component %d. Check your grompp.mdp file or make sure if you really want to simulate rod diffusion.",t);	  
	    if (num_atoms_in_type[t] % num_rods[t] != 0) LOG ("!Number of atoms (in group %d) should be an integer multiple of the number of rods!",t);
	    atoms_per_rod [t] = num_atoms_in_type[t] / num_rods[t];
	}
		
	
	LOG ("*FluorRods: dT = %lf, total_frames = %d, steps_per_frame = %d, TOTAL TIME = %lf, rod types = %d, total atoms in rods = %d",dT,total_frames, steps_per_frame, total_frames*dT, types, natoms);
//	#ifdef HARD_SPHERES
//	 LOG ("*Volume fraction (phi) = %0.3lf\%, nstlist = %d", natoms*4.0/3.0*3.14*MOL_SIZE*MOL_SIZE*MOL_SIZE/8.0/(8*X0*Y0*Z0)*100.0, nstlist);
//	#endif //HARD_SPHERES
}

FluorRods :: ~FluorRods () {
    delete [] ROD_LENGTH; delete [] num_rods; delete [] atoms_per_rod;
}

