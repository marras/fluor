#include "fluorescence.h"

//TODO dobrac odpowiedni rozmiar i gestosc przeszkod
// - na tyle drobnoziarniste, zeby zmiescilo sie n (~20-50) szescianikow w obj. konf
// - na tyle gruboziarniste, zeby nie przeskoczyc w jednym kroku!

FluorObstacles :: FluorObstacles () {
	simulate_diffusion = true; pinhole=0;
	current_frame = -1; fin = -1; natoms = -1;
	dark =-1;

	ReadGROMACSParameters ("config.dat");
	if (natoms <= 0) LOG ("!Incorrect number of molecules (natoms). Check your config.dat file.");

	//Allocate mem for molecules & set random starting positions
	mol = new Molecule [natoms];
	
	for (int i=0; i<natoms; i++) {
		for (int d=0; d<3; d++) {
			mol[i].x[d] = rng.rand (SIZE[d]);
			mol[i].init_x[d] = mol[i].x[d];
		}
			
//NOTE zmiana - elipsoida!!!
/* 	if ((mol[i].init_x[0]-CENTER[0])*(mol[i].init_x[0]-CENTER[0])/CENTER[0]/CENTER[0] + 
		    (mol[i].init_x[1]-CENTER[1])*(mol[i].init_x[1]-CENTER[1])/CENTER[1]/CENTER[1] + 
		    (mol[i].init_x[2]-CENTER[2])*(mol[i].init_x[2]-CENTER[2])/CENTER[2]/CENTER[2] > 1) {i--; continue;} //retry if molecule is outside the volume ///////////////// */
	}

	Init ();
	LOG ("Profile: %d",profil);
}


void FluorObstacles :: SimulateDiffusion () {
LOG("!Diffusion between obstacles is disabled in this version. Check src/diffusion_obstacles.cpp.");
  /*
   //DIFFUSION STEP ma zalezec od polozenia!
   const double DIFFUSION_STEP = sqrt (2*diff_coeff * dT/steps_per_frame*1e12); //dt = dT / steps_per_frame

	//Normalne zachowane (box szescienny)
	for (int t=0; t<steps_per_frame; t++) {	    //steps_per_frame = dT / dt (ile x gestsza rozdzielczosc sym. dyfuzji od sym. fluorescencji)
		for (int i=0; i<natoms; i++) {
//		    double a [] = {mol[i].x[0]/period, mol[i].x[1]/period, mol[i].x[2]/period};

	      //// Wersja dla podziału pół na pół: szary-biały ///
		double STEP = mol[i].x[2] > SIZE[2]/2 ? DIFFUSION_STEP / 10.0 : DIFFUSION_STEP;

		if (STEP == DIFFUSION_STEP) printf ("1\n");
		else printf ("0\n"); 	//FIXME! FAK FAK FAK! one utykaja na obszarach szarych!
		
		for (int d=0; d<3; d++) {
			mol[i].x[d] += STEP * rng.randNorm(0.0,1.0); 
// 			if (a[0]-((int)a[0]) < relative_size && a[1]-((int)a[1]) < relative_size && a[2]-((int)a[2]) < relative_size)
// 				mol[i].x[d] += (DIFFUSION_STEP * rng.randNorm(0.0,1.0) / 10.0) ; //jesli czasteczka w szarym szescianie, to chodzi 100x wolniej
// 			else 	mol[i].x[d] += DIFFUSION_STEP * rng.randNorm(0.0,1.0); //NOTE rozklad normalny - wolniejszy niz jednorodny, ale tak ma byc z row. Langevina...

  

			
			mol[i].CheckSimulationBoxLimits(d); //Periodic boundary conditions
 		 }
		}
	
	}

// 	for (int m=0; m<natoms; m++) { //save positions of molecules TODO change to binary for later use
// 		fprintf (fpos,"%d\t%f\t%f\t%f\n",m,mol[m].x[0], mol[m].x[1], mol[m].x[2]);
// 	}

//monitor x dimension of 1st molecule
//	LOG ("*Ruszylem sie o: %lf",mol[0].x[0] - mol[0].coll_test_x[0]);
*/
}


void FluorObstacles :: ReadGROMACSParameters (char * param_file) {
	ReadBasicGROMACSParameters (param_file);
	LOG ("Fluorescence: Reading further input parameters concerning Obstacles...");

	std::ifstream fpar (param_file);
	obst_size = 0; dim_obstacles = -1;

	char buf [200];
	
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

		if (zmienna == "obstacle_size") {obst_size = atof (wartosc.c_str()); LOG ("obst_size= %d",obst_size);}
		else if (zmienna == "dim_obstacles") {dim_obstacles = atoi (wartosc.c_str()); LOG ("dim_obstacles= %d",dim_obstacles);}  
	}

	if (dim_obstacles == -1) LOG ("!Insufficient data in config.dat: obstacles number not set.");
	if (obst_size == 0.0) LOG ("!Insufficient data in config.dat: obstacle size not set.");

	period = SIZE[0]/dim_obstacles; //odl. miedzy srodkami szecianow; zakladamy, ze sa rownomiernie rozm. w 3D
	relative_size = obst_size / period;
	LOG( "*period = %lf, rel_size=%lf",period, relative_size);
	
	LOG ("*FluorObstacles: dT = %lf, total_frames = %d, steps_per_frame = %d, TOTAL TIME = %lf, obstacles: size %f, num %d",dT,total_frames, steps_per_frame, total_frames*dT, obst_size, dim_obstacles);

	fpar.close();
}
