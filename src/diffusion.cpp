/* FluorRods calculation module - rods */
#include "fluorescence.h"
#include <iostream>
#include <fstream>

#include <typeinfo>

// double CENTER [3] = {-1.0, -1.0, -1.0};		//the middle of the simulation box  // V = 4*4*10 = 160 um^3  
double SIZE[3] = {-1.0, -1.0, -1.0};		//the size of the simulation box  // V = 4*4*10 = 160 um^3  
double WXY = -1.0, WZ = -1.0;			//size of the excitation volume (Radius!)
double SQR_KAPPA = -1.0; //(WZ/WXY)*(WZ/WXY);	//SQR_KAPPA = KAPPA^2 = (wxy / wz)^2 - shape factor of the excitation volume (kappa = 3).
double OBS_VOL_XY = WXY*2;	//1000.0	//observation volume dimensions (1x1x3 um)
double OBS_VOL_Z = WZ*2;	//3000.0

using namespace std;
/******************************************
 * Simulate diffusion for FREE small dyes *
 ******************************************/

void Fluorescence :: SimulateDiffusion () {

#ifdef ENABLE_GPU
MoveMolecules_GPU (); return;
#endif

#define _NO_COL -1
//#define _RE_CHECK -2

if (ile_babli == 0) {
	//Normalne zachowane (box szescienny)
	for (int t=0; t<steps_per_frame; t++) {	    //steps_per_frame = dT / dt (ile x gestsza rozdzielczosc sym. dyfuzji od sym. fluorescencji)
		for (int i=0; i<natoms; i++) {
		  const double DIFFUSION_STEP = mol[i].DIFFUSION_STEP();
		  //if (DIFFUSION_STEP == 0) LOG ("!DIFFUSION_STEP==0 dla czastki %d!",i); //DEBUG
		  
		  for (int d=0; d<3; d++) {
			#ifdef HARD_SPHERES
			mol[i].switched[d] = 0;	   //poki co nie wyszlismy za box w zadna strone
			mol[i].coll_test_x[d]= mol[i].x[d];
			#endif //HARD_SPHERES

// 			mol[i].x[d] += DIFFUSION_STEP * (2.0*rng.rand() - 1.0); //uniform od -1 do +1 TEST HACK
			mol[i].x[d] += DIFFUSION_STEP * rng.randNorm(0.0,1.0); //rozklad normalny - taki jaki powinien byc...
			mol[i].CheckSimulationBoxLimits(d); //Periodic boundary conditions
 		 }

		#ifdef HARD_SPHERES
		  UpdateMoleculeMovement(i);
		#endif //HARD_SPHERES
		}
	
  	#ifdef HARD_SPHERES
		for (int i=0; i<natoms; i++) {
			int m = CheckCollision (i);
			if (m != _NO_COL) UndoCollision (i,m); //if there was any collision, undo it
		}
	#endif //HARD_SPHERES
	}
} else if (ile_babli == 1) {
	/************************************
	* CONSTRICTION WITHIN AN ELLIPSOID *
	************************************/
	for (int t=0; t<steps_per_frame; t++) {	    //steps_per_frame = dT / dt (ile x gestsza rozdzielczosc sym. dyfuzji od sym. fluorescencji)
		for (int i=0; i<natoms; i++) {
			double new_x[3];
			for (int d=0; d<3; d++) new_x[d] = mol[i].x[d] + mol[i].DIFFUSION_STEP() * rng.randNorm(0.0,1.0);

			if (babel->IsInside(new_x)) for (int d=0; d<3; d++) mol[i].x[d] = new_x[d];
			// w przeciwnym razie (gdyby czasteczka miala wyjsc poza babel) nie wykonujemy kroku
		}
	}
} else LOG ("!Wrong number of ile_babli.");

//	for (int m=0; m<natoms; m++) { //save positions of molecules TODO change to binary for later use
//		fprintf (fpos,"%d\t%f\t%f\t%f\n",m,mol[m].x[0], mol[m].x[1], mol[m].x[2]);
//	}

//monitor x dimension of 1st molecule
//  	LOG ("*Ruszylem sie o: %lf",mol[0].x[0] - mol[0].coll_test_x[0]);
// LOG ("*Mol[1]_x: %lf",mol[0].x[0]); //pozycja 1-szej czasteczki
}


void Fluorescence :: ReadBasicGROMACSParameters (char * param_file) {
	LOG ("Fluorescence: Reading Basic GROMACS input parameters...");

	ifstream fpar (param_file, ios::in);
	if (fpar.is_open()) {LOG( "*Opened file %s",param_file);}
	int _nsteps=-1, _dt=-1; steps_per_frame =-1;
	types = -1;
	num_atoms_in_type = NULL; //number of molecules of in each group
	double * diff_coeff = NULL;
	int * _dark = NULL;

	natoms = 0; dark = 0; //just to make sure.

	double _babel_size = -1, _babel_shape = -1, _babel_diff_coeff = -1;
	double bx=0, by=0, bz=0; //pozycja babla
	
	char buf [200];
	
/*	cout << typeid(atof("1000")).name() << " " <<typeid(diff_coeff).name() << " " <<typeid(1).name() << " " <<typeid(10.0).name() <<endl;*/
	
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

		if (zmienna == "nsteps") { _nsteps = atoi (wartosc.c_str()); LOG ("nsteps= %d",_nsteps); }
		else if (zmienna == "types") {  
		    types = atoi (wartosc.c_str()); 
		    LOG ("types= %d",types);
		    num_atoms_in_type = new int [types];
		    diff_coeff = new double [types];
		    _dark = new int [types];
		}
		else if (zmienna == "dt") {  _dt = atoi (wartosc.c_str()); LOG ("dt= %d",_dt); }
		else if (zmienna == "nstxtcout" || zmienna == "steps_per_frame") { steps_per_frame = atoi (wartosc.c_str()); LOG ("xtcout (steps per frame) = %d",steps_per_frame);}
		else if (zmienna == "NATOMS" || zmienna=="diff_coeff" || zmienna=="dark") {LOG ("!%s: this convention is now deprecated. Please specify the molecule type id (starting from 0), such as: NATOMS_0, NATOMS_1 etc.", zmienna.c_str());break;} 

		else if (zmienna.find(std::string("NATOMS_"))!=NPOS && simulate_diffusion) { //to podejscie daje warning, ale chyba dziala HACK
		  int t = atoi (zmienna.substr(7).c_str());
		  num_atoms_in_type[t] = atoi (wartosc.c_str());
		  natoms += num_atoms_in_type[t]; //natoms - now TOTAL number of atoms.
		  LOG ("*SAVING PARAMS for atom NATOMS of type: %d = %d",t,num_atoms_in_type[t]);
		}
		else if (zmienna.find(std::string("diff_coeff_"))!=NPOS && simulate_diffusion) {
		  int t = atoi (zmienna.substr(11).c_str());
		  diff_coeff[t] = atof (wartosc.c_str());
		  LOG ("*diff coeff [%d] = %lf",t,diff_coeff[t]);
		}
		else if (zmienna.find(std::string("dark_"))!=NPOS && simulate_diffusion) {
		  int t = atoi (zmienna.substr(5).c_str());
		  _dark[t] = atoi (wartosc.c_str());
		  dark += _dark[t];
		  LOG ("*_dark [%d] = %d",t,_dark[t]);
		}
		  
		else if (zmienna == "epsilon") { EPSILON = atof(wartosc.c_str()); LOG ("EPSILON= %lf",EPSILON);}
		else if (zmienna == "Qf") { Qf = atof(wartosc.c_str()); LOG ("Qf= %lf",Qf);} 
 		else if (zmienna == "Qt") { Qt = atof(wartosc.c_str()); LOG ("Qt= %lf",Qt);} 
		else if (zmienna == "Qb") { Qb = atof(wartosc.c_str()); LOG ("Qb= %lf",Qb);}
		else if (zmienna == "tauT") { tauT = atof(wartosc.c_str()); LOG ("tauT= %lf",tauT);}
		else if (zmienna == "pinhole") { pinhole = atoi(wartosc.c_str()); LOG("pinhole= %d",pinhole);}
		else if (zmienna == "nstlist") { nstlist = atoi(wartosc.c_str()); LOG("nstlist= %d",nstlist);}
		else if (zmienna == "SIZE_X") {SIZE[0] = atof(wartosc.c_str()); LOG("SIZE_X= %lf",SIZE[0]);}
		else if (zmienna == "SIZE_Y") {SIZE[1] = atof(wartosc.c_str()); LOG("SIZE_Y= %lf",SIZE[1]);}		
		else if (zmienna == "SIZE_Z") {SIZE[2] = atof(wartosc.c_str()); LOG("SIZE_z= %lf",SIZE[2]);}
		else if (zmienna == "WXY") {WXY = atof(wartosc.c_str()); LOG("WXY= %lf",WXY);}
		else if (zmienna == "WZ") {WZ = atof(wartosc.c_str()); LOG("WZ= %lf",WZ);}
		else if (zmienna == "OBS_VOL_XY") {OBS_VOL_XY = atof(wartosc.c_str()); LOG("OBS_VOL_XY= %lf",OBS_VOL_XY);}
		else if (zmienna == "OBS_VOL_Z") {OBS_VOL_Z = atof(wartosc.c_str()); LOG("OBS_VOL_Z= %lf",OBS_VOL_Z);}
		else if (zmienna == "profile") {
			    if (wartosc == "gauss") profil = PR_GAUSS; 
			    else if (wartosc == "square") profil = PR_SQUARE; 
			    LOG("profil= %s=%d",wartosc.c_str(),profil);
		}
		else if (zmienna == "Focus") {num_foci = atoi(wartosc.c_str()); LOG("num_foci= %d",num_foci);}
		else if (zmienna == "f1x") {F1[0] = atof(wartosc.c_str()); LOG("F1[0]= %lf",F1[0]);}
		else if (zmienna == "f1y") {F1[1] = atof(wartosc.c_str()); LOG("F1[1]= %lf",F1[1]);}
		else if (zmienna == "f1z") {F1[2] = atof(wartosc.c_str()); LOG("F1[2]= %lf",F1[2]);}
		else if (zmienna == "f2x") {F2[0] = atof(wartosc.c_str()); LOG("F2[0]= %lf",F2[0]);}
		else if (zmienna == "f2y") {F2[1] = atof(wartosc.c_str()); LOG("F2[1]= %lf",F2[1]);}
		else if (zmienna == "f2z") {F2[2] = atof(wartosc.c_str()); LOG("F2[2]= %lf",F2[2]);}

		else if (zmienna == "babel") {ile_babli = atoi(wartosc.c_str()); LOG("Liczba babli = %d",ile_babli);}
		else if (zmienna == "babel_size") {_babel_size = atof(wartosc.c_str()); LOG("_babel_size = %lf",_babel_size);}
		else if (zmienna == "babel_shape") {_babel_shape = atof(wartosc.c_str()); LOG("_babel_shape = %lf",_babel_shape);}
		else if (zmienna == "babel_diff_coeff") {_babel_diff_coeff = atof(wartosc.c_str()); LOG("_babel_diff_coeff = %lf",_babel_diff_coeff);}
		else if (zmienna == "babel_x") {bx = atof(wartosc.c_str()); LOG("babel_x = %lf",bx);}
		else if (zmienna == "babel_y") {by = atof(wartosc.c_str()); LOG("babel_y = %lf",by);}
		else if (zmienna == "babel_z") {bz = atof(wartosc.c_str()); LOG("babel_z = %lf",bz);}
	#ifdef ENABLE_GPU
		else if (zmienna == "NThreads") {NThreads = atoi(wartosc.c_str()); LOG("GPU threads per block = %d",NThreads);}
		else if (zmienna == "NBlocks") {NBlocks = atoi(wartosc.c_str()); LOG("GPU blocks = %d",NBlocks);}
	#endif	
// 		else LOG ("!Surprised by a weird fault in the program! This should never happen!"); //nie plujmy sie, niektore zmienne sa przydatne tylko przy konkretnych rodzajach symulacji
	}

	fpar.close();

	if (natoms <= 0) LOG ("!Incorrect number of molecules (natoms). Check your grompp.mdp file.");
	if (SIZE[0] < 0 || SIZE[1] < 0 || SIZE[2] < 0) LOG ("!Missing data on simulation box size in grompp.mdp.");
	if (_nsteps == -1 || _dt == -1 || steps_per_frame == -1 || natoms == -1 || diff_coeff < 0) { LOG ("!Insufficient data in grompp.mdp"); exit(71);}
	if (dark < 0 || dark > natoms) LOG ("!Incorrect number of dark molecules.");
	if (_nsteps % steps_per_frame != 0) LOG ("!nsteps should be an integer multiple of steps_per_frame!");
	
	if (num_foci == 2) LOG ("*Starting dual-focus experiment.");
	else if (num_foci != 1) LOG ("!Incorrect foci number: %d.",num_foci);

	SQR_KAPPA = (WZ/WXY)*(WZ/WXY); //parametr struktury ogniska
	
	LOG("Wczytalismy plik!");
	
	dT = _dt * steps_per_frame;
	total_frames = _nsteps / steps_per_frame;

	if (ile_babli > 0) {
		LOG ("*Constricting diffusion to an ellipsoid bubble.");
		if (ile_babli > 1) LOG ("!More than one bubble is not implemented yet!");
		if (_babel_size < 0 || _babel_shape < 0 || _babel_diff_coeff != 0) LOG ("!Incorrect bubble size/shape/speed data!");
		babel = new Bubble (bx,by,bz,_babel_size,_babel_shape,_babel_diff_coeff, dT, steps_per_frame); //create new bubble TODO change to array
	}

	
	//----------------- KONIEC WCZYTYWANIA, POCZATEK INICJALIZACJI -----------------//
	
	int m = 0; //ktora czasteczka (liczac ogolnie)
	mol = new Molecule [natoms]; //Allocate mem for molecules & set random starting positions
	state = new enumStates [natoms];
	ile_wyswiecilam = new int [natoms];	/////

	for (int t=0;t<types;t++) {
	    for (int i=0;i<num_atoms_in_type[t];i++,m++) {
		ile_wyswiecilam[m] = 0; 
		mol[m].SetDiffusionCoeff (diff_coeff[t], dT, steps_per_frame); //set diff. speed for each molecule
		
		if (i<_dark[t]) state[m]=MS_DARK;  //allow for non-fluorescent molecules
		else state[m]=MS_GROUND; 	   //all the other molecules start in ground state	
	    }
	}	

	delete [] diff_coeff;
	//delete [] num_atoms_in_type;  //to sie jeszcze przyda (dla diff_rods, jesli sa)
}

/* -------------------- Wczytywanie parametrow z pliku -------------------*/
void Fluorescence :: ReadGROMACSParameters (char * param_file) {

	ReadBasicGROMACSParameters (param_file);
	LOG ("*Fluorescence: dT = %lf, total_frames = %d, steps_per_frame = %d, TOTAL TIME = %lf, NATOMS = %d",dT,total_frames, steps_per_frame, total_frames*dT, natoms);

	#ifdef HARD_SPHERES
	LOG ("*Volume fraction (phi) = %0.3lf\%, nstlist = %d", natoms*4.0/3.0*3.14*MOL_SIZE*MOL_SIZE*MOL_SIZE/8.0/(8*X0*Y0*Z0)*100.0, nstlist);
	#endif //HARD_SPHERES
}

/*******************
 * Calculate <x^2> *
 *******************/

double Fluorescence :: MeanSquareDisplacement () {
	double displacement = 0;
	for (int i=0; i<natoms; i++) 
		displacement += mol[i].MSD();

	return displacement/natoms;
}

/* Auxilliary function: square distance */
inline double SqrDistance (rvec i, rvec j) {		//square distance - zeby nie obliczac sqrt()
    double _x = i[0]-j[0], _y = i[1]-j[1], _z = i[2]-j[2];
    return _x*_x + _y*_y + _z*_z;
}


/***********************************************************
 * HARD SPHERE COLLISION ALGORITHMS (TODO not 100% secure) *
 ***********************************************************/
#ifdef HARD_SPHERES
void Fluorescence :: UndoCollision (int i, int m) {
	if (i<0 || m<0) LOG ("!Trying to undo an inexistent collision: %d & %d",i,m);

	for (int d=0;d<3;d++) {
//		mol[i].x[d] = mol[i].coll_test_x[d] + DIFFUSION_STEP * rng.randNorm(0.0,1.0);	//wykonujemy zamiast tego inny ruch obiema czasteczkami NOTE te o malych numerach sa faworyzowane do kolizji!
//		mol[m].x[d]=  mol[m].coll_test_x[d] + DIFFUSION_STEP * rng.randNorm(0.0,1.0);
		mol[i].x[d] = mol[i].coll_test_x[d];	//tylko cofamy czasteczki! Metropolis.
		mol[m].x[d] = mol[m].coll_test_x[d];

		//poprawka dla init_x przy kolizji po przejsciu przez sciane
		mol[i].init_x[d] += mol[i].switched[d] * 2*CENTER[d]; //switched = -1, 0, lub 1 w zal. w ktora strone
		mol[m].init_x[d] += mol[m].switched[d] * 2*CENTER[d]; //switched = -1, 0, lub 1 w zal. w ktora strone

		mol[i].switched[d] = 0; mol[m].switched[d] = 0;
	}
	UpdateMoleculeMovement (i); UpdateMoleculeMovement (m); //update sector data

	int n = CheckCollision (i);
	if (n != _NO_COL) UndoCollision (i,n); //HACK dalej juz rekurencyjnie, jezeli cos zajelo ich miejsce
	n = CheckCollision (m);
	if (n != _NO_COL) UndoCollision (m,n);

// static const double DIFFUSION_STEP = sqrt (2*1.381e-17*temperature*dT/steps_per_frame/(friction*1e-12)); //dt = dT / steps_per_frame NOTE static prevents multiple calculations

// 			for (int d=0;d<3;d++) {
// 				
// 				//Periodic boundary conditions
// 				mol[i].CheckSimulationBoxLimits(d);
// 				mol[m].CheckSimulationBoxLimits(d);
// 			}
}

inline int Fluorescence :: CheckCollision (const int & i) {
    int col = _NO_COL; //which molecule does #i collide with?

    int D [3]; //zmienna do iterowania po sektorach w 3 wymiarach (3D ;) )

    for (D[0]=-1; D[0]<=1; D[0]++) 
      for (D[1]=-1; D[1]<=1; D[1]++) 	//dla wszystkich 3 wym. bierzemy sekt. sasiadujace: poprzedni,obecny,nastepny
	 for (D[2]=-1; D[2] <=1; D[2]++) {
	    if (col != _NO_COL) break;

	    int _x[3];

	    for (int d=0; d<3; d++) { 	//enable through-the-wall interactions NEW
		//wylicz, z jakim sektorem sie stykamy przez sciane:
		//trzeba dodac NUM_SECTORS[0], bo tutaj modulo (%) dziala dla ujemnych inaczej niz
		//np. w Pythonie: (-1) % 5 = -1, a nie 4 !
		_x[d] = (mol[i].sec[d]+D[d]+NUM_SECTORS[d]) % NUM_SECTORS[d];

		//do sprawdzania kolizji przez sciane przesun czasteczki na druga strone pudla
		if (mol[i].sec[d] == 0 && D[d] == -1) mol[i].x[d] += 2*CENTER[d];
		if (mol[i].sec[d] == NUM_SECTORS[d]-1 && D[d] == 1) mol[i].x[d] -= 2*CENTER[d];
	    }

	    for (int j=0; j<num_sec[_x[0]][_x[1]][_x[2]]; j++) { //dla wszystkich czast. z sektora x1y1z1
		    int m = sector[_x[0]][_x[1]][_x[2]][j];
		    //do not compare a molecule with itself (or ones that have already been checked??)
		    if (i == m) continue;
//  			LOG(" #%d : Checking molecules %d & %d.\n ",current_frame,i,*j);
		    if (SqrDistance (mol[i].x,mol[m].x) < MOL_SIZE*MOL_SIZE) { //jezeli jest kolizja
// 			LOG (" #%d : Collision: %d & %d.\n ",current_frame,i,m);
			col = m;
			break;
		    }
	    }

	    for (int d=0; d<3; d++) { //reverse the changes KONIECZNIE!
		if (mol[i].sec[d] == 0 && D[d] == -1) mol[i].x[d] -= 2*CENTER[d];
		if (mol[i].sec[d] == NUM_SECTORS[d]-1 && D[d] == 1) mol[i].x[d] += 2*CENTER[d];
	    }
/*	    if (mol[i].sec[0] == 0 && d0 == -1) mol[i].x[0] -= 2*CENTER[0];  zamiast tego :P
	    if (mol[i].sec[1] == 0 && d1 == -1) mol[i].x[1] -= 2*CENTER[1];
	    if (mol[i].sec[2] == 0 && d2 == -1) mol[i].x[2] -= 2*CENTER[2];
	    if (mol[i].sec[0] == NUM_SECTORS[0]-1 && d0 == 1) mol[i].x[0] += 2*CENTER[0];
	    if (mol[i].sec[1] == NUM_SECTORS[1]-1 && d1 == 1) mol[i].x[1] += 2*CENTER[1];
	    if (mol[i].sec[2] == NUM_SECTORS[2]-1 && d2 == 1) mol[i].x[2] += 2*CENTER[2];*/
    }
    return col;	//nie bylo kolizji
}

void Fluorescence :: UpdateMoleculeMovement (const int & i) {
	const int s[3] = {mol[i].sec[0], mol[i].sec[1], mol[i].sec[2]};
	mol[i].UpdateSec();	//czasteczka sprawdza, w ktorym sektorze jest
	for (int d=0;d<3;d++) { //jesli chociaz w 1 wymiarze zmienil sie sektor
		if (mol[i].sec[d] != s[d]) {  //to trzeba uaktualnic pozycje czasteczki w liscie sektorow
			remove_from_sector(s[0],s[1],s[2],i);			//wersja dla array'a
			push_to_sector(mol[i].sec[0], mol[i].sec[1], mol[i].sec[2], i);
			break;
		}
	}
}
#endif //HARD_SPHERES


#ifdef HARD_SPHERES
/*
void Fluorescence :: UpdateNeighbourList () {
//	LOG ("*Updating neighbour list...");

//	int nei = 0;
     for (int i=0; i<natoms; i++) {
	int brzeg [3]; 
//	brzeg[0]=brzeg[1]=brzeg[2] = 0;	//0 - nie przy brzegu, 1 - z lewej, -1 - z prawej

	neighbours[i].clear();

	//compare the closest ones
	for (int j=0; j<i; j++) {
		if (SqrDistance (x[i],x[j]) < NSTLIST_CUTOFF*NSTLIST_CUTOFF) { //compare sqr distances - faster
			neighbours[i].push_back(j);
			neighbours[j].push_back(i);
		}
	}
	
	//Take care of molecules near simulation box boundaries
	for (int d=0;d<3;d++) {
		if (mol[i].x[d]<NSTLIST_CUTOFF)	//jestesmy blisko lewego brzegu
			brzeg[d] = 1;
		else if (2*CENTER[d] - mol[i].x[d] < NSTLIST_CUTOFF)
			brzeg[d] = -1;
		else brzeg[d] = 0;
	}	
			
	int ile_przy_brzegu
	
	
// 	for (int j=0; j<i; j++) {
// 		if (brzeg[0]) {
// 			x[i]d
// 			if (brzeg[1]) {
// 				if (brzeg[2]) {
// 			
// 			SqrDistance (x[i],x[j]) < NSTLIST_CUTOFF*NSTLIST_CUTOFF) { //compare sqr distances - faster
//			neighbours[i].push_back(j);
// 			neighbours[j].push_back(i);
//		}
//	}
	
	
	
	//////////////
			for (int j=0; j<i; j++) {
				if (2*CENTER[d] - mol[j].x[d] + mol[i].x[d] < NSTLIST_CUTOFF) { //NOTE niepotrzebnie tak duzo (2x cutoff range)
					neighbours[i].push_back(j);
					neighbours[j].push_back(i);
				}
			}
		}
		else if (2*CENTER[d]-mol[i].x[d] < NSTLIST_CUTOFF) { //jestesmy blisko prawego brzegu
// 			LOG ("*%d is much to the right (dim %d)!",i,d);
			for (int j=0; j<i; j++) {
				if (mol[j].x[d] < NSTLIST_CUTOFF) {
					neighbours[i].push_back(j);
					neighbours[j].push_back(i);
				}
			}
		}
	}

//	nei += neighbours[i].size();
     }
//	LOG ("*AVG nei = %5.2f",((float)nei) / natoms);
}
*/

// for -1,0,1: check sector i % NUM_SECTORS!


// bool Fluorescence :: CheckAndAvoidCollisions() {
//     //HARD SPHERES
// /*    ktora_proba++;
//     if (ktora_proba >2) LOG ("*[#%d] Proba: %d",current_frame, ktora_proba);*/
//     bool bylo_zderzenie = false;
// 
//     for (int i=0; i<natoms; i++) { 
// 	if (CheckCollision (i)) bylo_zderzenie = true;
//     }
//     
//     return bylo_zderzenie;
// }
#endif //HARD_SPHERES

