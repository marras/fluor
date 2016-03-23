/* FluorRods calculation module - rods */
#include "fluorescence.h"
#include <iostream>
#include <fstream>
#include <typeinfo>

// double CENTER [3] = {-1.0, -1.0, -1.0};	//the middle of the simulation box  // V = 4*4*10 = 160 um^3  
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

if (ile_babli == 0) {
	//Normalne zachowane (box szescienny)
	for (int t=0; t<steps_per_frame; t++) {	    //steps_per_frame = dT / dt (ile x gestsza rozdzielczosc sym. dyfuzji od sym. fluorescencji)
		for (int i=0; i<natoms; i++) {
		  const double DIFFUSION_STEP = mol[i].DIFFUSION_STEP();
		  //if (DIFFUSION_STEP == 0) LOG ("!DIFFUSION_STEP==0 dla czastki %d!",i); //DEBUG
		  
		  for (int d=0; d<3; d++) {
// 			mol[i].x[d] += DIFFUSION_STEP * (2.0*rng.rand() - 1.0); //uniform od -1 do +1 TEST HACK
			mol[i].x[d] += DIFFUSION_STEP * rng.randNorm(0.0,1.0); //rozklad normalny - taki jaki powinien byc...
			mol[i].CheckSimulationBoxLimits(d); //Periodic boundary conditions
 		 }
		}

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

	if (natoms <= 0) LOG ("!Incorrect number of molecules (natoms). Check your config.dat file.");
	if (SIZE[0] < 0 || SIZE[1] < 0 || SIZE[2] < 0) LOG ("!Missing data on simulation box size in config.dat.");
	if (_nsteps == -1 || _dt == -1 || steps_per_frame == -1 || natoms == -1 || diff_coeff < 0) { LOG ("!Insufficient data in config.dat"); exit(71);}
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


