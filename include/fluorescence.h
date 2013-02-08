#ifndef _FCS_FLUORESCENCE_H_
#define _FCS_FLUORESCENCE_H_

#include <fstream>
#include "MersenneTwister.h"
#include <stdio.h>
#include <vector>
#include "log.h"
#include <string>

#include <assert.h>

#ifndef NPOS
  #define NPOS -1
#endif

// #define ENABLE_GPU //turn on CUDA! ### THIS is now done via the MAKEFILE!!! NOTE DO NOT TURN ON MANUALLY!

// #ifdef ENABLE_GPU
//  #error Miałem tego nie wlaczac!
// #endif

extern MTRand rng;

typedef double rvec[3];

//NOTE default units as in GROMACS: picoseconds, nanometers

// PARAMETERS FOR FLUORESCENCE CALCULATION
// NOTE x[i][j] - i-ta czasteczka, j-ty skladnik wektora(x/y/z)
// NOTE all times in [ps], all distances in [nm]

extern double SIZE [3];			//size of simulation box
extern double WXY, WZ;			//size of the excitation volume (Radius!)
extern double SQR_KAPPA; //(WZ/WXY)*(WZ/WXY);	//SQR_KAPPA = KAPPA^2 = (wxy / wz)^2 - shape factor of the excitation volume (kappa = 3).

extern double EPSILON;  //1e6 / 1e12; 	//Absorption coefficient = specific brightness / 10^12 ( units: counts/(molecule*ps) )(?)
extern double Qf, Qt, Qb;		//Quantum yields of Fluorescence, ISC(Triplet), and Photobleaching
extern double tauT;			//exponential time constant for Triplet->Ground conversion, 5 us = 5*10^6 ps	

extern double OBS_VOL_XY;	//1000.0	//observation volume dimensions (1x1x3 um)
extern double OBS_VOL_Z;	//3000.0


//#define CALC_MSD	//for diffusion coefficient calculations
#define PROGRESS_IND_FREQ 5000
// #define HARD_SPHERES
const double MOL_SIZE = 0; //10; // 2*233; 	//size (diameter) of a dye molecule NOTE only matters for HARD_SPHERES collision checking!

// For Neighboursearching algorithm
//const double NST_MAGIC_FACTOR = 3.5;
//NOTE this should actually depend also on the avg brownian movement distance and nstlist
//const double NSTLIST_CUTOFF = NST_MAGIC_FACTOR * MOL_SIZE; //3.0 - magic constant (2.0 + rozmiar samych czasteczek)

#ifdef HARD_SPHERES
// For Sector algorithm
const double SECTOR_SIZE = 500.0; //rozmiar sektora do sprawdzania kolizji
const int MAX_LIST_LEN = 9; //int (SECTOR_SIZE*SECTOR_SIZE*SECTOR_SIZE / 3.141*6.0/(MOL_SIZE*MOL_SIZE*MOL_SIZE) ) + 1; //2.0, b0 4/3piR^3, ale MOL_SIZE = 2R
//const int NUM_SECTORS [3] = {int (X0*2/SECTOR_SIZE), int (Y0*2/SECTOR_SIZE), int (Z0*2/SECTOR_SIZE)};
const int NUM_SEC_X = 4;
const int NUM_SEC_Y = 4;
const int NUM_SEC_Z = 12;
const int NUM_SECTORS [3] = {NUM_SEC_X , NUM_SEC_Y , NUM_SEC_Z };
#endif //HARD_SPHERES

enum enumProfiles {
	PR_GAUSS = 0,
	PR_SQUARE,
};

enum enumErrorNumbers {
	ERR_BAD_DATA = 1,
	ERR_EOF,
};

enum enumStates {	// Molecule Excitation States
  MS_GROUND = 0,	//ground state
  MS_EXC_1,		//excited by laser 1
  MS_EXC_2,		//excited by laser 2
  MS_TRI,		//triplet
  MS_PHB,		//photobleached
  MS_DARK,		//not fluorescent at all (HS "filler")
};

/***********************************
 * Fluorescence computation module *
 ***********************************/
class Molecule {	//teraz bedzie dluzej: "mol[m].x[d]" zamiast  "x[m][d]", ale przejrzysciej!
   private:
     	double _diffusion_step;	//calculated from diff_coeff separately for each molecule (read only for external functions!)

   public:
	rvec x;			//molecule position vectors, x[j] -> j-x/y/z
	rvec old_x;		//for monitoring if a photobleached molecule has crossed the periodic box boundaries
	rvec init_x;		//for calculating <x^2>

	double DIFFUSION_STEP () {return _diffusion_step;} //getter method
	void SetDiffusionCoeff (double diff_coeff, double dT, int steps_per_frame) {
	    this->_diffusion_step = sqrt (2*diff_coeff * dT/steps_per_frame*1e-6); //dt = dT / steps_per_frame
// 	    LOG("*Hello! I'm setting diff coeff from %lf to DIFF_STEP = %lf",diff_coeff,DIFFUSION_STEP());
	}

	Molecule () {_diffusion_step = 0;} //W razie czego zauwazymy, ze nie jest zainicjalizowane, bo nie beda sie ruszac :)


      #ifdef HARD_SPHERES
	rvec coll_test_x;	//for storing the molecule positions to revert to after a collision
	int sec [3];		//w ktorym sektorze jestesmy
	inline void UpdateSec() {	//funkcja sprawdza, w ktorym sektorze jestesmy i uaktualnia zmienna sec
		for (int d=0;d<3;d++) {
			sec[d] = int (x[d]/SECTOR_SIZE);
			//NOTE Uratuj sie przed niepewnosciami dodawania double'i np. -1e-12 + 2000.0 = 2000.0 czyli sektor 20 zamiast 19!
			if (sec[d]==NUM_SECTORS[d]) {/*LOG ("*INVALID SECTOR: %lf %lf %lf!",x[0],x[1],x[2]); */sec[d]=NUM_SECTORS[d]-1;}  
		}
	}

	int switched[3];
      #endif //HARD_SPHERES

	inline int CheckSimulationBoxLimits (int d) { //d- ktory wmiar
		if (x[d] > SIZE[d]) {
			x[d] -= SIZE[d];
		#ifdef CALC_MSD
			init_x[d] -= SIZE[d];
		#endif
	  #ifdef HARD_SPHERES
		switched[d] = 1;	//wyszedl w prawo
	  #endif //HARD_SPHERES		
		}
		else if (x[d] < 0) {
			x[d] += SIZE[d];
		#ifdef CALC_MSD
			init_x[d] += SIZE[d];
		#endif
	  #ifdef HARD_SPHERES
		switched[d] = -1;	//wyszedl w lewo
	  #endif //HARD_SPHERES		
		}
	  #ifdef HARD_SPHERES
		return switched[d];
	  #else 
		return 0;
	  #endif //HARD_SPHERES
	}
	
	double MSD () {
		#ifndef CALC_MSD
			LOG("!MSD() called while CALC_MSD was not defined!");
		#endif
	   	return (x[0]-init_x[0])*(x[0]-init_x[0]) + (x[1]-init_x[1])*(x[1]-init_x[1]) 
			+ (x[2]-init_x[2])*(x[2]-init_x[2]);
	}
  // #endif //HARD_SPHERES
};

class Bubble {
   private:
	double pos [3];
	double sqr_size;
	double sqr_shape;
	double _diffusion_step;
	   
   public:
	Bubble(double x,double y,double z,double _size,double _shape,double diff_coeff, double dT, int steps_per_frame) {
	    this->_diffusion_step = sqrt (2*diff_coeff * dT/steps_per_frame*1e-6);
	    pos[0]=x;pos[1]=y;pos[2]=z;
	    this->sqr_shape = _shape*_shape; this->sqr_size = _size*_size;
	}
	
	~Bubble() {}
	
	bool IsInside (const double x [3]) {
	    return ((x[0]-pos[0])*(x[0]-pos[0])+ 
			(x[1]-pos[1])*(x[1]-pos[1])+ 
			(x[2]-pos[2])*(x[2]-pos[2])/sqr_shape <= sqr_size);	//Elipsoida
	}
	
	void Move();
	
	
};

class Fluorescence {
  private:
	bool bOK;
	
  protected:
    #ifdef ENABLE_GPU
	int NBlocks, NThreads; //dzielimy macierz docelowa poziomo na NBlokow, a pionowo na NWatkow :) -do GPU
	unsigned * seed;  
    #endif //ENABLE_GPU

	Bubble * babel;
	int ile_babli;

        double F1[3], F2[3];	//focus centers
	int num_foci;		//ile ognisk (1 lub 2)

	int types; //ile rodzajow czastek
	int * num_atoms_in_type ;

	int ile_wybielonych;
		int * ile_wyswiecilam; /////// ile kazda czasteczka wyswiecila fotonow do czasu PhBL //DEBUG
		// dla wolnej dyfuzji te liczby sa wieksze - czemu?
		//problem I : zliczam tylko fotowybielone
		//problem II : resetuje przy przejsciu przez granice pudla tylko fotowybielone


	/*** data stored in the trajectory file, see GROMACS documentation for details ***/
	int fin;		// input file pointer
	FILE * fout, * fout2;		// output file pointers (fout2 - for cross-correlation in 2 foci experiments)
	FILE * fp2, * fp3, *fmsd, *fp4;	// file with the avg. number of molecules in the obs. volume, MSD traj. file, number of photons before bleaching
	FILE * fpos; //output for positions (diffusion)
	int magic;
	int natoms;		// TOTAL numbers of atoms (here: fluorophore molecules) of each type
	int dark;
	int step;
	int nstlist;		// how often should we update neighbour list (every nstlist frames)
// 	double czas;
	double prec;
	//double box[3][3];
	enumProfiles profil;
	
	Molecule * mol;
	
//	list<int> * neighbours;	//array of neighbour lists for every molecule (for faster collision checks)
//	list<int> *** sector; //nie moge automatycznie przez NUM_SECTORS, bo to array - dla C++ 'dynamiczny' :/
#ifdef HARD_SPHERES
	int sector [NUM_SEC_X][NUM_SEC_Y][NUM_SEC_Z] [MAX_LIST_LEN];
	int num_sec [NUM_SEC_X][NUM_SEC_Y][NUM_SEC_Z];		//how many molecules are in each sector
	void remove_from_sector (int x, int y, int z, int m);
	void push_to_sector (int x, int y, int z, int m);
	void UpdateMoleculeMovement (const int & i);
#endif //HARD_SPHERES

	int total_frames;
	int steps_per_frame; 	//interesting only when simulating diffusion on our own
	bool simulate_diffusion; //should we simulate the diffusion or read the trajectory from a file?	
	int pinhole;		// Fluorescence cutoff: 0 - off, 1 - on

	enumStates * state;	//excitation states of the molecules HACK aaargh such a bad practice!!!! Move to Molecule class later on!!!!!!!
	int current_photons1, current_photons2;	//the number of photon counts in the current frame

	void Init ();				//initialize the molecule states and auxiliary variables
	void InitCUDA ();			//allocat mem, copy basic stuff etc
	void DeInitCUDA ();
	void ReadFirstFrame ();			//reads the first frame from the molecule trajectory (.xtc) file
	double Distance (rvec i, rvec j);
	
	virtual void ReadGROMACSParameters (char * param_filename);
	void ReadBasicGROMACSParameters (char * param_filename);
	virtual void SimulateDiffusion ();	//simulate the movement of the molecules

	void MoveMolecules_GPU (); //NOTE CUDA only!
	void TryExcitationALL_GPU ();
	
	void ApplyStateChanges ();		//checks for excitation, fluorescence etc.
	void TryExcitationALL ();
	void TryExcitation (int molecule_index);
	void PhotoConversion (int molecule_index);
	void UpdateMoleculePositions ();	//reads the next frame of the molecule trajectory

	double MeanSquareDisplacement ();

	//bool CheckAndAvoidCollisions();		//global colision checking procedure
	int CheckCollision (const int & i);	//check if molecule #i is colliding with anything
	void UndoCollision (int i, int m);
	//void UpdateNeighbourList ();
	
  public:
	double dT;		//trajectory resolution, read from the .xtc file
	int current_frame;

	Fluorescence ();		//use this constructor to simulate diffusion from scratch
//	Fluorescence (char* traj_file); //use this one if diff. trajectory is already calculated NOTE DISABLED
	virtual ~Fluorescence () { 
		if (fin != -1) close_xtc (fin); 
		
		#ifdef ENABLE_GPU
		  DeInitCUDA();
		#endif
		
		if (ile_babli > 0) delete babel;
		
		delete [] mol; delete [] state; delete [] ile_wyswiecilam; delete [] num_atoms_in_type;
/*		if (sector == NULL) return;
		for (int i=0; i<NUM_SECTORS[0]; i++) {
			for (int j=0; j<NUM_SECTORS[1]; j++)
				delete [] sector[i][j];
			delete [] sector[i];
		}	
		delete [] sector;*/
	}

	//Opens molecule diffusion trajectory file & reads the first frame
	//void OpenTrajectoryFile (const char* traj_file);

	//Generates the photon count trajectory
	void CalculateFluorescenceTrajectory (const char * output_file);
};

class FluorRods : public Fluorescence {
  protected:
	int * atoms_per_rod;
	double * ROD_LENGTH;
	int * num_rods;
	
	void Normalize (double v[3]);

	virtual void ReadGROMACSParameters (char * param_filename);
// 	virtual void ReadBasicGROMACSParameters (FILE * fp); //already defined 
	virtual void SimulateDiffusion ();	//simulate the movement of the molecules

  public:
	FluorRods ();
	~FluorRods ();
};

class FluorObstacles : public Fluorescence {
  protected:
// 	void Normalize (double v[3]);

	virtual void ReadGROMACSParameters (char * param_filename);
	virtual void SimulateDiffusion ();	//simulate the movement of the molecules

	double period;
	double relative_size;
	
  public:
	int dim_obstacles;
	double obst_size;
	
	FluorObstacles ();
	~FluorObstacles () {};
};


#endif //_FCS_FLUORESCENCE_H_
