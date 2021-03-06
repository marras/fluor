/* Fluorescence calculation module */
#include <fstream>
#include "fluorescence.h"

double EPSILON = -1.0;			//excitation intensity
double Qf = -1.0, Qt=-1.0, Qb=-1.0;	//Quantum yields of Fluorescence, ISC(Triplet), and Photobleaching
double tauT = -1.0;			//exponential time constant for Triplet->Ground conversion, 5 us = 5*10^6 ps	

double TRIPLET_TRESHOLD = -1;
double BLEACH_TRESHOLD = -1;

/********************************************
 * Opens molecule diffusion trajectory file *
 ********************************************/
void Fluorescence :: Init () {
	for (int i=0;i<natoms;i++)
	   for (int d=0; d<3; d++) mol[i].old_x[d] = mol[i].x[d];

	current_frame = 0;
	ile_wybielonych = 0;
}


/* Simulates both diffusion and fluorescence */
Fluorescence :: Fluorescence () {
	//LOG ("*Fluorescence () constructor called."); //kazda klasa potomna tak wywola tę funkcję
  
	simulate_diffusion = true; pinhole=0;
	current_frame = -1; fin = -1; natoms = 0;
	bOK = true; dark =-1; num_foci=-1;
	F1[0]=-1;F1[1]=-1;F1[2]=-1;
	ile_babli = 0;

	ReadGROMACSParameters ("config.dat");

	for (int i=0; i<natoms; i++) {
		for (int d=0; d<3; d++) {
			mol[i].x[d] = rng.rand (SIZE[d]);
			mol[i].init_x[d] = mol[i].x[d];
		}

		//jesli nie trafilismy we wnetrze babla, to probujemy ponownie
		if (ile_babli > 0) { //TODO change to allow many bubbles
			if (!babel->IsInside(mol[i].x)) {printf ("Retry placing molecule %d.\n",i);i--;}
		}
	}

	Init (); //generuj startowe pozycje czastek itp.

	#ifdef ENABLE_GPU
	  InitCUDA (); //uruchom GPU, wgraj pozycje czastek i zasiej ziarno randomizera GPU
	#endif 
	
	LOG ("Profile: %d",profil);
}

/* Auxilliary functions */

double Fluorescence :: Distance (rvec i, rvec j) {
    return sqrt ( (i[0]-j[0])*(i[0]-j[0]) + (i[1]-j[1])*(i[1]-j[1]) + (i[2]-j[2])*(i[2]-j[2]) );
}

void FluorRods :: Normalize (double v[3]) {
	double length = sqrt (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	for (int i=0;i<3;i++) v[i] /= length;	//make the resulting length equal 1.0
}
inline double abs (double x) {
 return x>=0 ? x : -x;
}

/***********************************************************************************
 * Calculates the fluorescence trajectory, storing the results in the output file. *
 ***********************************************************************************/
void Fluorescence :: CalculateFluorescenceTrajectory (const char * file_out) {
  if (!simulate_diffusion) {assert (fin != -1);}

  TRIPLET_TRESHOLD = Qf + Qt;		//accumulated probability for triplet crossing (auxilliary constant for optimization purposes)
  BLEACH_TRESHOLD = Qf + Qt + Qb;	//accumulated probability for photobleaching (auxilliary constant)

  fout = fopen64 (file_out, "wt");	//Fluorescence trajectory file - NOTE uses large file support
  if (num_foci == 2) fout2 = fopen64 ("test2.flu", "wt");	//Fluorescence trajectory file - NOTE uses large file support TODO stala nazwa test2.flu
  
  fp2 = fopen ("Nv.txt","wt");	//"Number of molecules in observation volume" trajectory file
  fp3 = fopen ("PHB.txt","wt");

  // DEBUG (sort of)
  fp4 = fopen ("photons_to_bleach.txt","wt");////
  fpos = fopen64 ("positions.txt","wt"); //do obliczenia profilu fluorescencji :)
  
#ifdef CALC_MSD
  fmsd = fopen ("msd.dat", "wt");
#endif //CALC_MSD

   try {
     while (1) {
	UpdateMoleculePositions ();	//Read the new frame from the molecule position trajectory file	
	ApplyStateChanges ();		//Check for photophysical events
	current_frame++;
      }
  }
  catch (enumErrorNumbers e) {
	if (e == ERR_BAD_DATA) {printf ("ERROR: Corrupted trajectory file!\n");exit(1);}
	else if (e == EOF) printf ("End of file.\n");
  }

  int num_phb = 0;
  for (int i=0;i<natoms;i++) if (state[i]==MS_PHB) num_phb++;
  printf ("Number of photobleached molecules: %d.\n",num_phb);

  fclose (fout); fclose (fp2); fclose (fp3); //close output files
  if (num_foci == 2) fclose (fout2);
  fclose(fp4);fclose(fpos);//// 
#ifdef CALC_MSD
  fclose (fmsd);
#endif //CALC_MSD

}

/***************************************************************************
 * For each molecule, check whether it undergoes excitation state changes, *
 * and count the emitted photons.                                          *
 ***************************************************************************/
void Fluorescence:: ApplyStateChanges () {
  current_photons1 = 0; current_photons2 = 0;						//set the number of photon counts in the current frame to 0
  int ile_w_obj=0;									//the number of molecules in the observation volume

#ifndef ENABLE_GPU
  TryExcitationALL ();						//check for excitation
#else //ENABLE_GPU
  TryExcitationALL_GPU ();					//check for excitation (using CUDA)
#endif //ENABLE_GPU

  for (int m = 0; m < natoms; m++) {
	if (mol[m].x[0]>(F1[0] - OBS_VOL_XY/2) && mol[m].x[0]<(F1[0] + OBS_VOL_XY/2) &&
	    mol[m].x[1]>(F1[1]- OBS_VOL_XY/2) && mol[m].x[1]<(F1[1] + OBS_VOL_XY/2) &&
	    mol[m].x[2]>(F1[2] - OBS_VOL_Z/2) && mol[m].x[2]<(F1[2] + OBS_VOL_Z/2)) ile_w_obj++;	//check the number of molecules in the observation volume

	if (state[m] == MS_EXC_1 || state[m] == MS_EXC_2) PhotoConversion (m);		//if the molecule is in the excited state, check for ISC, fluorescence, or photobleaching
	else if (state[m] == MS_TRI)	{						//check for phosphorescence //NOTE teraz mamy ELSE!
		double r3 = rng.rand();
		if (r3 < (1/tauT)*dT) {
			state[m] = MS_GROUND;
		}
	}
  }

  fprintf (fout, "%d\n",current_photons1);						//print the number of photon counts to the .flu file
  if (num_foci==2) fprintf (fout2, "%d\n",current_photons2);				//print the number of photon counts to the second .flu file
  fprintf (fp2, "%d\n",ile_w_obj);
//   fprintf (fp3, "%d\n",ile_wybielonych);
  
  static int frame_counter = 0;								//every 100 frames, write output to logfile
  if (frame_counter >= PROGRESS_IND_FREQ)
	{
	 frame_counter -= PROGRESS_IND_FREQ; LOG ("Frame %d: %d photons in channel 1.",current_frame,current_photons1);
	 printf ("\n%2.2f%%",(double)current_frame/total_frames*100.0);
	#ifdef CALC_MSD
	 fprintf (fmsd, "%d %le\n", current_frame, MeanSquareDisplacement()/1e18);
	 fflush (fmsd);
	#endif //CALC_MSD
	} 
  frame_counter ++;
} 
	
/***************************************************************************
 * This function checks whether a molecule will absorp a photon and switch *
 * to the excited state, using a Gaussian excitation profile.              *
 ***************************************************************************/
void Fluorescence:: TryExcitationALL () {
    for (int m = 0; m < natoms; m++) if (state[m] == MS_GROUND) TryExcitation (m);
}

void Fluorescence:: TryExcitation (int m) {						//m - molecule index
  assert (state[m] == MS_GROUND);	// DEBUG

/*  double prob_ex = EPSILON *dT*exp(-((mol[m].x[0]-F1[0])*(mol[m].x[0]-F1[0])			//probability of excitation (Gaussian profile) [Dix]
					+(mol[m].x[1]-F1[1])*(mol[m].x[1]-F1[1])
					+(mol[m].x[2]-F1[2])*(mol[m].x[2]-F1[2])/SQR_KAPPA ) / (2*WXY*WXY));*/
 
 //PROFIL GAUSSA
 if (profil == PR_GAUSS) {
  double prob_ex = EPSILON *dT*exp(-2*((mol[m].x[0]-F1[0])*(mol[m].x[0]-F1[0])			//probability of excitation (Gaussian profile) [Winkler] - chyba poprawny
					+(mol[m].x[1]-F1[1])*(mol[m].x[1]-F1[1])
					+(mol[m].x[2]-F1[2])*(mol[m].x[2]-F1[2])/SQR_KAPPA ) / (WXY*WXY));
  double r1 = rng.rand();
  if (r1 < prob_ex) state[m] = MS_EXC_1;
  
  if (num_foci==2) {	//dual-focus setup
    prob_ex = EPSILON *dT*exp(-2*((mol[m].x[0]-F2[0])*(mol[m].x[0]-F2[0])			//[Winkler]
					+(mol[m].x[1]-F2[1])*(mol[m].x[1]-F2[1])
					+(mol[m].x[2]-F2[2])*(mol[m].x[2]-F2[2])/SQR_KAPPA ) / (WXY*WXY));
    r1 = rng.rand();
    if (r1 < prob_ex && (state[m] == MS_GROUND || (state[m] == MS_EXC_1 && rng.rand() > 0.5))) state[m] = MS_EXC_2; //jesli oba wzbudzily, to szansa na zliczenie w 1 lub 2 idzie pol na pol.
  } 
 }

 //PROFIL PROSTOKATNY (czyli "ostro obcięta" elipsoida - tylko do porownania)
 else if (profil == PR_SQUARE) {
   if ( (mol[m].x[0]-F1[0])*(mol[m].x[0]-F1[0]) + (mol[m].x[1]-F1[1])*(mol[m].x[1]-F1[1]) + (mol[m].x[2]-F1[2])*(mol[m].x[2]-F1[2])/SQR_KAPPA < WXY*WXY) state[m] = MS_EXC_1;
   
   if (num_foci == 2) { //dual-focus
     if ( (mol[m].x[0]-F2[0])*(mol[m].x[0]-F2[0]) + (mol[m].x[1]-F2[1])*(mol[m].x[1]-F2[1]) + (mol[m].x[2]-F2[2])*(mol[m].x[2]-F2[2])/SQR_KAPPA < WXY*WXY) 
       if (state[m] == MS_GROUND || (state[m] == MS_EXC_1 && rng.rand() > 0.5)) state[m] = MS_EXC_2; //jesli oba wzbudzily, to szansa na zliczenie w 1 lub 2 idzie pol na pol.
   }
 }
  
 else LOG ("!Profile not defined!");
}

/****************************************************************************************
 * On reaching the excited state, molecules undergo one of the possible photophysical   *
 * conversions, with respective probabilities given by the quantum yields Qf, Qt and Qb.*
 ****************************************************************************************/
void Fluorescence:: PhotoConversion (int m) {
  assert (state[m] == MS_EXC_1 || state[m] == MS_EXC_2);
  assert (Qf + Qt + Qb <= 1.0);

  double r2 = rng.rand();

  if (r2 < Qf) {	//Fluorescence
	if (pinhole) {	//only count the photons from inside confocal volume
		if (num_foci == 2) LOG ("*pinhole for 2 foci not implemented!");
		if ((mol[m].x[0]-F1[0])*(mol[m].x[0]-F1[0])/WXY/WXY + 
		   (mol[m].x[1]-F1[1])*(mol[m].x[1]-F1[1])/WXY/WXY + 
	  	   (mol[m].x[2]-F1[2])*(mol[m].x[2]-F1[2])/WZ/WZ <= 1) {current_photons1++;}  
	}
	else {
		if (state[m]==MS_EXC_1) current_photons1++;	//if pinhole is not taken into account, we get all the fluorescence signal
		else if (state[m]==MS_EXC_2) current_photons2++; //jesli wzbudzil drugi laser, to foton wpada do drugiego detektora (NOTE no crosstalk!)
	}

	state[m] = MS_GROUND;
	ile_wyswiecilam[m]++;
// 	fprintf (fpos, "%lf %lf %lf\n",mol[m].x[0],mol[m].x[1],mol[m].x[2]); //Changed: wypisz do pliku positions.txt miejsca, w ktorych nastapila fluorescencja
  }
  else if (r2 < TRIPLET_TRESHOLD) {	//Intersystem Crossing, TRIPLET_TRESHOLD = Qf + Qt
	state[m] = MS_TRI;
  }
  else if (r2 < BLEACH_TRESHOLD) {	//Photobleaching, BLEACH_TRESHOLD = Qf + Qt + Qb
	state[m] = MS_PHB;

	mol[m].old_x[0]= mol[m].x[0];
	mol[m].old_x[1]= mol[m].x[1];
	mol[m].old_x[2]= mol[m].x[2];
	ile_wybielonych ++;

	fprintf(fp4,"%d\n",ile_wyswiecilam[m]);
  }
  else {				//Relax
	state[m] = MS_GROUND;
  }
}

/************************************************************************
 * Read molecule positions from the next frame of diffusion trajectory. *
 ************************************************************************/
void Fluorescence:: UpdateMoleculePositions () {

	SimulateDiffusion ();	//NOTE reading trajectory is DISABLED
//	if (simulate_diffusion) SimulateDiffusion ();		//if we aren't using a diffusion trajectory file, simulate movement
//	else read_next_xtc(fin, natoms, &step, &czas, box, x, &prec,&bOK);

	//Regenerate the photobleached molecules which moved crossed the box boundary
	for (int i=0;i<natoms;i++) {
		if (state[i]!=MS_PHB) continue;		//tylko dla fotowybielonych

		if ( (abs(mol[i].old_x[0] - mol[i].x[0]) >= SIZE[0]/2) ||
		     (abs(mol[i].old_x[1] - mol[i].x[1]) >= SIZE[1]/2) ||
		     (abs(mol[i].old_x[2] - mol[i].x[2]) >= SIZE[2]/2) ) {
			ile_wybielonych --;
			ile_wyswiecilam[i] = 0;
			state[i]=MS_GROUND;
		}

		mol[i].old_x[0]= mol[i].x[0];
		mol[i].old_x[1]= mol[i].x[1];
		mol[i].old_x[2]= mol[i].x[2];
	}

	if (!bOK) throw ERR_BAD_DATA;
	if (current_frame >= total_frames) throw ERR_EOF;
}

/* Starts fluorescence simulation basing on a movement trajectory */ //NOTE DISABLED! to use, convert to Molecule class
// Fluorescence :: Fluorescence (char * traj_file) {
// 	if (traj_file == NULL) {
// 		LOG("*Fluorescence::Fluorescence(): NULL trajectory given, skipping contructor.");
// 		return; //DIRTY HACK! just an option for the derived classes to skip this constructor ;)
// 	}
// 
//      simulate_diffusion = false; pinhole=0;
// 
// 	fin = open_xtc(traj_file,"r");
// 	read_first_(fin, &natoms, &step, &czas, box, &x,&prec, &bOK);
// 
// 	ReadGROMACSParameters ("config.dat");
// 
// 	Init ();
// }
