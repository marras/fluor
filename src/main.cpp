#include "fluorescence.h"
#include "correlator.h"
#include "photonHistogram.h"
#include <string>

// const int RNG_SEED = 12;	//zarodek do generatora liczb losowych
// MTRand rng(RNG_SEED);
MTRand rng;

C_log Log ("log.txt");

enum tag_Options{
 LINEAR = 1,
 LOGARITHMIC,
};

/***************************************************************************************************
 * Generate fluorescence trajectory for molecules trapped in a channel between low-diffusive areas *
 ***************************************************************************************************/
void CalcObstaclesFluorescence() {
   LOG ("Starting simulation of diffusion with obstacles.");
   Log.debug = HEAVY_DEBUG;

    FluorObstacles fluor;
    fluor.CalculateFluorescenceTrajectory ("test.flu");
   LOG ("Fluorescence simulation finished.");
    double TAU_MAX = 0.1*fluor.dT*fluor.current_frame;	//default max delay for correlation function calculation: 0.1 * trajectory length

   LOG ("*Number of frames in fluorescence trajectory: %d, dT = %f",fluor.current_frame,fluor.dT);

   FILE * fpo = fopen ("params.txt","wt");			//write down parameters for the next stages
//Print:Total number of frames, Frame size, Maximum correlation time, default bin size for linear correlation calculations: 20 frames, some extra data
    fprintf (fpo,"NUM_FRAMES %d\ndT %f\nTAU_MAX %f\nBIN_SIZE 1\n\n#Q_t=%f\n#Q_b=%f\n#Epsilon=%e\n#Obstacles=%d^3\n#Obstacle size=%lf\n",
		 fluor.current_frame,fluor.dT,TAU_MAX,Qt,Qb,EPSILON,fluor.dim_obstacles,fluor.obst_size);
   fclose (fpo);
}


/***************************************************************************
 * Generate fluorescence trajectory for rigid, non-rotating diffusing rods *
 ***************************************************************************/
void CalcRodsFluorescence() {
   LOG ("Starting rods diffusion and fluorescence simulation.");
   Log.debug = HEAVY_DEBUG;

   #ifdef ENABLE_GPU
    LOG ("*Simulating rods diffusion with CUDA not implemented. The program will utilize the CPU only.");
   #endif
   
    FluorRods fluor;
    fluor.CalculateFluorescenceTrajectory ("test.flu");
   LOG ("Fluorescence simulation finished.");
    double TAU_MAX = 0.1*fluor.dT*fluor.current_frame;	//default max delay for correlation function calculation: 0.1 * trajectory length

   LOG ("*Number of frames in fluorescence trajectory: %d, dT = %f",fluor.current_frame,fluor.dT);

   FILE * fpo = fopen ("params.txt","wt");			//write down parameters for the next stages
//Print:Total number of frames, Frame size, Maximum correlation time, default bin size for linear correlation calculations: 20 frames, some extra data
    fprintf (fpo,"NUM_FRAMES %d\ndT %f\nTAU_MAX %f\nBIN_SIZE 1\n\n#Q_t=%f\n#Q_b=%f\n#Epsilon=%e\n",
		 fluor.current_frame,fluor.dT,TAU_MAX,Qt,Qb,EPSILON);
   fclose (fpo);
}

/*********************************************************************************
 * Generate fluorescence trajectory and the number of mol. in observation volume *
 *********************************************************************************/
void CalcFluorescence () {
   LOG ("Starting fluorescence simulation.");
   Log.debug = HEAVY_DEBUG;
//    Fluorescence fluor ("trajektoria.xtc");
    Fluorescence fluor;
    fluor.CalculateFluorescenceTrajectory ("test.flu");
   LOG ("Fluorescence simulation finished.");
    double TAU_MAX = 0.1*fluor.dT*fluor.current_frame;	//default max delay for correlation function calculation: 0.1 * trajectory length

   LOG ("*Number of frames in fluorescence trajectory: %d, dT = %f",fluor.current_frame,fluor.dT);

   FILE * fpo = fopen ("params.txt","wt");			//write down parameters for the next stages
//Print:Total number of frames, Frame size, Maximum correlation time, default bin size for linear correlation calculations: 20 frames, some extra data
    fprintf (fpo,"NUM_FRAMES %d\ndT %f\nTAU_MAX %f\nBIN_SIZE 1\n\n#Q_t=%f\n#Q_b=%f\n#Epsilon=%e\n",
		 fluor.current_frame,fluor.dT,TAU_MAX,Qt,Qb,EPSILON);
   fclose (fpo);
}

/*******************************************************************************************************
 * Calculate correlation function, using Simple or PAT correlator, in a Linear or Logarithmic fashion. *
 *******************************************************************************************************/
void CalcCorrelation (int method, char * filename = "test.flu") {
    FILE * fp_params = fopen ("params.txt","rt");
    if (fp_params == NULL) {LOG("!Could not find parameters file!"); return;}

    std::string fname (filename);
    std::string out_name ("corr");
    
    if (fname != "test.flu") //in case of custom filename, change the output name too
	out_name += fname.substr(8,100);

    out_name += ".dat";
    LOG ("*Output correlation to %s",out_name.c_str());


    int num_steps = -1;
    double dT = -1.0, TAU_MAX = -1.0;
    int STEPS_PER_BIN = -1;

    if (fscanf (fp_params,"NUM_FRAMES %d\ndT %lf\nTAU_MAX %lf\nBIN_SIZE %d\n",&num_steps, &dT, &TAU_MAX, &STEPS_PER_BIN) != 4)
     { LOG ("!Error reading parameters file (params.txt)."); exit (2); }
    
    LOG ("*Calculating correlation for %s.\nNUM_FRAMES %d, dT %f, TAU_MAX %f, BIN_SIZE %d\n",filename,num_steps, dT, TAU_MAX, STEPS_PER_BIN);

    FILE * fp;

    if (method == LINEAR) {
//      fp = fopen ("correlationLIN.dat","wt");
     fp = fopen (out_name.c_str(),"wt");
     
     if (STEPS_PER_BIN <= 0) {LOG ("!Need a fixed number of bins for LINEAR method. Currently params.txt defines BIN_SIZE as %d.",STEPS_PER_BIN); exit (54);} 

     LOG ("Calculating correlation (Linear).");
     fprintf (fp, "@title \"G(tau), %d Steps/bin, dT = %f\"\n",STEPS_PER_BIN,dT);

     Correlator * corr1 = new CorrelatorSimple (filename, dT, num_steps, STEPS_PER_BIN);

     for (double time=STEPS_PER_BIN*dT; time <= TAU_MAX; time+= STEPS_PER_BIN*dT) {
	fprintf (fp,"%le %lf\n",time/1e12, corr1->Calculate_correlation (time));

	static int counter = 0;
	if (counter > 100) {printf ("%2.2f%%\n",time/TAU_MAX*100);counter-=100;}
	counter++;
     }
    }

    else if (method == LOGARITHMIC) {
	fp = fopen (out_name.c_str(),"wt");
        LOG ("Calculating correlation (Logarithmic).");
        fprintf (fp, "@title \"G(tau), dT = %f\"\n",dT);

        int counter = 0; //counter=1;

	if (STEPS_PER_BIN <= 0)	{	//STEPS_PER_BIN undefined = varying bin size
	  STEPS_PER_BIN = 1;		//start form bin size = frame size, then increase
	  double time = 0;//STEPS_PER_BIN*dT;
	  Correlator * corr1;

	  while (time <= TAU_MAX) {
	      corr1 = new CorrelatorSimple (filename, dT, num_steps, STEPS_PER_BIN);

	      while (counter < 4) {
		time += STEPS_PER_BIN*dT;
		printf ("STEPS_PER_BIN=%d.\n",STEPS_PER_BIN);

		double result = corr1->Calculate_correlation (time);
		fprintf (fp,"%5.9le %5.9le\n",time/1e12, result);
		printf ("g(%5.9le)=%5.9le\n", time/1e12, result);

		counter++;
	      }
	      
	      delete corr1;
	      counter=0; STEPS_PER_BIN*=2;	//every 4 iterations, the bin size (and time change) is doubled
	  }
	}
	else {				//STEPS_PER_BIN constant, as defined in params.txt
	  LOG ("Calculating LOG correlation for BIN_SIZE = %d.",STEPS_PER_BIN);

	  double time = STEPS_PER_BIN*dT;
	  double TIME_STEP = STEPS_PER_BIN;
	  Correlator * corr1 = new CorrelatorSimple (filename, dT, num_steps, STEPS_PER_BIN);

	  printf ("STEPS_PER_BIN=%d.\n",STEPS_PER_BIN);

	  while (time <= TAU_MAX) {
	      double result = corr1->Calculate_correlation (time);
	      fprintf (fp,"%5.9le %5.9le\n",time/1e12, result);
	      printf ("g(%5.9le)=%5.9le\n", time/1e12, result);

	      counter++;
	      if (counter >= 4) {counter=0; TIME_STEP*=2;}	//every 4 iterations, the correlation time increase is doubled
	      time += TIME_STEP*dT;
	  }
	  delete corr1;
	}
    }
    else {LOG ("!Incorrect method choice!"); return;}		//should never happen // DEBUG
    fclose (fp);
    LOG ("Finished calculating correlation.");
}

/******************************
 * Cross-correlation analysis *
 ******************************/
void CalcCrossCorrelation (char * filename1 = "test.flu",char * filename2 = "test2.flu") {
    FILE * fp_params = fopen ("params.txt","rt");
    if (fp_params == NULL) {LOG("!Could not find parameters file!"); return;}

    std::string fname1 (filename1);
    std::string fname2 (filename2);
    std::string out_name ("Xcorr");
    
    if (fname1 != "test.flu" || fname2 != "test2.flu") //in case of custom filename, change the output name too
	out_name += fname1.substr(8,100) += fname2.substr(8,100);

    out_name += ".dat";
    LOG ("*Output correlation to %s",out_name.c_str());

    int num_steps = -1;
    double dT = -1.0, TAU_MAX = -1.0;
    int STEPS_PER_BIN = -1;

    if (fscanf (fp_params,"NUM_FRAMES %d\ndT %lf\nTAU_MAX %lf\nBIN_SIZE %d\n",&num_steps, &dT, &TAU_MAX, &STEPS_PER_BIN) != 4)
     { LOG ("!Error reading parameters file (params.txt)."); exit (2); }
    
    LOG ("*Calculating correlation for %s & %s.\nNUM_FRAMES %d, dT %f, TAU_MAX %f, BIN_SIZE %d\n",filename1,filename2,num_steps, dT, TAU_MAX, STEPS_PER_BIN);

    FILE * fp;

    fp = fopen (out_name.c_str(),"wt");
    LOG ("Calculating cross correlation (Logarithmic).");
    fprintf (fp, "@title \"G(tau), dT = %f\"\n",dT);

    int counter = 0; //counter=1;

    if (STEPS_PER_BIN <= 0)	{	//STEPS_PER_BIN undefined = varying bin size
      STEPS_PER_BIN = 1;		//start form bin size = frame size, then increase
      double time = STEPS_PER_BIN*dT;
      Correlator * corr1;

      while (time <= TAU_MAX) {
	  corr1 = new CrossCorrelator (filename1, filename2, dT, num_steps, STEPS_PER_BIN);
	  printf ("STEPS_PER_BIN=%d.\n",STEPS_PER_BIN);

	  double result = corr1->Calculate_correlation (time);
	  fprintf (fp,"%5.9le %5.9le\n",time/1e12, result);
	  printf ("g(%5.9le)=%5.9le\n", time/1e12, result);

	  delete corr1;

	  counter++;
	  if (counter >= 4) {counter=0; STEPS_PER_BIN*=2;}	//every 4 iterations, the bin size (and time change) is doubled
	  time += STEPS_PER_BIN*dT;
      }
    }
    else {				//STEPS_PER_BIN constant, as defined in params.txt
      LOG ("Calculating CROSS-correlation for BIN_SIZE = %d.",STEPS_PER_BIN);

      double time = STEPS_PER_BIN*dT;
      double TIME_STEP = STEPS_PER_BIN;
      Correlator * corr1 = new CrossCorrelator (filename1, filename2, dT, num_steps, STEPS_PER_BIN);

      printf ("STEPS_PER_BIN=%d.\n",STEPS_PER_BIN);

      while (time <= TAU_MAX) {
	  double result = corr1->Calculate_correlation (time);
	  fprintf (fp,"%5.9le %5.9le\n",time/1e12, result);
	  printf ("g(%5.9le)=%5.9le\n", time/1e12, result);

	  counter++;
	  if (counter >= 4) {counter=0; TIME_STEP*=2;}	//every 4 iterations, the correlation time increase is doubled
	  time += TIME_STEP*dT;
      }
      delete corr1;
    }

    fclose (fp);
    LOG ("Finished calculating cross-correlation.");
}


void CalcHistograms (int res, int range, char * input = "test.flu") {
   LOG ("Creating histograms.");
    Histogram pch (res, input);
    pch.HistogramInterdetectionTimes (range,"histogramTimes.txt");
    pch.HistogramPhotonCounts (range,"histogramPCH.txt");
}

//----------------------------------- M A I N -------------------------------------//
int main (int argc, char * argv []) {

  if (argc<2 || argc >5)
  { LOG ("!\nIncorrect number of command-line parameters.\nUsage: %s [-f|-a|-b|-h] [filename]",argv[0]); return 1;}

switch (argv[1][1]) {	//drugi znak w pierwszym parametrze  ('f','c', lub 'h') (HACK)
  case 'f':
	CalcFluorescence(); break;
  case 'k':	//kanaliki
	CalcObstaclesFluorescence(); break;
  case 'r':
	CalcRodsFluorescence(); break;
  case 'h':
	if (argc == 4) 	CalcHistograms(atoi(argv[2]), atoi(argv[3]));  
	else if (argc == 5) CalcHistograms(atoi(argv[2]), atoi(argv[3]), argv[4]);  
	else LOG ("!Usage: %s -h binning range. Default: %s -h 200 50 [input_file]\n",argv[0],argv[0]);
	break;
  case 'a':
	CalcCorrelation(LINEAR);  break;
  case 'b':
	if (argc == 3) CalcCorrelation(LOGARITHMIC, argv[2]);  
	else if (argc == 2) CalcCorrelation(LOGARITHMIC);
	else LOG ("!Incorrect number of command-line parameters.");
	break;
  case 'x':	//cross correlation (for 2 foci analysis)
	if (argc == 4) CalcCrossCorrelation(argv[2], argv[3]);  
	else LOG ("!Please supply two trajectory files to calculate cross-correlation.");
	break;
  default:
	LOG ("!The option %s does not stand for any command. Please use -f/k/r, -a/b/x or -h.");
	return 2;
 }

  LOG ("DONE.");
  return 0;
}
