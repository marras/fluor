#include "correlator.h"
#include "log.h"

/****************************************************************************************
 * Constructor, to be used if the fluorescence trajectory is given directly as int[].   *
 *   Parameters:                                                                        *
 * . int* traj - fluorescence trajectory (array of photon counts),                      *
 * . const double step_size - trajectory resolution in [ps],                            *
 * . const int num_steps: length of the trajectory in step_size units,                  *
 * . const int _steps_per_bin: binning resolution for correlation function calculation. *
 ****************************************************************************************/
CorrelatorSimple :: CorrelatorSimple (int * traj, const double step_size, const int num_steps, 
	 const int _steps_per_bin) : Correlator (step_size, num_steps, _steps_per_bin) {

	photons = new int [num_steps/_steps_per_bin];
	for (int i=0;i<B;i++) {
		photons[i]=0;
		for (int j=0;j<_steps_per_bin;j++) {
			photons[i]+=traj[i*_steps_per_bin+j];
 			K+=traj[i*_steps_per_bin+j];
		}
	}
	avg_intensity = (double)K/(double)B;
	LOG ("K= %lld, B= %d, AVG_INT = %f\n",K,B,avg_intensity);
}

/*******************************************************************************
 * Constructor, to be used if the fluorescence trajectory is stored in a file. *
 * Arguments as above.                                                         *
 *******************************************************************************/
CorrelatorSimple:: CorrelatorSimple (const char * fluor_file, const double step_size, const int num_steps, const int _steps_per_bin) : Correlator (step_size, num_steps, _steps_per_bin) 
{
	LOG ("Creating correlatorSIMPLE object with step_size = %f [ps], num_steps = %d, steps_per_bin = %d.\n",step_size,num_steps,_steps_per_bin);

	assert (B == num_steps/_steps_per_bin);	//DEBUG
	printf ("Using %d frames.\n",B*_steps_per_bin);
	photons = new int [B];			//allocate memory for numbers of photon counts in each bin
	FILE * fp;
	assert ((fp = fopen64 (fluor_file, "rt")) != NULL);	//NOTE use large file support

	for (int i=0; i<B; i++) {
		photons[i]=0;			// clear the array
		for (int j=0;j<_steps_per_bin;j++) {
			int phot;
			fscanf (fp, "%d\n",&phot);
			assert (phot>=0);	// //DEBUG
			photons[i]+=phot;
			assert (photons[i]>=0); // //DEBUG
 			K+=phot;
		 }
	}
	fclose(fp);
	avg_intensity = (double)K/(double)B;
	LOG ("K= %lld, B= %d, AVG_INT = %f\n",K,B,avg_intensity);
}

CorrelatorSimple :: ~CorrelatorSimple () {
	delete [] photons;
}

double CorrelatorSimple :: Calculate_correlation (double time) {

	int k = int (time / bin_size);		//k - delay (bin size units)
	double sum = 0;

	assert (k>=0);

	for (int i=0;i<B-k;i++) {
		if (photons[i]<0 || photons[i+k]<0) {	//DEBUG
			printf ("ERROR! ph[%d]= %d, ph[%d] = %d!\n",i,photons[i],i+k,photons[i+k]);
		}

		assert (photons[i]>=0);
		assert (photons[i+k]>=0);
 		sum += photons[i]*photons[i+k];
	}
 	double correlation = sum/((B - k)*(avg_intensity*avg_intensity)) - 1;	//calculate the value of the correlation
	return correlation;
}

/********************
 * EXAMPLE OF USAGE *
 ********************/
/*
const int NUM_STEPS = 1500;	//liczba krokow trajektorii fluorescencji (jedn. dT)
const double LAMBDA = 0.1;	//f(t) = lambda * e^( -lambda*t )
const int dT	= 1;		//dlugosc kroku do obliczania F(t)	(step size)
const int DT	= 2;		//szerokosc kanalu do zliczania G(tau)  (bin size)
const double TAU = 6;		//odstep czasowy, dla jakiego obliczamy funkcje korelacji

void Gen_random_times (const int num_steps, const double step_size, int * & dist) {
 	dist = new int [num_steps];
	memset (dist,0,num_steps);
	double t = 0.0f;
	int i=0;

	FILE * fp = fopen ("times.txt","wt");

	while (i<num_steps) {
		//printf ("i=%d, t= %f, dist[%d]=%d\n",i,t,i,dist[i]); //////
		if (t>step_size) { i++; t-=step_size; continue; }

		double rnd = rng.rand();
		t+= (-log(1-rnd))/LAMBDA;
		//t += 10 + rnd/5;	//f. okresowa!
		fprintf (fp,"%f\n",t);
		//t += rnd;	//jednostajna?
		dist[i]++;
	}
	fclose (fp);
}

int main () {
  srand (RAND_SEED);

  int * intensity=NULL;					//trajektoria intensywnosci fluorescencji
  Gen_random_times (NUM_STEPS, dT, intensity);		//zapisana z rozdzielczoscia dT

  for (int i=0;i<NUM_STEPS;i++) printf ("F(%d) = %d \n",i,intensity[i]); ////////

  CorrelatorSimple corr (intensity, dT, NUM_STEPS, DT/dT); 	//create a correlator object

  FILE * fp = fopen ("correlation.txt","wt");

//  fprintf (fp,"%d %f\n",0, corr.Calculate_correlation (0));
    for (int k=0; k<NUM_STEPS*dT/DT*3/4; k++) {	//// * 3/4 - nie bierzemy calej traj. bo pod koniec sa szumy
      fprintf (fp,"%d %f\n",k, corr.Calculate_correlation (k));
    }

  fclose (fp);

 return 0;
}
*/
