#include "photonHistogram.h"
/*****************************************************************************************************
 * res - resolution (how many frames are grouped in one "fluorescence bin" for calculation purposes) *
 *****************************************************************************************************/
Histogram :: Histogram (int r, const char * input_file) {
 std::ifstream fin (input_file);

 this->res = r; //save the resolution
 
 int tmp;
 bool ok = true;
 while (1) {
   int suma = 0;
   for (int i=0; i<res; i++) {
      fin >> tmp;
      if (!fin.good()) {ok = false; printf ("*** file over!");break;}
      suma += tmp;
   }
   
   if (!ok) break;
   traj.push_back (suma);
   printf ("%d\n",suma);
 }
}

/**********************************************
 * Create a histogram of interdetection times *
 * in frame length units (dT)                 *
 **********************************************/
//num_bins - liczba kanalow w histogramie (a liczba klatek traj.fluorescencji = traj.size() )
void Histogram :: HistogramInterdetectionTimes (int num_bins, const char * output_file) {
 int k=0;  //k = time between photons

 printf ("\nGenerating interdetection times histogram...\n");

 int * T = new int [num_bins+1];
 for (int i=0;i<num_bins+1;i++) T[i]=0;

 for (unsigned int i=0;i<traj.size();i++) {
	assert (traj[i]>=0);	//sanity check - czy nie jest ujemne

	k++;

	if (traj[i]) {
		T[0]+=traj[i]-1;
		if (k>num_bins) T[num_bins] ++;
		else T[k]++;
// 		printf("Adding to bin %d\n",k);
		k=0;
	}
 }

 FILE * fp = fopen (output_file, "wt");
 fprintf (fp, "@title \"Interdetection times histogram\"\n"); 
 for (int i = 0; i <= num_bins; i++) fprintf (fp, "%d %d\n",i,T[i]);

 fclose (fp);
}

/*****************************************************************************************************
 * num_bins - number of output bins in the histogram                                                 *
 *****************************************************************************************************/
void Histogram :: HistogramPhotonCounts (int num_bins, const char * output_file) {
 printf ("\nGenerating photon count histogram...\n");

 int sum = 0;				 		//sum of photon counts in a fluorescence bin (width = res * dT)

 int * P = new int [num_bins+1];
 for (int i=0;i<num_bins+1;i++) P[i]=0;			//create a clear array of photon count numbers

 for (unsigned int i=0;i<traj.size();i++) {
	assert (traj[i]>=0);			//sanity check (DEBUG)
	if (traj[i]< num_bins) P[traj[i]] ++;		//if the number of photons is out of range, count it into the last bin of the histogram
	else P[num_bins] ++;			//else, increase the count value in the adequate histogram bin
 }

 FILE * fp = fopen (output_file, "wt");			//write out results to a file
 fprintf (fp, "@title \"Photon count histogram (resolution: %d)\"\n",res); 
 for (int i = 0; i <= num_bins; i++) fprintf (fp, "%d %d\n",i,P[i]);

 fclose (fp);
}
