#include "correlator.h"
#include "log.h"

/******************************************************************************
 * Constructor, to be used if the fluorescence trajectory is stored in files. *
 * Arguments as above.                                                        *
 ******************************************************************************/
CrossCorrelator:: CrossCorrelator (const char * fluor_file1, const char * fluor_file2, const double step_size, const int num_steps, const int _steps_per_bin) : Correlator (step_size, num_steps, _steps_per_bin) 
{
	LOG ("Creating CrossCorrelator object with step_size = %f [ps], num_steps = %d, steps_per_bin = %d.\n",step_size,num_steps,_steps_per_bin);
	K2 = 0; avg_intensity2=0;
	
	assert (B == num_steps/_steps_per_bin);	//DEBUG
	printf ("Using %d frames.\n",B*_steps_per_bin);
	photons1 = new int [B]; photons2 = new int [B];		//allocate memory for numbers of photon counts in each bin
	
	FILE * fp1, *fp2;
	assert ((fp1 = fopen64 (fluor_file1, "rt")) != NULL);	//NOTE use large file support
	assert ((fp2 = fopen64 (fluor_file2, "rt")) != NULL);	//NOTE use large file support

	for (int i=0; i<B; i++) {
		photons1[i]=0;	photons2[i]=0;		// clear the array
		for (int j=0;j<_steps_per_bin;j++) {
			int phot;
			fscanf (fp1, "%d\n",&phot);
			assert (phot>=0);	// //DEBUG
			photons1[i]+=phot;
 			K+=phot;

			fscanf (fp2, "%d\n",&phot); //read 2nd channel
			assert (phot>=0);	// //DEBUG
			photons2[i]+=phot;
 			K2+=phot;
		}
	}
	fclose(fp1);fclose(fp2);
	avg_intensity = (double)K/(double)B;
	avg_intensity2 = (double)K2/(double)B;
	LOG ("K1= %lld, B= %d, AVG_INT1 = %f\n. K2= %lld, B= %d, AVG_INT2 = %f\n",K,B,avg_intensity,K2,B,avg_intensity2);
}

CrossCorrelator :: ~CrossCorrelator () {
	delete [] photons1; delete [] photons2;
}

double CrossCorrelator :: Calculate_correlation (double time) {

	int k = int (time / bin_size);		//k - delay (bin size units)
	double sum = 0;

	assert (k>=0);

	for (int i=0;i<B-k;i++) {
		if (photons1[i]<0 || photons2[i+k]<0) {	//DEBUG
			printf ("ERROR! ph[%d]= %d, ph[%d] = %d!\n",i,photons1[i],i+k,photons2[i+k]);
		}

		assert (photons1[i]>=0);
		assert (photons2[i+k]>=0);
 		sum += photons1[i]*photons2[i+k];
	}
 	double correlation = sum/((B - k)*(avg_intensity*avg_intensity2)) - 1;	//calculate the value of the correlation
	return correlation;
}
