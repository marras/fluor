#include<stdio.h> 

int main (int argc, char ** argv) {
	if (argc != 2) {
		printf ("Wrong number of arguments!\n");
		return 1;
	}

	FILE * fp = fopen (argv[1],"rt");

	long long total = 0;
	int lines = 0;

	int a = 0;

	while (!feof(fp)) {
		fscanf (fp,"%d",&a);
		total += a;
		lines ++;
	}

	lines--;
	total -=a;

	printf ("Srednia z %s: %lld / %d = %lf\n", argv[1],total, lines, (double)total/(double)lines);
	return 0;
}
