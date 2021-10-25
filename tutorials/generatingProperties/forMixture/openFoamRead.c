#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define NX 1001  // KSJ : number of internalFields in openFOAM
//#define NV 3 // KSJ : number of data points that you want to extract
#define NV 8 // KSJ : number of data points that you want to extract


char line[4096];
double data[NV*NX];
char casename[1024];
char msg[1024];

int i,j;
int tecfileread();

int KK_RD, KK_SK, II_SK;

int main(int argc, char* argv[])
{
	FILE* fp;
	double *buf;
	int n, m, iflag;


	if (argc!=3){
		printf("Wrong\n");
		return -1;
	}

	strncpy(casename, argv[1], sizeof(casename));

	if(tecfileread(argv[1])){
		printf("error reading tec file: %s\n", argv[1]);
		return -1;
	}

	fp = fopen(argv[2], "wt+");
	if(fp==NULL){
		printf("error opening %s\n", argv[2]);
		return -1;
	}
       
	for(n=0; n<NX; n++){
        for (j=0; j<NV-1; j++){
            fprintf(fp, "%.6e   ", data[j+NV*n]);
        }
            fprintf(fp, "%.6e\n", data[NV-1+NV*n]);                      
    }

//	for(n=0; n<NX*NY*NV; n++){
//	fprintf(fp, "%.6g\n", data[n]);
//	}

	fclose(fp);

	
	return 0;
}

int tecfileread()
{
int i, j, n;
i=0;
double val;
FILE* fp;

    memset(data,0,sizeof(data));
	
char s1[20000]="T p rho mu kappa Cv Cp HE";   
// KSJ : specify the variables that you want to extract. Should be same as NV

char *species= strtok(s1, " ");
char fraction[50];

while (species !=NULL)
{

	sprintf(fraction,"./1e-08/""%s",species);   
	fp= fopen(fraction, "rt");
   
     if(fp==NULL){
	    printf("failed to open %s\n", fp);
    return -1;
        }
	
			 
//skip the header
    while(1){
	if (fgets(line, sizeof(line), fp) ==NULL){
		printf("error reading header\n");
		fclose(fp);
		return -1;
	    }
	if (!strncmp(line, "(",1)) break;
	} 
	
//read data

     sprintf(msg, "reading variable %d...", i);

     for (j=0; j<NX; j++){
        n = fscanf(fp, "%lf", &val);
        if (n!=1){
             printf("error reading data\n");
             fclose(fp);
             return -1;                   
			}
    data[j*NV+i] = val;
    }
	
	i++;
    species =strtok(NULL, " ");
}

  printf("finished reacing data.\n"); 

//try to read more data
n=0;
while ( fscanf(fp, "lf", &val) >0 ){
n++;
}
if(n>0){
	printf("contains more data than expected : %d\n",n);
	sprintf(msg, "contains more data than expected: %d\n", n);
}

	
    fclose(fp);
    return 0;
}

