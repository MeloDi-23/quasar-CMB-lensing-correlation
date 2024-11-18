#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*** determine the dimension of the data vector ***/
int DataDimen(char *filename)
{
 FILE   *fp;
 int    Ndimen;
 double x,y,err;

 if(!(fp=fopen(filename,"rt"))) {
   printf("cannot open file %s\n",filename);
   exit(-1);
 }
 fscanf(fp,"%lf %lf %lf",&x,&y,&err);
 Ndimen=0;
 while(!(feof(fp))) {
    Ndimen=Ndimen+1;
    fscanf(fp,"%lf %lf %lf",&x,&y,&err);
 }
 fclose(fp);

 return Ndimen;
}

/*** read in data                              ***/
int ReadData(double *x, double *y, int n, char *filename)
{
 FILE   *fp;
 int    i;
 double err;

 if(!(fp=fopen(filename,"rt"))) {
   printf("cannot open file %s\n",filename);
   exit(-1);
 }
 for(i=1;i<=n;i++) fscanf(fp,"%lf %lf %lf",&(x[i]),&(y[i]),&err);
 fclose(fp);

 return 0;  
}

/*** read in covar matrix                      ***/
int ReadCovar(double **B, int n, char *filename)
{
 FILE   *fp;
 int    i,j;

 if(!(fp=fopen(filename,"rt"))) {
   printf("cannot open file %s\n",filename);
   exit(-1);
 }
 for(i=1;i<=n;i++)
    for(j=1;j<=n;j++) 
       fscanf(fp,"%lf",&(B[i][j]));
 for(i=1;i<=n;i++)
    for(j=i+1;j<=n;j++)
       if(B[i][j]!=B[j][i]) 
         printf("Warning: asymmetric covar in %s\n",filename);
 fclose(fp);
 
 return 0;
}

void CheckFile(char *filename)
{
 FILE *fp;

 if(!(fp=fopen(filename,"rt"))) { 
   printf("cannot open file %s\n",filename);
   exit(-1);
 }
 fclose(fp);
}

void FO_error(char *filename)
{
  printf("Cannot open file %s\n",filename);
  exit(-1);
}

