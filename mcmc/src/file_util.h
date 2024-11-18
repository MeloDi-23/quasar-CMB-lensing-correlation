#ifndef FILEUTIL_H
#define FILEUTIL_H
 
int DataDimen(char *filename);
int ReadData(double *x, double *y, int n, char *filename);
int ReadCovar(double **B, int n, char *filename);
void CheckFile(char *filename);
void FO_error(char *filename);

#endif
