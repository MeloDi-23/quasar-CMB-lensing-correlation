#ifndef TABLEUTIL_H
#define TABLEUTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "typedef.h"

void load_table(char *filename,table_struct *table);
void free_table(table_struct *table);
void CheckScale(double *rp_data,int Nobs,table_struct table,int *BinFlag);

#endif

