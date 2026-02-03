#ifndef _BWT_H_
#define _BWT_H_

int fc_transform(unsigned char * T, int n, unsigned char * num_indexes, int * indexes);
int fc_untransform(unsigned char * T, int n, int index, unsigned char num_indexes, int * indexes);
#endif