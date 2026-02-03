#ifndef _FC_H_
#define _FC_H_

int fcinit();
 int fc_encode(const unsigned char * input, unsigned char * output, int inputSize, int outputSize);
 int fc_decode(const unsigned char * input, unsigned char * output);

#endif