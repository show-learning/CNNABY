#include <stdio.h>
#include <stdlib.h>

void double_2_uint8(double *double_img, unsigned char *uint8_img, int cols, int rows);


int main(int argc, char *argv[])
{
	char *inFile = argv[1];
	char *outFile = argv[2];

	//Upsampler parameters
	int scale = 2;

	//Compressed Assault Cube
	int num = 1; //Number of frames to interpolate
	int inCols = 176; //Width of input (downsampled) video
	int inRows = 144; //Height of input (downsampled) video

	int outCols = inCols*scale;
	int outRows = inRows*scale;

	FILE *inFp, *outFp;

	inFp = fopen(inFile, "rb");
	if (inFp == NULL)
	{
		printf("\n We have null pointer \n");
	}
	outFp = fopen(outFile, "wb");
	if (outFp == NULL)
	{
		printf("\n We have null pointer \n");
	}

	// To read and write each frame in an unsigned character format
	unsigned char *inBuf = (unsigned char *)malloc(inCols*inRows*sizeof(unsigned char));
	unsigned char *outBuf = (unsigned char *)malloc(outCols*outRows*sizeof(unsigned char));
	// To work with each pixel in the range of 0~1
	double *inBuf_tmp = (double *)malloc(inCols*inRows*sizeof(double));
	double *outBuf_tmp = (double *)malloc(outCols*outRows*sizeof(double));
	FILE *Y_Fp;
	Y_Fp = fopen("FSR_Y.txt", "r");
	for (int fcnt = 0; fcnt < num; fcnt++)
	{
		//////// Interpolate each frame using FSRCNN for Y component and simple repitition for U and V components
		// Pointer to obtain value of each tpixel of input frame
		unsigned char *inP = inBuf;
		double *inP_tmp = inBuf_tmp;
		// Pointer to obtain value of each pixel of output frame
		unsigned char *outP = outBuf;
		double *outP_tmp = outBuf_tmp;

		//Y Component
		
		int i, j;
		for (int i = 0; i < outRows * outCols; i++)
		{
			fscanf(Y_Fp, "%lf", &outP_tmp[i]);
		}
		outP_tmp = outBuf_tmp;
		double_2_uint8( outP_tmp, outP, outCols, outRows);

		fwrite(outBuf, sizeof(unsigned char), outCols*outRows, outFp);
		fread(inBuf, sizeof(unsigned char), inCols*inRows, inFp);
		//U Component
		fread(inBuf, sizeof(unsigned char), inCols*inRows / 4, inFp);

		inP = inBuf;
		outP = outBuf;

		for (i = 0; i < inRows / 2; i++)
		for (j = 0; j < inCols / 2; j++) {

			int cnt = 2 * (i * outCols / 2 + j);

			unsigned char x = *inP++;

			*(outP + cnt) = x;
			*(outP + cnt + 1) = x;
			*(outP + cnt + outCols / 2) = x;
			*(outP + cnt + outCols / 2 + 1) = x;

		}

		fwrite(outBuf, sizeof(unsigned char), outCols*outRows / 4, outFp);

		// V Component
		fread(inBuf, sizeof(unsigned char), inCols*inRows / 4, inFp);
		inP = inBuf;
		outP= outBuf;

		for (i = 0; i < inRows / 2; i++)
		for (j = 0; j < inCols / 2; j++) {

			int cnt = 2 * (i*outCols / 2 + j);

			unsigned char x = *inP++;

			*(outP + cnt) = x;
			*(outP + cnt + 1) = x;
			*(outP + cnt + outCols / 2) = x;
			*(outP + cnt + outCols / 2 + 1) = x;

		}

		fwrite(outBuf, sizeof(unsigned char), outCols*outRows / 4, outFp);
		
	}
	free(inBuf);
	inBuf = NULL;
	free(inBuf_tmp);
	inBuf_tmp = NULL;
	free(outBuf);
	outBuf = NULL;
	free(outBuf_tmp);
	outBuf_tmp = NULL;
}
void double_2_uint8(double *double_img, unsigned char *uint8_img, int cols, int rows)
{
	int i, j, cnt, k;

	for (i = 0; i < rows;i++)
	for (j = 0; j < cols; j++)
	{
		cnt = i*cols + j;

		if (*(double_img + cnt) < 0)
			* (uint8_img + cnt) = 0;
		if (*(double_img + cnt) > 255)
			* (uint8_img + cnt) = 255;
			
		k = (int)*(double_img + cnt)+0.5;
		*(uint8_img + cnt) =  k;

		// for (k = 0; k < 255; k++)
		// {
		// 	if (*(double_img + cnt) >= k && *(double_img + cnt) < (k+0.5))
		// 	*(uint8_img + cnt) =  k;

		// 	if (*(double_img + cnt) >= (k+0.5) && *(double_img + cnt) < (k+1))
		// 		*(uint8_img + cnt) = k + 1;
		// }

	}
}