#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <ENCRYPTO_utils/crypto/crypto.h>
#include <abycore/aby/abyparty.h>
#include <abycore/circuit/share.h>
#include <abycore/circuit/booleancircuits.h>
#include <abycore/circuit/arithmeticcircuits.h>
#include <abycore/sharing/sharing.h>
#include<bits/stdc++.h>
#include <aby.class.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <string>
#include <netinet/in.h>

using namespace std;
ABYParty* party;
BooleanCircuit* bc;
ArithmeticCircuit* ac;

void buildparty(std::string address,std::string p,std::string r,std::string path){
    //default argument
    uint16_t port=(uint16_t)std::stoi(p);
    e_role role=(r=="SERVER")?SERVER:CLIENT;
    uint32_t secparam=128;
    seclvl seclvl = get_sec_lvl(secparam);
    e_mt_gen_alg mt_alg = MT_OT;
    //std::string circuit_dir = "/home/wei/ABY/bin/circ/";
    uint32_t nthreads=1,bitlen=64,out_n,out_l;
    //build party
    party=new ABYParty(role,address,port,seclvl,bitlen,nthreads,mt_alg, 100000, path);
    std::vector<Sharing*>& sharings = party->GetSharings();
    bc = (BooleanCircuit*) sharings[S_BOOL]->GetCircuitBuildRoutine();
    ac = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine(); 
}

void pad_image(double *img, double *img_pad, int rows, int cols, int padsize)
{ // This function receives an image and paddes its border in a replicative manner
	int cols_pad = cols + 2 * padsize;
	int rows_pad = rows + 2 * padsize;
	int i, j, k, cnt, cnt_pad, k1, k2;
	// Centeral pixels
	for (i = padsize; i < rows_pad - padsize; i++)
	for (j = padsize; j < cols_pad - padsize; j++)
	{
		cnt_pad = i * cols_pad + j;
		cnt = (i - padsize)*(cols) + j - padsize;
		double x = *(img + cnt);
		*(img_pad + cnt_pad) = x;
	}
	// Top and Bottom Rows
	for (j = padsize; j < cols_pad - padsize; j++)
	for (k = 0; k < padsize; k++)
	{
		// Top Rows 
		cnt_pad = j + k*cols_pad;
		cnt = j - padsize;
		*(img_pad + cnt_pad) = *(img + cnt);
		// Bottom Rows
		cnt_pad = j + (rows_pad - 1 - k)* cols_pad;
		cnt = (j - padsize) + (rows - 1)*cols;
		*(img_pad + cnt_pad) = *(img + cnt);
	}
	// Left and Right Columns
	for (i = padsize; i < rows_pad - padsize; i++)
	for (k = 0; k < padsize; k++)
	{
		// Left Columns
		cnt = (i - padsize)*cols;
		cnt_pad = i*cols_pad + k;
		*(img_pad + cnt_pad) = *(img + cnt);
		// Right Columns
		cnt = (i - padsize)*cols + cols - 1;
		cnt_pad = i*cols_pad + cols_pad - 1 - k;
		*(img_pad + cnt_pad) = *(img + cnt);
	}
	// Corner Pixels
	for (k1 = 0; k1 < padsize; k1++)
	for (k2 = 0; k2 < padsize; k2++)
	{
		// Upper Left Corner
		cnt_pad = k1*cols_pad + k2;
		*(img_pad + cnt_pad) = *(img);
		// Upper Right Corner
		cnt_pad = k1*cols_pad + cols_pad - 1 - k2;
		*(img_pad + cnt_pad) = *(img + cols - 1);
		// Lower Left Corner
		cnt_pad = (rows_pad - 1 - k1)*cols_pad + k2;
		*(img_pad + cnt_pad) = *(img + (rows - 1)*cols);
		// Lower Right Corner
		cnt_pad = (rows_pad - 1 - k1)*cols_pad + cols_pad - 1 - k2;
		*(img_pad + cnt_pad) = *(img + (rows - 1)*cols + cols - 1);
	}
}

void imfilter(double *img, double *kernel, double *img_fltr, int rows, int cols, int padsize)
{
	// img_pad is the pointer to padded image
	// kernel is the pointer to the kernel which used for convolution
	// img_fltr is the pointer to the filtered image by applying convolution
	int cols_pad = cols + 2 * padsize;
	int rows_pad = rows + 2 * padsize;
	int i, j, cnt, cnt_pad, cnt_krnl, k1, k2;
	double sum;
	double *img_pad = (double *)malloc(rows_pad * cols_pad * sizeof(double));
	pad_image(img, img_pad, rows, cols, padsize);
	for (i = padsize; i < rows_pad - padsize; i++)
	for (j = padsize; j < cols_pad - padsize; j++)
	{
		cnt = (i - padsize)*cols + (j - padsize); // counter which shows current pixel in filtered image (central pixel in convolution window)
		sum = 0;
		cnt_krnl = 0; // counter which determines kernel elements
		for (k1 = -padsize; k1 <= padsize; k1++)
		for (k2 = -padsize; k2 <= padsize; k2++)
		{
			cnt_pad = (i + k1)*cols_pad + j + k2; // counter which shows each neighbouring pixel of padded image used for convolution with kernel
			sum = sum + (*(img_pad + cnt_pad))*(*(kernel + cnt_krnl));
			cnt_krnl++;
		}
		*(img_fltr + cnt) = sum;
	}
	free(img_pad);
	img_pad = NULL;
}

double Max(double a, double b)
{
	double c;
	c = a > b ? a : b;
	return c;
}

double Min(double a, double b)
{
	double c;
	c = a > b ? b : a;
	return c;
}

void PReLU(double *img_fltr,int rows, int cols, double bias, double prelu_coeff)
{
	int cnt = 0;
	for (int i = 0; i < rows;i++)
	for (int j = 0; j < cols; j++)
	{
		cnt = i*cols + j;
		*(img_fltr + cnt) = Max(*(img_fltr + cnt) + bias, 0) + prelu_coeff * Min(*(img_fltr + cnt) + bias, 0);
	}
}

void imadd(double *img_fltr_sum, double *img_fltr_crnt, int cols, int rows)
{
	// *img_fltr_crnt ==> pointer to current feature map
	// *img_fltr_sum ==> pointer to the cumulutive feature map

	int cnt = 0;
	for (int i = 0; i < rows;i++)
	for (int j = 0; j < cols; j++)
	{
		cnt = i*cols + j;
		*(img_fltr_sum + cnt) = *(img_fltr_sum + cnt) + *(img_fltr_crnt + cnt);
	}
}

void deconv(double *img_input, double *img_output, double *kernel, int cols, int rows, int stride)
{
	int border = 1;
	int fsize = 9;
	int rows_pad = rows + 2 * border;
	int cols_pad = cols + 2 * border;
	double *img_input_padded = (double *)malloc(rows_pad * cols_pad * sizeof(double));
	pad_image(img_input, img_input_padded, rows, cols, border);
	
	int rows_out_pad = rows_pad * stride;
	int cols_out_pad = cols_pad * stride;
	double *img_output_tmp = (double *)calloc((rows_out_pad + fsize - 1)* (cols_out_pad + fsize - 1), sizeof(double));
	double *kernel_modif = (double *)malloc(fsize * fsize * sizeof(double));

	int idx, idy;
	for (int i = 0; i < rows_pad; i++)
	for (int j = 0; j < cols_pad; j++)
	{
		int cnt_img = i*cols_pad + j;
		idx = i*stride;
		idy = j*stride;
		int cnt_img_output = idx*(cols_out_pad + fsize - 1) + idy; // (idx,idy) coordinate in temporal output image
		int cnt_kernel = 0;
		for (int k_r = 0; k_r < fsize; k_r++)
		{
		for (int k_c = 0; k_c < fsize; k_c++)
		{
			cnt_kernel = k_r*fsize + k_c;
			*(kernel_modif + cnt_kernel) = (*(kernel + cnt_kernel))*(*(img_input_padded + cnt_img));
			*(img_output_tmp + cnt_img_output + k_c) = *(img_output_tmp + cnt_img_output + k_c) + *(kernel_modif + cnt_kernel);
			
		}
		cnt_img_output = cnt_img_output + (cols_out_pad + fsize - 1);
	    }
		
	}

	int rows_out = rows*stride;
	int cols_out = cols*stride;

	for (int i = 0; i < rows_out; i++)
	for (int j = 0; j < cols_out; j++)
	{
		int i_tmp = i + ((fsize + 1) / 2) + stride*border - 1;
		int j_tmp = j + ((fsize + 1) / 2) + stride*border - 1;
		int cnt_img_out = i*cols_out + j;
		int cnt_img_out_tmp = i_tmp*(cols_out_pad + fsize - 1) + j_tmp; // (cols-pad+fsize-1) is the number of columns in the img_out_tmp
		*(img_output + cnt_img_out) = *(img_output_tmp + cnt_img_out_tmp);

	}

	free(img_input_padded); img_input_padded = NULL;
	free(img_output_tmp); img_output_tmp = NULL;
	free(kernel_modif); kernel_modif = NULL;
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
	}
}

int main(int argc,char** argv){
	//Upsampler parameters
	int scale = 2;

	//Compressed Assault Cube
	int inCols = 176; //Width of input (downsampled) video
	int inRows = 144; //Height of input (downsampled) video

	int rows=144,cols=176;

	int outCols = inCols*scale;
	int outRows = inRows*scale;

	FILE *inFp, *outFp;

	inFp = fopen("akiyo_qcif.yuv", "rb");
	if (inFp == NULL)
	{
		printf("\n We have null pointer \n");
	}
	outFp = fopen("akiyo_test.yuv", "wb");
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

	unsigned char *inP = inBuf;
	double *inP_tmp = inBuf_tmp;
	// Pointer to obtain value of each pixel of output frame
	unsigned char *outP = outBuf;
	double *outP_tmp = outBuf_tmp;

	//Y Component
	fread(inBuf, sizeof(unsigned char), inCols*inRows, inFp);
	int i, j;
	*inP=*(inP+176*144);

	FILE *L5Fp;
    L5Fp = fopen("L5_output.txt", "r");
    int num_filters5 = 12;
    double *img_fltr_5_tmp = (double *)calloc(rows * cols * num_filters5 , sizeof(double));
    double *img_fltr_5 = (double *)calloc(rows * cols * num_filters5 , sizeof(double));
    double *img_fltr_p5 = img_fltr_5;
	std:string address(argv[1]),port(argv[2]),r(argv[3]),path(argv[4]);
	buildparty(address,port,r,path);
	SHR *s_img_fltr_5 = (SHR *)malloc(rows * cols * num_filters5*sizeof(SHR));
	e_role role=(r=="SERVER")?SERVER:CLIENT;

    for (int i = 0; i < rows * cols * num_filters5; i++)
	{
		fscanf(L5Fp, "%lf", (img_fltr_5_tmp+i));
		s_img_fltr_5[i] = SHR((uint64_t *)&img_fltr_5_tmp[i], bc, party);
        img_fltr_5[i] = s_img_fltr_5[i].check();
	}
    free(img_fltr_5_tmp);
    /////////// Layer6
	// Reading weights of 6th layer
	FILE *weights_layer6_ptr;
	weights_layer6_ptr = fopen("weights_layer6.txt", "r");
	if (weights_layer6_ptr == NULL) { printf("Error in the reading weights of 6th layer\n"); };
	double weights_layer6[1296];
	for (int i = 0; i < 1296; i++)
	{
		fscanf(weights_layer6_ptr, "%lf", &weights_layer6[i]);
	}
	fclose(weights_layer6_ptr);
	// Reading biases of 6th layer
	FILE *biases_layer6_ptr;
	biases_layer6_ptr = fopen("biasess_layer6.txt", "r");
	if (biases_layer6_ptr == NULL) { printf("Error in the reading biases of 6th layer\n"); };
	double biases_layer6[12];
	for (int i = 0; i < 12; i++)
	{
		fscanf(biases_layer6_ptr, "%lf", &biases_layer6[i]);
	}
	fclose(biases_layer6_ptr);
	// Other parameters
	int filtersize6 = 9; //3X3
	int patchsize6 = 3;
	int padsize6 = (patchsize6 - 1) / 2;
	int num_filters6 = 12;
	int num_channels6 = 12;
	double prelu_coeff_layer6 = 0.7806;
	// Convolution
	double *img_fltr_6 = (double *)calloc(rows * cols * num_filters6 , sizeof(double));
	double *kernel6 = (double *)malloc(filtersize6*sizeof(double));
	double *img_fltr_p6 = img_fltr_6; // Pointer to img_fltr6
	double *img_fltr_6_tmp = (double *)malloc(rows * cols * sizeof(double));

	int cnt_weight = 0;

	for (int i = 0; i < num_filters6; i++)
	{
		img_fltr_p5 = img_fltr_5;
		for (int j = 0; j < num_channels6; j++)
		{
			// reading corresponding weights to kernel
			for (int cnt_kernel = 0; cnt_kernel < filtersize6; cnt_kernel++)
			{
				*(kernel6 + cnt_kernel) = weights_layer6[cnt_weight + cnt_kernel];
			}

			imfilter(img_fltr_p5, kernel6, img_fltr_6_tmp, rows, cols, padsize6);
			imadd(img_fltr_p6, img_fltr_6_tmp, cols, rows);

			cnt_weight = cnt_weight + filtersize6;
			img_fltr_p5 = img_fltr_p5 + rows*cols;
		}
		double bias_tmp = biases_layer6[i];
		PReLU(img_fltr_p6, rows, cols, bias_tmp, prelu_coeff_layer6);
		img_fltr_p6 = img_fltr_p6 + rows*cols;
	}

	free(img_fltr_5);
	img_fltr_5 = NULL;

	/////////// Layer7
	// Reading weights of 7th layer
	FILE *weights_layer7_ptr;
	weights_layer7_ptr = fopen("weights_layer7.txt", "r");
	if (weights_layer7_ptr == NULL) { printf("Error in the reading weights of 7th layer\n"); };
	double weights_layer7[672];
	for (int i = 0; i < 672; i++)
	{
		fscanf(weights_layer7_ptr, "%lf", &weights_layer7[i]);
	}
	fclose(weights_layer7_ptr);
	// Reading biases of 7th layer
	FILE *biases_layer7_ptr;
	biases_layer7_ptr = fopen("biasess_layer7.txt", "r");
	if (biases_layer7_ptr == NULL) { printf("Error in the reading biases of 7th layer\n"); };
	double biases_layer7[56];
	for (int i = 0; i < 56; i++)
	{
		fscanf(biases_layer7_ptr, "%lf", &biases_layer7[i]);
	}
	fclose(biases_layer7_ptr);
	// Other parameters
	int filtersize7 = 1; //1X1
	int patchsize7 = 1;
	int padsize7 = (patchsize7 - 1) / 2;
	int num_filters7 = 56;
	int num_channels7 = 12;
	double prelu_coeff_layer7 = 0.0087;
	// Convolution
	double *img_fltr_7 = (double *)calloc(rows * cols * num_filters7 , sizeof(double));
	double *kernel7 = (double *)malloc(filtersize7*sizeof(double));
	double *img_fltr_p7 = img_fltr_7; // Pointer to img_fltr7
	double *img_fltr_7_tmp = (double *)malloc(rows * cols * sizeof(double));

	cnt_weight = 0;

	for (int i = 0; i < num_filters7; i++)
	{
		img_fltr_p6 = img_fltr_6;
		for (int j = 0; j < num_channels7; j++)
		{
			// reading corresponding weights to kernel
			for (int cnt_kernel = 0; cnt_kernel < filtersize7; cnt_kernel++)
			{
				*(kernel7 + cnt_kernel) = weights_layer7[cnt_weight + cnt_kernel];
			}

			imfilter(img_fltr_p6, kernel7, img_fltr_7_tmp, rows, cols, padsize7);
			imadd(img_fltr_p7, img_fltr_7_tmp, cols, rows);

			cnt_weight = cnt_weight + filtersize7;
			img_fltr_p6 = img_fltr_p6 + rows*cols;
		}
		double bias_tmp = biases_layer7[i];
		PReLU(img_fltr_p7, rows, cols, bias_tmp, prelu_coeff_layer7);
		img_fltr_p7 = img_fltr_p7 + rows*cols;
	}

	free(img_fltr_6);
	img_fltr_6 = NULL;
	free(img_fltr_6_tmp);
	img_fltr_6_tmp = NULL;
	free(kernel6);
	kernel6 = NULL;

	/////////// Layer8
	// Reading weights of 8th layer
	FILE *weights_layer8_ptr;
	weights_layer8_ptr = fopen("weights_layer8.txt", "r");
	if (weights_layer8_ptr == NULL) { printf("Error in the reading weights of 8th layer\n"); };
	double weights_layer8[4536];
	for (int i = 0; i < 4536; i++)
	{
		fscanf(weights_layer8_ptr, "%lf", &weights_layer8[i]);
	}
	fclose(weights_layer8_ptr);
	// Reading biases of 8th layer
	double biases_layer8 = - 0.03262640000;
	// Other parameters
	int filtersize8 = 81; //9x9
	int patchsize8 = 9;
	int num_filters8 = 1;
	int num_channels8 = 56;

	// Decvolution ==> output is outP_tmp
	double *img_fltr_8 = (double *)calloc((rows*scale) *(cols*scale) * num_filters8 , sizeof(double));
	double *kernel8 = (double *)malloc(filtersize8*sizeof(double));
	double *img_fltr_8_tmp = (double *)malloc((rows*scale) *(cols*scale) * sizeof(double));
	
	cnt_weight = 0;
	img_fltr_p7 = img_fltr_7;

	for (int j = 0; j < num_channels8; j++)
	{
		// reading corresponding weights to kernel
		for (int cnt_kernel = 0; cnt_kernel < filtersize8; cnt_kernel++)
		{
			*(kernel8 + cnt_kernel) = weights_layer8[cnt_weight + cnt_kernel];
		}
		cnt_weight = cnt_weight + filtersize8;

		deconv(img_fltr_p7, img_fltr_8_tmp, kernel8, cols, rows, scale);

		imadd(img_fltr_8, img_fltr_8_tmp, cols*scale, rows*scale);
		
		img_fltr_p7 = img_fltr_p7 + rows*cols;
	}


	for (int i=0;i<rows*scale;i++)
	for (int j = 0;j<cols*scale; j++)
	{
		int cnt_fnl = i*cols*scale + j;
		*(outP_tmp + cnt_fnl) = *(img_fltr_8 + cnt_fnl) + biases_layer8;
	}

	free(img_fltr_7);
	img_fltr_7 = NULL;
	free(img_fltr_7_tmp);
	img_fltr_7_tmp = NULL;
	free(kernel7);
	kernel7 = NULL;

	free(img_fltr_8);
	img_fltr_8 = NULL;
	free(img_fltr_8_tmp);
	img_fltr_8_tmp = NULL;
	free(kernel8);
	kernel8 = NULL;
	free(s_img_fltr_5);


	outP_tmp = outBuf_tmp;
		
	for (i = 0; i<inRows*scale; i++)
	for (j = 0; j<inCols*scale; j++)
	{
		int cnt = i*inCols*scale + j;
		*(outP_tmp + cnt) = *(outP_tmp + cnt) * 255;
	}
	double_2_uint8( outP_tmp, outP, outCols, outRows);

	fwrite(outBuf, sizeof(unsigned char), outCols*outRows, outFp);

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
		// V COmponent
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