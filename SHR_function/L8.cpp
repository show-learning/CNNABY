#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <ENCRYPTO_utils/crypto/crypto.h>
#include <abycore/aby/abyparty.h>
#include <abycore/circuit/share.h>
#include <abycore/circuit/booleancircuits.h>
#include <abycore/circuit/arithmeticcircuits.h>
#include <abycore/sharing/sharing.h>
#include <iostream>
#include<bits/stdc++.h>
#include <aby.class.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <string>
#include <netinet/in.h>
// #include <R_ext/Print.h>
// #include <R.h>
// #include <Rinternals.h>
// #include <Rembedded.h>
// #include <Rdefines.h> 
#define inCols 176
#define inRows 144
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
// Replicate image padding by the factor of "padsize"
void s_pad_image(SHR *img, SHR*img_pad, int rows, int cols, int padsize)
//void pad_image(double *img, double *img_pad, int rows, int cols, int padsize)
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
		// double x = *(img + cnt);
		// *(img_pad + cnt_pad) = x;
		*(img_pad + cnt_pad) = *(img + cnt);
		//*(img_pad + cnt_pad) = *(img + cnt);
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
// void imfilter(double *img, double *kernel, double *img_fltr, int rows, int cols, int padsize)
// {
// 	// img_pad is the pointer to padded image
// 	// kernel is the pointer to the kernel which used for convolution
// 	// img_fltr is the pointer to the filtered image by applying convolution
// 	int cols_pad = cols + 2 * padsize;
// 	int rows_pad = rows + 2 * padsize;
// 	int i, j, cnt, cnt_pad, cnt_krnl, k1, k2;
// 	double sum;

// 	double *img_pad = (double *)malloc(rows_pad * cols_pad * sizeof(double));
	
// 	pad_image(img, img_pad, rows, cols, padsize);

// 	for (i = padsize; i < rows_pad - padsize; i++)
// 	for (j = padsize; j < cols_pad - padsize; j++)
// 	{
// 		cnt = (i - padsize)*cols + (j - padsize); // counter which shows current pixel in filtered image (central pixel in convolution window)
// 		sum = 0;
// 		cnt_krnl = 0; // counter which determines kernel elements
// 		for (k1 = -padsize; k1 <= padsize; k1++)
// 		for (k2 = -padsize; k2 <= padsize; k2++)
// 		{
// 			cnt_pad = (i + k1)*cols_pad + j + k2; // counter which shows each neighbouring pixel of padded image used for convolution with kernel
// 			sum = sum + (*(img_pad + cnt_pad))*(*(kernel + cnt_krnl));
// 			cnt_krnl++;
// 		}
// 		*(img_fltr + cnt) = sum;
// 	}

// 	free(img_pad);
// 	img_pad = NULL;
// }



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

void s_deconv(SHR *img_input, SHR *img_output, SHR *kernel, int cols, int rows, int stride)
{
	int border = 1;
	int fsize = 9;
	int rows_pad = rows + 2 * border;
	int cols_pad = cols + 2 * border;
	// double *img_input_padded = (double *)malloc(rows_pad * cols_pad * sizeof(double));
	// pad_image(img_input, img_input_padded, rows, cols, border);
	SHR *s_img_input_padded = (SHR *)malloc(rows_pad * cols_pad * sizeof(SHR));
	s_pad_image(img_input, s_img_input_padded, rows, cols, border);
	cout <<"spad complete" <<endl;
	int rows_out_pad = rows_pad * stride;
	int cols_out_pad = cols_pad * stride;
	//double *img_output_tmp = (double *)calloc((rows_out_pad + fsize - 1)* (cols_out_pad + fsize - 1), sizeof(double));
	//double *kernel_modif = (double *)malloc(fsize * fsize * sizeof(double));
	SHR *s_kernel_modif = (SHR *)malloc(fsize * fsize * sizeof(SHR));
	SHR *s_img_output_tmp = (SHR *)malloc((rows_out_pad + fsize - 1)* (cols_out_pad + fsize - 1)* sizeof(SHR));
	double_t z=0.0;
	SHR zero( (uint64_t*)&z ,bc,party);
	for (int i = 0; i < (rows_out_pad + fsize - 1)* (cols_out_pad + fsize - 1); i++){
		*(s_img_output_tmp + i) = zero;
	}
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
			//*(kernel_modif + cnt_kernel) = (*(kernel + cnt_kernel))*(*(img_input_padded + cnt_img));
			//*(img_output_tmp + cnt_img_output + k_c) = *(img_output_tmp + cnt_img_output + k_c) + *(kernel_modif + cnt_kernel);
			//*(s_kernel_modif + cnt_kernel) = (*(kernel + cnt_kernel))*(*(s_img_input_padded + cnt_img));
			*(s_img_output_tmp + cnt_img_output + k_c) = ((*(s_img_output_tmp + cnt_img_output + k_c)) + (*(kernel + cnt_kernel))*(*(s_img_input_padded + cnt_img)));
		}
		cnt_img_output = cnt_img_output + (cols_out_pad + fsize - 1);
	    }
		cout << i << " " << j <<endl;
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
		//*(img_output + cnt_img_out) = *(img_output_tmp + cnt_img_out_tmp);
		*(img_output + cnt_img_out) = *(s_img_output_tmp + cnt_img_out_tmp);
	}

	//free(img_input_padded); img_input_padded = NULL;
	// free(img_output_tmp); img_output_tmp = NULL;
	// free(kernel_modif); kernel_modif = NULL;
	free(s_img_input_padded); s_img_input_padded = NULL;
	free(s_img_output_tmp); s_img_output_tmp = NULL;
	free(s_kernel_modif); s_kernel_modif = NULL;
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

int main(int argc, char *argv[]){
    std::string address(argv[1]),port(argv[2]),r(argv[3]),path(argv[4]);
	e_role role=(r=="SERVER")?SERVER:CLIENT;
	buildparty(address,port,r,path);
    int cols = 176; //Width of input (downsampled) video
	int rows = 144; //Height of input (downsampled) video
	int scale =2;
	int outCols = cols * scale;
	int outRows = rows * scale;
    FILE *inFp, *outFp,*L7Fp;
	if(role==SERVER){
		L7Fp = fopen("L7_output.txt", "r");
	}
	int num_filters7 = 56;
    double *img_fltr_7 = (double *)calloc(rows * cols * num_filters7 , sizeof(double));
    double *img_fltr_p7 = img_fltr_7;
	//rows * cols * num_filters7 = 1419264
	SHR *s_img_fltr_7 = (SHR*)calloc(rows * cols * num_filters7 , sizeof(SHR));
	SHR *s_img_fltr_p7 = s_img_fltr_7;
	//std::array<SHR, 1419264> s_img_fltr_7;
    for (int i = 0; i < rows * cols * num_filters7; i++)
	{
		//SHR assign
		uint64_t *shared;
		if(role==SERVER){
			fscanf(L7Fp, "%lf", (img_fltr_7+i));
			shared = new uint64_t(*(uint64_t *)&img_fltr_7[i]);
		}
		else{
			shared = NULL;
		}
		*(s_img_fltr_7 + i) = SHR(shared, bc, party, role);
	}
	if(role==SERVER){
		fclose(L7Fp);
	}
    /////////// Convolution3 ------------------- Layer 8
	cout << "if success" << endl;
	/////////// Layer8
	// Reading weights of 8th layer
	FILE *weights_layer8_ptr;
	if(role==SERVER){
		weights_layer8_ptr = fopen("weights_layer8.txt", "r");
	}
	//if (weights_layer8_ptr == NULL) { printf("Error in the reading weights of 8th layer\n"); };
	double weights_layer8[4536];

	SHR* s_weights_layer8= (SHR*)calloc(4536 , sizeof(SHR));
	//std::array<SHR, 4536> s_weights_layer8;
	for (int i = 0; i < 4536; i++)
	{
		uint64_t *shared;
		if(role==SERVER){
			fscanf(weights_layer8_ptr, "%lf", &weights_layer8[i]);
			shared = new uint64_t(*(uint64_t *)&weights_layer8[i]);
		}
		else{
			shared = NULL;
		}
		*(s_weights_layer8+i) = SHR(shared, bc, party, role);
	}
	if(role==SERVER){
		fclose(weights_layer8_ptr);
	}
	
	// Reading biases of 8th layer
	double biases_layer8 = - 0.03262640000;
	// Other parameters
	int filtersize8 = 81; //9x9
	int patchsize8 = 9;
	int num_filters8 = 1;
	int num_channels8 = 56;

	// Decvolution ==> output is the img_hr
	double *img_fltr_8 = (double *)calloc((rows*scale) *(cols*scale) * num_filters8 , sizeof(double));
	double *kernel8 = (double *)malloc(filtersize8*sizeof(double));
	double *img_fltr_8_tmp = (double *)malloc((rows*scale) *(cols*scale) * sizeof(double));
	
	SHR *s_img_fltr_8 = (SHR*)calloc(101376 , sizeof(SHR));
	SHR *s_kernel8 = (SHR*)malloc(filtersize8*sizeof(SHR));
	SHR *s_img_fltr_8_tmp = (SHR*)calloc(101376 , sizeof(SHR));
	// std::array<SHR, 101376> s_img_fltr_8;
	// std::array<SHR, 81> s_kernel8;
	// std::array<SHR, 101376> s_img_fltr_8_tmp;
	int cnt_weight = 0;
	img_fltr_p7 = img_fltr_7;
	s_img_fltr_p7 = s_img_fltr_7;
	for (int j = 0; j < num_channels8; j++)
	{
		// reading corresponding weights to kernel
		for (int cnt_kernel = 0; cnt_kernel < filtersize8; cnt_kernel++)
		{
			//origin
			*(kernel8 + cnt_kernel) = weights_layer8[cnt_weight + cnt_kernel];
			//SHR
			SHR x = *(s_weights_layer8+cnt_weight + cnt_kernel);
			*(s_kernel8+cnt_kernel) = x;
		}
		cnt_weight = cnt_weight + filtersize8;

		deconv(img_fltr_p7, img_fltr_8_tmp, kernel8, cols, rows, scale);
		
		s_deconv(s_img_fltr_p7, s_img_fltr_8_tmp, s_kernel8, cols, rows, scale);
		cout << "s_ok" << endl;
		cout << img_fltr_8_tmp[1] << s_img_fltr_8_tmp[1].check() << endl;

		imadd(img_fltr_8, img_fltr_8_tmp, cols*scale, rows*scale);
		
		img_fltr_p7 = img_fltr_p7 + rows*cols;
	}
	double *outBuf_tmp = (double *)malloc(outCols*outRows*sizeof(double));
	double *outP_tmp = outBuf_tmp;
	for (int i=0;i<rows*scale;i++)
	for (int j = 0;j<cols*scale; j++)
	{
		int cnt_fnl = i*cols*scale + j;
		*(outP_tmp + cnt_fnl) = *(img_fltr_8 + cnt_fnl) + biases_layer8;
	}
	/*for ( int i = 0; i < 10; i++)
	{
		printf("%f\n", *(outP_tmp + i));
	}*/
	free(img_fltr_8);
	img_fltr_8 = NULL;
	free(img_fltr_8_tmp);
	img_fltr_8_tmp = NULL;
	free(kernel8);
	kernel8 = NULL;
	unsigned char *inBuf = (unsigned char *)malloc(cols*rows*sizeof(unsigned char));
	unsigned char *inP = inBuf;
	unsigned char *outBuf = (unsigned char *)malloc(outCols*outRows*sizeof(unsigned char));
	unsigned char *outP = outBuf;

	outP_tmp = outBuf_tmp;
		
		for (int i = 0; i<rows*scale; i++)
		for (int j = 0; j<cols*scale; j++)
		{
			int cnt = i*cols*scale + j;
			*(outP_tmp + cnt) = *(outP_tmp + cnt) * 255;
		}
	double_2_uint8( outP_tmp, outP, outCols, outRows);
	int i,j;
	inFp = fopen("akiyo_qcif.yuv", "rb");
	outFp = fopen("test_cif.yuv", "wb");
	fwrite(outBuf, sizeof(unsigned char), outCols*outRows, outFp);
	//U Component
	fread(inBuf, sizeof(unsigned char), cols*rows, inFp);
	fread(inBuf, sizeof(unsigned char), cols*rows / 4, inFp);
	inP = inBuf;
	outP = outBuf;
	for (i = 0; i < rows / 2; i++)
	for (j = 0; j < cols / 2; j++) {
		int cnt = 2 * (i * outCols / 2 + j);
		unsigned char x = *inP++;
		*(outP + cnt) = x;
		*(outP + cnt + 1) = x;
		*(outP + cnt + outCols / 2) = x;
		*(outP + cnt + outCols / 2 + 1) = x;
	}
	fwrite(outBuf, sizeof(unsigned char), outCols*outRows / 4, outFp);
	// V Component
	fread(inBuf, sizeof(unsigned char), cols*rows / 4, inFp);
	inP = inBuf;
	outP= outBuf;
	for (i = 0; i < rows / 2; i++)
	for (j = 0; j < cols / 2; j++) {
		int cnt = 2 * (i*outCols / 2 + j);
		unsigned char x = *inP++;
		*(outP + cnt) = x;
		*(outP + cnt + 1) = x;
		*(outP + cnt + outCols / 2) = x;
		*(outP + cnt + outCols / 2 + 1) = x;
	}
	fwrite(outBuf, sizeof(unsigned char), outCols*outRows / 4, outFp);
	
	free(inBuf);
	inBuf = NULL;

	free(outBuf);
	outBuf = NULL;
	free(outBuf_tmp);
	outBuf_tmp = NULL;
}