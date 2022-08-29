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


/////////// Layer4
int main(int argc, char *argv[]){
	std::string address(argv[1]),port(argv[2]),r(argv[3]),path(argv[4]);
	buildparty(address,port,r,path);
	e_role role=(r=="SERVER")?SERVER:CLIENT;
	// Reading weights of 4th layer
	FILE *weights_layer4_ptr;
	weights_layer4_ptr = fopen("weights_layer4.txt", "r");
	if (weights_layer4_ptr == NULL) { printf("Error in the reading weights of 4th layer\n"); };
	double weights_layer4[1296];
	for (int i = 0; i < 1296; i++)
	{
		fscanf(weights_layer4_ptr, "%lf", &weights_layer4[i]);
	}
	fclose(weights_layer4_ptr);
	// Reading biases of 4th layer
	FILE *biases_layer4_ptr;
	biases_layer4_ptr = fopen("biasess_layer4.txt", "r");
	if (biases_layer4_ptr == NULL) { printf("Error in the reading biases of 4th layer\n"); };
	double biases_layer4[12];
	for (int i = 0; i < 12; i++)
	{
		fscanf(biases_layer4_ptr, "%lf", &biases_layer4[i]);
	}
	fclose(biases_layer4_ptr);
	// Other parameters
	int filtersize4 = 9; //3X3
	int patchsize4 = 3;
	int padsize4 = (patchsize4 - 1) / 2;
	int num_filters4 = 12;
	int num_channels4 = 12;
	double prelu_coeff_layer4 = 0.2476;
	int cols = 176; 
	int rows = 144;
	// Convolution
	double *img_fltr_4 = (double *)calloc(rows * cols * num_filters4 , sizeof(double));
	double *kernel4 = (double *)malloc(filtersize4*sizeof(double));
	double *img_fltr_p4 = img_fltr_4; // Pointer to img_fltr4
	double *img_fltr_4_tmp = (double *)malloc(rows * cols * sizeof(double));

	int cnt_weight = 0;
	double bias_tmp;
    /////////////////////////////////////////////////
    int num_filters3 = 12;
    double *img_fltr_3 = (double *)calloc(rows * cols * num_filters3 , sizeof(double));
	double *img_fltr_p3 = img_fltr_3;
    FILE *layer3_ptr;
	layer3_ptr = fopen("L3_output.txt", "r");
    double tmp;
    int cnt=0;
	for (int i = 0; i < rows * cols * num_filters3; i++)
	{
		fscanf(layer3_ptr, "%lf", (img_fltr_3+i));
	}
    fclose(layer3_ptr);

	for (int i = 0; i < num_filters4; i++)
	{
		img_fltr_p3 = img_fltr_3; 
		for (int j = 0; j < num_channels4; j++)
		{
			// reading corresponding weights to kernel
			for (int cnt_kernel = 0; cnt_kernel < filtersize4; cnt_kernel++)
			{
				*(kernel4 + cnt_kernel) = weights_layer4[cnt_weight + cnt_kernel];
			}

			imfilter(img_fltr_p3, kernel4, img_fltr_4_tmp, rows, cols, padsize4);
			imadd(img_fltr_p4, img_fltr_4_tmp, cols, rows);

			cnt_weight = cnt_weight + filtersize4;
			img_fltr_p3 = img_fltr_p3 + rows*cols;
		}
		double bias_tmp = biases_layer4[i];
		PReLU(img_fltr_p4, rows, cols, bias_tmp, prelu_coeff_layer4);
		img_fltr_p4 = img_fltr_p4 + rows*cols;
	}

	SHR *s_img_fltr_4 = (SHR *)malloc(rows * cols * num_filters4*sizeof(SHR));
	// FILE *fp = fopen("L4_output.txt", "wb");
	// for(int i=0;i<rows * cols * num_filters4;i++){
	// 	uint64_t *shared;
	// 	if(role==SERVER){
	// 		shared=new uint64_t(*(uint64_t*)&img_fltr_4[i]);
	// 	}
	// 	else{
	// 		shared=NULL;
	// 	}
	// 	*(s_img_fltr_4+i)= SHR(shared, bc, party, role);
	// 	delete shared;
	// 	fprintf(fp,"%lf\n", s_img_fltr_4[i].get_d());
	// }
    // fclose(fp);
	for(int i=0;i<rows * cols * num_filters4;i++){
		uint64_t *shared;
		if(role==SERVER){
			shared=new uint64_t(*(uint64_t*)&img_fltr_4[i]);
		}
		else{
			shared=NULL;
		}
		*(s_img_fltr_4+i)= SHR(shared, bc, party, role);
		delete shared;
		cout<<"getd: "<< s_img_fltr_4[i].get_d()<<endl;;
	}
	free(s_img_fltr_4);
	// for (int i = 0; i < 56; i++)
    // {
    //     uint64_t *shared;
    //     if (role == SERVER)
    //     {
    //         shared = new uint64_t(*(uint64_t *)&biases_layer7[i]);
    //     }
    //     else
    //     {
    //         shared = NULL;
    //     }
    //     s_biases_layer7[i] = SHR(shared, bc, party, role);
    // }

	free(img_fltr_3);
	img_fltr_3 = NULL;
	// free(img_fltr_3_tmp);
	// img_fltr_3_tmp = NULL;
	// free(kernel3);
	// kernel3 = NULL;
	free(img_fltr_4_tmp);
	img_fltr_4_tmp = NULL;
}


