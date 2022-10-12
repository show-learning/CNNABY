#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <iostream>
#include <string.h>
#include <fstream>
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
using namespace std;
//since hardware limit, we partition L5 for_loop to 12 times execute in makefile

ABYParty *party;
BooleanCircuit *bc;
ArithmeticCircuit *ac;

void buildparty(std::string address, std::string p, std::string r, std::string path)
{
	// default argument
	uint16_t port = (uint16_t)std::stoi(p);
	e_role role = (r == "SERVER") ? SERVER : CLIENT;
	uint32_t secparam = 128;
	seclvl seclvl = get_sec_lvl(secparam);
	e_mt_gen_alg mt_alg = MT_OT;
	// std::string circuit_dir = "/home/wei/ABY/bin/circ/";
	uint32_t nthreads = 1, bitlen = 64, out_n, out_l;
	// build party
	party = new ABYParty(role, address, port, seclvl, bitlen, nthreads, mt_alg, 100000, path);
	std::vector<Sharing *> &sharings = party->GetSharings();
	bc = (BooleanCircuit *)sharings[S_BOOL]->GetCircuitBuildRoutine();
	ac = (ArithmeticCircuit *)sharings[S_ARITH]->GetCircuitBuildRoutine();
}
void s_pad_image(SHR *img, SHR *img_pad, int rows, int cols, int padsize)
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
		img_pad[cnt_pad] = img[cnt];
	}
	// Top and Bottom Rows
	for (j = padsize; j < cols_pad - padsize; j++)
	for (k = 0; k < padsize; k++)
	{
		// Top Rows 
		cnt_pad = j + k*cols_pad;
		cnt = j - padsize;
		img_pad[cnt_pad] = img[cnt];
		// Bottom Rows
		cnt_pad = j + (rows_pad - 1 - k)* cols_pad;
		cnt = (j - padsize) + (rows - 1)*cols;
		img_pad[cnt_pad] = img[cnt];
	}
	// Left and Right Columns
	for (i = padsize; i < rows_pad - padsize; i++)
	for (k = 0; k < padsize; k++)
	{
		// Left Columns
		cnt = (i - padsize)*cols;
		cnt_pad = i*cols_pad + k;
		img_pad[cnt_pad] = img[cnt];
		// Right Columns
		cnt = (i - padsize)*cols + cols - 1;
		cnt_pad = i*cols_pad + cols_pad - 1 - k;
		img_pad[cnt_pad] = img[cnt];
	}
	// Corner Pixels
	for (k1 = 0; k1 < padsize; k1++)
	for (k2 = 0; k2 < padsize; k2++)
	{
		// Upper Left Corner
		cnt_pad = k1*cols_pad + k2;
		img_pad[cnt_pad] = img[0];
		// Upper Right Corner
		cnt_pad = k1*cols_pad + cols_pad - 1 - k2;
		img_pad[cnt_pad] = img[cols - 1];
		// Lower Left Corner
		cnt_pad = (rows_pad - 1 - k1)*cols_pad + k2;
		img_pad[cnt_pad] = img[(rows - 1)*cols];
		// Lower Right Corner
		cnt_pad = (rows_pad - 1 - k1)*cols_pad + cols_pad - 1 - k2;
		img_pad[cnt_pad] = img[(rows - 1)*cols + cols - 1];
	}
}

void s_imfilter(SHR *img, SHR *kernel, SHR *img_fltr, int rows, int cols, int padsize)
{
	int cols_pad = cols + 2 * padsize;//176+2
	int rows_pad = rows + 2 * padsize;//144+2
	int i, j, cnt, cnt_pad, cnt_krnl, k1, k2;
	double_t z=0.0;
	SHR sum( (uint64_t*)&z ,bc,party);
	SHR *img_pad=(SHR *)malloc(rows_pad * cols_pad*sizeof(SHR));
	s_pad_image(img, img_pad, rows, cols, padsize);//pass
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
		*(img_fltr+cnt) = sum;
	}
	free(img_pad);
}

SHR s_Max(SHR a, SHR b){
	SHR cond=a>b;
	SHR result=cond.mux(a,b);
	return result;
}

SHR s_Min(SHR a, SHR b){
	SHR cond=a<b;
	SHR result=cond.mux(a,b);
	return result;
}

void s_PReLU(SHR *img_fltr,int rows, int cols, SHR bias, double prelu_coeff)
{
	double_t z=0.0;
	SHR zero((uint64_t*)&z ,bc,party, (bool) 1);
	int cnt = 0;
	for (int i = 0; i < rows;i++)
	for (int j = 0; j < cols; j++)
	{
		cnt = i*cols + j;
		//*(img_fltr + cnt) = s_Max(*(img_fltr + cnt) + bias, zero) + prelu_coeff * s_Min(*(img_fltr + cnt) + bias, 0);
		*(img_fltr+cnt) = s_Max(*(img_fltr+cnt) + bias, zero) +  s_Min(*(img_fltr+cnt) + bias, zero) * prelu_coeff;
	}
}

void s_imadd(SHR *img_fltr_sum, SHR *img_fltr_crnt, int cols, int rows)
{

	int cnt = 0;
	for (int i = 0; i < rows;i++)
	for (int j = 0; j < cols; j++)
	{
		cnt = i*cols + j;
		//*(img_fltr_sum + cnt) = *(img_fltr_sum + cnt) + *(img_fltr_crnt + cnt);
		*(img_fltr_sum+cnt) = *(img_fltr_sum+cnt) + *(img_fltr_crnt+cnt);
	}
}

int main(int argc, char **argv)
{

	int cols = 176; // Width of input (downsampled) video
	int rows = 144; // Height of input (downsampled) video
	FILE *inFp, *outFp, *filter_cnt;
	int filter; //filter means times
    if((filter_cnt = fopen("filter_cnt.txt", "r"))!=NULL){
        fscanf(filter_cnt,"%d",&filter);
        fclose(filter_cnt);
		if(filter >=12){//do 12 times(0~11 OK) -> end
			return 0;
		}
    }
    else
        filter = 0;

	std::string address(argv[1]), port(argv[2]), r(argv[3]), path(argv[4]);
	e_role role = (r == "SERVER") ? SERVER : CLIENT;
	buildparty(address, port, r, path);
    
    inFp = fopen("L4_output.txt", "r");
    int num_filters4 = 12;
    double *img_fltr_4 = (double *)calloc(rows * cols * num_filters4 , sizeof(double));
    double *img_fltr_p4 = img_fltr_4;
	SHR *s_img_fltr_4 = (SHR *)malloc(rows * cols * num_filters4*sizeof(SHR));
	SHR *s_img_fltr_p4 = s_img_fltr_4;
	for (int i = 0; i < rows * cols * num_filters4; i++)
	{
		fscanf(inFp, "%lf", (img_fltr_4+i));
		//double stemp=(double)img_fltr_4[i];
		s_img_fltr_4[i] = SHR((uint64_t *)&img_fltr_4[i], bc, party);
	}
	//cout<<"read output4 done"<<endl;
    /////////// Layer5
	// Reading weights of 5th layer
	FILE *weights_layer5_ptr;
	weights_layer5_ptr = fopen("weights_layer5.txt", "r");
	if (weights_layer5_ptr == NULL) { printf("Error in the reading weights of 5th layer\n"); };
	double weights_layer5[1296];
	SHR s_weights_layer5[1296];
	for (int i = 0; i < 1296; i++)
	{
		uint64_t *shared;
		if(role==SERVER){
			fscanf(weights_layer5_ptr, "%lf", &weights_layer5[i]);
			shared=new uint64_t(*(uint64_t*)&weights_layer5[i]);
		}
		else{
			shared=NULL;
		}
		*(s_weights_layer5+i)= SHR(shared, bc, party, role);
		delete shared;
	}
	fclose(weights_layer5_ptr);
	// Reading biases of 5th layer
	FILE *biases_layer5_ptr;
	biases_layer5_ptr = fopen("biasess_layer5.txt", "r");
	if (biases_layer5_ptr == NULL) { printf("Error in the reading biases of 5th layer\n"); };
	double biases_layer5[12];
	SHR s_biases_layer5[12];
	for (int i = 0; i < 12; i++)
	{
		uint64_t *shared;
		if(role==SERVER){
			fscanf(biases_layer5_ptr, "%lf", &biases_layer5[i]);
			shared=new uint64_t(*(uint64_t*)&biases_layer5[i]);
		}
		else{
			shared=NULL;
		}
		*(s_biases_layer5+i) = SHR(shared, bc, party, role);
		delete shared;
	}
	fclose(biases_layer5_ptr);
	// Other parameters
	int filtersize5 = 9; //3X3
	int patchsize5 = 3;
	int padsize5 = (patchsize5 - 1) / 2;
	int num_filters5 = 12;
	int num_channels5 = 12;
	double prelu_coeff_layer5 = 0.3495;
    
    
	// Convolution
	double *img_fltr_5 = (double *)calloc(rows * cols * num_filters5 , sizeof(double));
	/*double *kernel5 = (double *)malloc(filtersize5*sizeof(double));
	double *img_fltr_p5 = img_fltr_5; // Pointer to img_fltr5
	double *img_fltr_5_tmp = (double *)malloc(rows * cols * sizeof(double));*/
	SHR *s_kernel5 = (SHR *)malloc(filtersize5 * sizeof(SHR));;
	SHR *s_img_fltr_5 = (SHR *)malloc(rows * cols * num_filters5*sizeof(SHR));
	SHR *s_img_fltr_p5 = s_img_fltr_5;
	SHR *s_img_fltr_5_tmp = (SHR *)malloc(rows * cols * sizeof(SHR));
    double zero = 0.0;
    FILE* fltr_5_tmp;
    if(filter==0){
        for(int i=0;i<rows * cols * num_filters5;i++){
		    uint64_t *shared;
		    if(role==SERVER){
                shared=new uint64_t(*(uint64_t*)&zero);
		    	
		    }
		    else{
		    	shared=NULL;
		    }
		    *(s_img_fltr_p5+i)= SHR(shared, bc, party, role);
		    delete shared;
	    }
    }
    else{
		char fname[20];
		sprintf(fname,"L5_output_%d.txt", filter);
        fltr_5_tmp = fopen(fname, "r");
        for (int i = 0; i < rows * cols * num_filters5; i++)
	    {
            fscanf(fltr_5_tmp, "%lf", (img_fltr_5+i));
			//double stemp=(double)img_fltr_4[i];
			s_img_fltr_5[i] = SHR((uint64_t *)&img_fltr_5[i], bc, party);
	    }
		fclose(fltr_5_tmp);
    }
	
	SHR s_bias_tmp( (uint64_t*)&zero ,bc,party);
	int cnt_weight = filter*filtersize5*num_channels5;
    //double bias_tmp;
	//cout<<"Layer5 START"<<endl;

    s_img_fltr_p4 = s_img_fltr_4;
    s_img_fltr_p5 = s_img_fltr_5 + rows*cols*filter;
	for (int j = 0; j < num_channels5; j++)
	{
		// reading corresponding weights to kernel
		for (int cnt_kernel = 0; cnt_kernel < filtersize5; cnt_kernel++)
		{
			//*(kernel5 + cnt_kernel) = weights_layer5[cnt_weight + cnt_kernel];
			*(s_kernel5+cnt_kernel)= *(s_weights_layer5 + cnt_weight + cnt_kernel);
		}
		//imfilter(img_fltr_p4, kernel5, img_fltr_5_tmp, rows, cols, padsize5);
		//imadd(img_fltr_p5, img_fltr_5_tmp, cols, rows);
		//cnt_weight = cnt_weight + filtersize5;
		//img_fltr_p4 = img_fltr_p4 + rows*cols;
		s_imfilter(s_img_fltr_p4, s_kernel5, s_img_fltr_5_tmp, rows, cols, padsize5);
		//cout <<"imfilter done"<<endl;	
		s_imadd(s_img_fltr_p5, s_img_fltr_5_tmp, cols, rows);
		//cout<< "imadd done" <<endl;
		cnt_weight = cnt_weight + filtersize5;
		s_img_fltr_p4 = s_img_fltr_p4 + rows*cols;
		//cout <<"filiters = "<<filter<<" channels = "<<j<<endl;
    }
	//bias_tmp = biases_layer5[i];
	//PReLU(img_fltr_p5, rows, cols, bias_tmp, prelu_coeff_layer5);
	//img_fltr_p5 = img_fltr_p5 + rows*cols;
	s_bias_tmp = *(s_biases_layer5+filter);
	s_PReLU(s_img_fltr_p5, rows, cols, s_bias_tmp, prelu_coeff_layer5);
	//cout<<"prelu done"<<endl;
	
	//cout<<"num_filters "<<filter<<" finished"<<endl;

    //outFp = fopen("L5_output.txt", "wb");
	for (int i = 0; i < rows * cols * num_filters5; i++)
	{
		//fprintf(outFp,"%lf ", s_img_fltr_5[i].get_d());
		cout <<s_img_fltr_5[i].get_d()<<endl;
	}
   	//fclose(outFp);
    filter++;
    FILE* filter_cnt_p = fopen("filter_cnt.txt", "wb");	
	fprintf(filter_cnt_p,"%d \n",filter);
   	fclose(filter_cnt_p);

	free(img_fltr_4);
	free(s_img_fltr_4);
	free(s_img_fltr_5);
	free(s_img_fltr_5_tmp);
	free(s_kernel5);
	return 0;
}