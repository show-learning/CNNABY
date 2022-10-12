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
	//cout <<"spad complete" <<endl;
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
		//cout << i << " " << j <<endl;
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


	}
}

int main(int argc, char *argv[])
{   
    std::string address(argv[1]), port(argv[2]), r(argv[3]), path(argv[4]);
	e_role role = (r == "SERVER") ? SERVER : CLIENT;
	buildparty(address, port, r, path);
	char *inFile = argv[5];
	char *outFile = argv[6];

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
	//double *inBuf_tmp = (double *)malloc(inCols*inRows*sizeof(double));
    SHR *s_inBuf_tmp = (SHR *)calloc(inCols * inRows, sizeof(SHR));
	double *outBuf_tmp = (double *)malloc(outCols*outRows*sizeof(double));

	for (int fcnt = 0; fcnt < num; fcnt++)
	{
		//////// Interpolate each frame using FSRCNN for Y component and simple repitition for U and V components
		// Pointer to obtain value of each tpixel of input frame
		unsigned char *inP = inBuf;
		//double *inP_tmp = inBuf_tmp;
        SHR *s_inP_tmp = s_inBuf_tmp;
		// Pointer to obtain value of each pixel of output frame
		unsigned char *outP = outBuf;
		double *outP_tmp = outBuf_tmp;

		//Y Component
		fread(inBuf, sizeof(unsigned char), inCols*inRows, inFp);
		int i, j;

		for (i = 0; i<inRows; i++){
            for (j = 0; j<inCols; j++)
            {
                int cnt = i*inCols + j;
                int x = *inP++;
                *(s_inP_tmp + cnt) = (double)(x / 255.0);

                // SHR assign
                uint64_t *shared;
                if (role == CLIENT){
                    shared = new uint64_t(*(uint64_t *)&s_inP_tmp[cnt]);
                }
                else{
                    shared = NULL;
                }
                s_inP_tmp[cnt] = SHR(shared, bc, party, role);
                delete shared;
            }
        }

//====== FSRCNN start
		//FSRCNN(outP_tmp, inP_tmp, inRows, inCols, scale);
        // General Settings
        int num_layers = 8;
        /////////// Convolution1 -------- Layer1
        // Reading weights of first layer
        FILE *weights_layer1_ptr;
        weights_layer1_ptr = fopen("weights_layer1.txt", "r");
        if (weights_layer1_ptr == NULL){
            printf("Error in the reading weights of first layer\n");
        };
        double weights_layer1[1400];
		SHR s_weights_layer1[1400];

        for (int i = 0; i < 1400; i++)
        {
            uint64_t *shared;
			if (role == SERVER){
				fscanf(weights_layer1_ptr, "%lf", &weights_layer1[i]);
				shared = new uint64_t(*(uint64_t *)&weights_layer1[i]);
			}
			else{
				shared = NULL;
			}
			s_weights_layer1[i] = SHR(shared, bc, party, role);
			delete shared;
        }
        fclose(weights_layer1_ptr);

        // Reading biases of first layer
        FILE *biases_layer1_ptr;
        biases_layer1_ptr = fopen("biasess_layer1.txt", "r");
        if (biases_layer1_ptr == NULL){
            printf("Error in the reading biases of first layer\n");
        };
        double biases_layer1[56];
        SHR s_biases_layer1[56];
        for (int i = 0; i < 56; i++){
            uint64_t *shared;
			if (role == SERVER)
			{
				fscanf(biases_layer1_ptr, "%lf", &biases_layer1[i]);
				shared = new uint64_t(*(uint64_t *)&biases_layer1[i]);
			}
			else
			{
				shared = NULL;
			}
			s_biases_layer1[i] = SHR(shared, bc, party, role);
			delete shared;
        }
        fclose(biases_layer1_ptr);

        // other parameters
        int filtersize = 25; //5X5
        int patchsize = 5;
        int padsize = (patchsize - 1) / 2;
        int num_filters = 56;
        double prelu_coeff_layer1 = -0.8986;

        // Convolution
        //double *img_fltr_1 = (double *)malloc(rows * cols * num_filters * sizeof(double));
        //double *kernel = (double *)malloc(filtersize*sizeof(double));
        //double *img_fltr_p1 = img_fltr_1; // Pointer to img_fltr1 ==>> Using this way to be able to shift it to access data

        SHR *s_kernel = (SHR *)malloc(filtersize * sizeof(SHR));
		SHR *s_img_fltr_1 = (SHR *)calloc(inRows * inCols * num_filters, sizeof(SHR));
		SHR *s_img_fltr_p1 = s_img_fltr_1;

        double zero = 0.0;
		for (int i = 0; i < inRows * inCols * num_filters; i++){
			uint64_t *shared;
			if (role == SERVER)
			{
				shared = new uint64_t(*(uint64_t *)&zero);
			}
			else
			{
				shared = NULL;
			}
			*(s_img_fltr_p1 + i) = SHR(shared, bc, party, role);
			delete shared;
		}

        int cnt_weight = 0;
        SHR s_bias_tmp((uint64_t *)&zero, bc, party);
        //double bias_tmp;
        for (int i = 0; i < num_filters; i++)
        {
            // reading corresponding weights to kernel pointer
            for (int cnt_kernel = 0; cnt_kernel < filtersize; cnt_kernel++)
            {
                *(s_kernel + cnt_kernel) = s_weights_layer1[cnt_weight + cnt_kernel];
            }
            cnt_weight = cnt_weight + filtersize;
            s_imfilter(s_inP_tmp, s_kernel, s_img_fltr_p1, inRows, inCols, padsize);

            s_bias_tmp = s_biases_layer1[i];
			s_PReLU(s_img_fltr_p1, inRows, inCols, s_bias_tmp, prelu_coeff_layer1);

            s_img_fltr_p1 = s_img_fltr_p1 + inCols * inRows;
        }

        /////////// Convolution2 ------------------- Layer 2~7

        /////////// Layer2
        // Reading weights of 2nd layer
        FILE *weights_layer2_ptr;
        weights_layer2_ptr = fopen("weights_layer2.txt", "r");
        if (weights_layer2_ptr == NULL) { printf("Error in the reading weights of 2nd layer\n"); };
        // Note: weights must be saved in a way which that corresponding weights of each channel can be read by pointer concept ==>> for this layer 12X56 matrix is reshaped to (12X56)*1 vector
        double weights_layer2[672];
        SHR s_weights_layer2[672];
        for (int i = 0; i < 672; i++)
        {
            uint64_t *shared;
            if (role == SERVER)
            {
                fscanf(weights_layer2_ptr, "%lf", &weights_layer2[i]);
                shared = new uint64_t(*(uint64_t *)&weights_layer2[i]);
            }
            else
            {
                shared = NULL;
            }
            s_weights_layer2[i] = SHR(shared, bc, party, role);
            delete shared;
        }
        fclose(weights_layer2_ptr);

        // Reading biases of 2nd layer
        FILE *biases_layer2_ptr;
        biases_layer2_ptr = fopen("biasess_layer2.txt", "r");
        if (biases_layer2_ptr == NULL) { printf("Error in the reading biases of 2nd layer\n"); };
        double biases_layer2[12];
        SHR s_biases_layer2[12];
        for (int i = 0; i < 12; i++)
        {
            uint64_t *shared;
            if (role == SERVER)
            {
                fscanf(biases_layer2_ptr, "%lf", &biases_layer2[i]);
                shared = new uint64_t(*(uint64_t *)&biases_layer2[i]);
            }
            else
            {
                shared = NULL;
            }
            s_biases_layer2[i] = SHR(shared, bc, party, role);
            delete shared;
        }
        fclose(biases_layer2_ptr);

        // Other parameters
        int filtersize2 = 1; //1X1
        int patchsize2 = 1;
        int padsize2 = (patchsize2 - 1) / 2;
        int num_filters2 = 12;
        int num_channels2 = 56;
        double prelu_coeff_layer2 = 0.3236;
        // Convolution
        //double *img_fltr_2 = (double *)calloc(rows * cols * num_filters2 , sizeof(double)); // use calloc to initialize all variables to zero
        //double *img_fltr_2_tmp = (double *)malloc(rows * cols * sizeof(double));
        //double *kernel2 = (double *)malloc(filtersize2*sizeof(double));
        //double *img_fltr_p2 = img_fltr_2; // Pointer to img_fltr2
        
        SHR *s_img_fltr_2 = (SHR *)calloc(inRows * inCols * num_filters2, sizeof(SHR));
        SHR *s_kernel2 = (SHR *)malloc(filtersize2 * sizeof(SHR));
        SHR *s_img_fltr_2_tmp = (SHR *)malloc(inRows * inCols * sizeof(SHR));
        SHR *s_img_fltr_p2 = s_img_fltr_2;
        zero = 0.0;
        for (int i = 0; i < inRows * inCols * num_filters2; i++)
        {
            uint64_t *shared;
            if (role == SERVER)
            {
                shared = new uint64_t(*(uint64_t *)&zero);
            }
            else
            {
                shared = NULL;
            }
            *(s_img_fltr_p2 + i) = SHR(shared, bc, party, role);
            delete shared;
        }

        cnt_weight = 0;
        //SHR s_bias_tmp((uint64_t *)&zero, bc, party);

        for (int i = 0; i < num_filters2; i++)
        {
            s_img_fltr_p1 = s_img_fltr_1; // Return pointer to the first of array which contains feature map of previous layer
            for (int j = 0; j < num_channels2; j++){
            // reading corresponding weights to kernel
                for (int cnt_kernel = 0; cnt_kernel < filtersize2; cnt_kernel++){
                    *(s_kernel2 + cnt_kernel) = s_weights_layer2[cnt_weight + cnt_kernel];
                }
                s_imfilter(s_img_fltr_p1, s_kernel2, s_img_fltr_2_tmp, inRows, inCols, padsize2);
                s_imadd(s_img_fltr_p2, s_img_fltr_2_tmp, inCols, inRows);

                cnt_weight = cnt_weight + filtersize2;
                s_img_fltr_p1 = s_img_fltr_p1 + inRows * inCols;
            }
            // double bias_tmp;
            s_bias_tmp = s_biases_layer2[i];

            // SHR
            // PReLU(s_img_fltr_p2, rows, cols, bias_tmp, prelu_coeff_layer2);
            s_PReLU(s_img_fltr_p2, inRows, inCols, s_bias_tmp, prelu_coeff_layer2);
            s_img_fltr_p2 = s_img_fltr_p2 + inRows * inCols;
        }

        free(s_img_fltr_1);
        free(s_kernel);
        free(s_img_fltr_2_tmp);
        
        
        /////////// Layer3
        // Reading weights of 3rd layer
        FILE *weights_layer3_ptr;
        weights_layer3_ptr = fopen("weights_layer3.txt", "r");
        if (weights_layer3_ptr == NULL) { printf("Error in the reading weights of 3rd layer\n"); };
        double weights_layer3[1296];
        SHR s_weights_layer3[1296];
        for (int i = 0; i < 1296; i++)
        {
            //fscanf(weights_layer3_ptr, "%lf", &weights_layer3[i]);
            uint64_t *shared;
            if (role == SERVER)
            {
                fscanf(weights_layer3_ptr, "%lf", &weights_layer3[i]);
                shared = new uint64_t(*(uint64_t *)&weights_layer3[i]);
            }
            else
            {
                shared = NULL;
            }
            s_weights_layer3[i] = SHR(shared, bc, party, role);
            delete shared;
        }
        fclose(weights_layer3_ptr);

        // Reading biases of 3rd layer
        FILE *biases_layer3_ptr;
        biases_layer3_ptr = fopen("biasess_layer3.txt", "r");
        if (biases_layer3_ptr == NULL) { printf("Error in the reading biases of 3rd layer\n"); };
        double biases_layer3[12];
        SHR s_biases_layer3[12];
        for (int i = 0; i < 12; i++)
        {
            //fscanf(biases_layer3_ptr, "%lf", &biases_layer3[i]);
            uint64_t *shared;
            if (role == SERVER)
            {
                fscanf(biases_layer3_ptr, "%lf", &biases_layer3[i]);
                shared = new uint64_t(*(uint64_t *)&biases_layer3[i]);
            }
            else
            {
                shared = NULL;
            }
            s_biases_layer3[i] = SHR(shared, bc, party, role);
            delete shared;
        }
        fclose(biases_layer3_ptr);
        // Other parameters
        int filtersize3 = 9; //3X3
        int patchsize3 = 3;
        int padsize3 = (patchsize3 - 1) / 2;
        int num_filters3 = 12;
        int num_channels3 = 12;
        double prelu_coeff_layer3 = 0.2288;
        // Convolution
        //double *img_fltr_3 = (double *)calloc(rows * cols * num_filters3 , sizeof(double));
        //double *kernel3 = (double *)malloc(filtersize3*sizeof(double));
        //double *img_fltr_p3 = img_fltr_3; // Pointer to img_fltr2
        //double *img_fltr_3_tmp = (double *)malloc(rows * cols * sizeof(double));
        SHR *s_img_fltr_3 = (SHR *)calloc(inRows * inCols * num_filters3, sizeof(SHR));
        SHR *s_kernel3 = (SHR *)malloc(filtersize3 * sizeof(SHR));
        SHR *s_img_fltr_p3 = s_img_fltr_3;
        SHR *s_img_fltr_3_tmp = (SHR *)malloc(inRows * inCols * sizeof(SHR));
        
        zero = 0.0;
        for (int i = 0; i < inRows * inCols * num_filters3; i++)
        {
            uint64_t *shared;
            if (role == SERVER)
            {
                shared = new uint64_t(*(uint64_t *)&zero);
            }
            else
            {
                shared = NULL;
            }
            *(s_img_fltr_p3 + i) = SHR(shared, bc, party, role);
            delete shared;
        }

        cnt_weight = 0;
        
        for (int i = 0; i < num_filters3; i++)
        {
            s_img_fltr_p2 = s_img_fltr_2; // Return pointer to the first cell of array which contains feature maps of previous layer
            for (int j = 0; j < num_channels3; j++)
            {
                // reading corresponding weights to kernel
                for (int cnt_kernel = 0; cnt_kernel < filtersize3; cnt_kernel++)
                {
                    *(s_kernel3 + cnt_kernel) = s_weights_layer3[cnt_weight + cnt_kernel];
                }

                s_imfilter(s_img_fltr_p2, s_kernel3, s_img_fltr_3_tmp, inRows, inCols, padsize3);
                s_imadd(s_img_fltr_p3, s_img_fltr_3_tmp, inCols, inRows);

                cnt_weight = cnt_weight + filtersize3;
                s_img_fltr_p2 = s_img_fltr_p2 + inRows*inCols;
            }
            s_bias_tmp = s_biases_layer3[i];
            s_PReLU(s_img_fltr_p3, inRows, inCols, s_bias_tmp, prelu_coeff_layer3);
            s_img_fltr_p3 = s_img_fltr_p3 + inRows*inCols;
        }

        free(s_img_fltr_2);
        free(s_img_fltr_2_tmp);
        free(s_kernel2);
        free(s_img_fltr_3_tmp);

        /////////// Layer4
        // Reading weights of 4th layer
        FILE *weights_layer4_ptr;
        weights_layer4_ptr = fopen("weights_layer4.txt", "r");
        if (weights_layer4_ptr == NULL) { printf("Error in the reading weights of 4th layer\n"); };
        double weights_layer4[1296];
        SHR s_weights_layer4[1296];
        for (int i = 0; i < 1296; i++)
        {
            //fscanf(weights_layer4_ptr, "%lf", &weights_layer4[i]);
            uint64_t *shared;
            if (role == SERVER){
                fscanf(weights_layer4_ptr, "%lf", &weights_layer4[i]);
                shared = new uint64_t(*(uint64_t *)&weights_layer4[i]);
            }
            else{
                shared = NULL;
            }
            s_weights_layer4[i] = SHR(shared, bc, party, role);
            delete shared;
        }
        fclose(weights_layer4_ptr);

        // Reading biases of 4th layer
        FILE *biases_layer4_ptr;
        biases_layer4_ptr = fopen("biasess_layer4.txt", "r");
        if (biases_layer4_ptr == NULL) { printf("Error in the reading biases of 4th layer\n"); };
        double biases_layer4[12];
        SHR s_biases_layer4[12];
        for (int i = 0; i < 12; i++)
        {
            //fscanf(biases_layer4_ptr, "%lf", &biases_layer4[i]);
            uint64_t *shared;
            if (role == SERVER){
                fscanf(biases_layer4_ptr, "%lf", &biases_layer4[i]);
                shared = new uint64_t(*(uint64_t *)&biases_layer4[i]);
            }
            else{
                shared = NULL;
            }
            s_biases_layer4[i] = SHR(shared, bc, party, role);
            delete shared;
        }
        fclose(biases_layer4_ptr);
        // Other parameters
        int filtersize4 = 9; //3X3
        int patchsize4 = 3;
        int padsize4 = (patchsize4 - 1) / 2;
        int num_filters4 = 12;
        int num_channels4 = 12;
        double prelu_coeff_layer4 = 0.2476;
        // Convolution
        // double *img_fltr_4 = (double *)calloc(rows * cols * num_filters4 , sizeof(double));
        // double *kernel4 = (double *)malloc(filtersize4*sizeof(double));
        // double *img_fltr_p4 = img_fltr_4; // Pointer to img_fltr4
        // double *img_fltr_4_tmp = (double *)malloc(rows * cols * sizeof(double));
        SHR *s_img_fltr_4 = (SHR *)calloc(inRows * inCols * num_filters4, sizeof(SHR));
        SHR *s_kernel4 = (SHR *)malloc(filtersize4 * sizeof(SHR));
        SHR *s_img_fltr_p4 = s_img_fltr_4;
        SHR *s_img_fltr_4_tmp = (SHR *)malloc(inRows * inCols * sizeof(SHR));

        zero = 0.0;
        for (int i = 0; i < inRows * inCols * num_filters4; i++)
        {
            uint64_t *shared;
            if (role == SERVER)
            {
                shared = new uint64_t(*(uint64_t *)&zero);
            }
            else
            {
                shared = NULL;
            }
            *(s_img_fltr_p4 + i) = SHR(shared, bc, party, role);
            delete shared;
        }

        cnt_weight = 0;

        for (int i = 0; i < num_filters4; i++)
        {
            s_img_fltr_p3 = s_img_fltr_3; 
            for (int j = 0; j < num_channels4; j++)
            {
                // reading corresponding weights to kernel
                for (int cnt_kernel = 0; cnt_kernel < filtersize4; cnt_kernel++)
                {
                    *(s_kernel4 + cnt_kernel) = s_weights_layer4[cnt_weight + cnt_kernel];
                }

                s_imfilter(s_img_fltr_p3, s_kernel4, s_img_fltr_4_tmp, inRows, inCols, padsize4);
                s_imadd(s_img_fltr_p4, s_img_fltr_4_tmp, inCols, inRows);

                cnt_weight = cnt_weight + filtersize4;
                s_img_fltr_p3 = s_img_fltr_p3 + inRows*inCols;
            }
            s_bias_tmp = s_biases_layer4[i];
            s_PReLU(s_img_fltr_p4, inRows, inCols, s_bias_tmp, prelu_coeff_layer4);
            s_img_fltr_p4 = s_img_fltr_p4 + inRows*inCols;
        }

        free(s_img_fltr_3);
        free(s_img_fltr_3_tmp);
        free(s_kernel3);
        free(s_img_fltr_4_tmp);

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
            s_weights_layer5[i]= SHR(shared, bc, party, role);
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
            s_biases_layer5[i] = SHR(shared, bc, party, role);
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
        // double *img_fltr_5 = (double *)calloc(rows * cols * num_filters5 , sizeof(double));
        // double *kernel5 = (double *)malloc(filtersize5*sizeof(double));
        // double *img_fltr_p5 = img_fltr_5; // Pointer to img_fltr5
        // double *img_fltr_5_tmp = (double *)malloc(rows * cols * sizeof(double));
        SHR *s_kernel5 = (SHR *)malloc(filtersize5 * sizeof(SHR));;
        SHR *s_img_fltr_5 = (SHR *)malloc(inRows * inCols * num_filters5*sizeof(SHR));
        SHR *s_img_fltr_p5 = s_img_fltr_5;
        SHR *s_img_fltr_5_tmp = (SHR *)malloc(inRows * inCols * sizeof(SHR));
        zero = 0.0;
        for(int i=0;i<inRows * inCols * num_filters5;i++){
            uint64_t *shared;
            if(role==SERVER){
                shared=new uint64_t(*(uint64_t*)&zero);
            }
            else{
                shared=NULL;
            }
            s_img_fltr_p5[i]= SHR(shared, bc, party, role);
            delete shared;
        }

        cnt_weight = 0;

        for (int i = 0; i < num_filters5; i++)
        {
            s_img_fltr_p4 = s_img_fltr_4;
            for (int j = 0; j < num_channels5; j++)
            {
                // reading corresponding weights to kernel
                for (int cnt_kernel = 0; cnt_kernel < filtersize5; cnt_kernel++)
                {
                    *(s_kernel5 + cnt_kernel) = s_weights_layer5[cnt_weight + cnt_kernel];
                }

                s_imfilter(s_img_fltr_p4, s_kernel5, s_img_fltr_5_tmp, inRows, inCols, padsize5);
                s_imadd(s_img_fltr_p5, s_img_fltr_5_tmp, inCols, inRows);
                cnt_weight = cnt_weight + filtersize5;
                s_img_fltr_p4 = s_img_fltr_p4 + inRows*inCols;
            }
            s_bias_tmp = s_biases_layer5[i];
            s_PReLU(s_img_fltr_p5, inRows, inCols, s_bias_tmp, prelu_coeff_layer5);
            s_img_fltr_p5 = s_img_fltr_p5 + inRows*inCols;
        }

        free(s_img_fltr_4);
        free(s_img_fltr_4_tmp);
        free(s_kernel4);
        free(s_img_fltr_5_tmp);

        /////////// Layer6
        // Reading weights of 6th layer
        FILE *weights_layer6_ptr;
        weights_layer6_ptr = fopen("weights_layer6.txt", "r");
        if (weights_layer6_ptr == NULL) { printf("Error in the reading weights of 6th layer\n"); };
        double weights_layer6[1296];
        SHR s_weights_layer6[1296];
        for (int i = 0; i < 1296; i++)
        {
            uint64_t *shared;
            if(role==SERVER){
                fscanf(weights_layer6_ptr, "%lf", &weights_layer6[i]);
                shared=new uint64_t(*(uint64_t*)&weights_layer6[i]);
            }
            else{
                shared=NULL;
            }
            s_weights_layer6[i]= SHR(shared, bc, party, role);
            delete shared;
        }
        fclose(weights_layer6_ptr);
        // Reading biases of 6th layer
        FILE *biases_layer6_ptr;
        biases_layer6_ptr = fopen("biasess_layer6.txt", "r");
        if (biases_layer6_ptr == NULL) { printf("Error in the reading biases of 6th layer\n"); };
        double biases_layer6[12];
        SHR s_biases_layer6[12];
        for (int i = 0; i < 12; i++)
        {
            uint64_t *shared;
            if(role==SERVER){
                fscanf(biases_layer6_ptr, "%lf", &biases_layer6[i]);
                shared=new uint64_t(*(uint64_t*)&biases_layer6[i]);
            }
            else{
                shared=NULL;
            }
            s_biases_layer6[i] = SHR(shared, bc, party, role);
            delete shared;
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
        // double *img_fltr_6 = (double *)calloc(rows * cols * num_filters6 , sizeof(double));
        // double *kernel6 = (double *)malloc(filtersize6*sizeof(double));
        // double *img_fltr_p6 = img_fltr_6; // Pointer to img_fltr6
        // double *img_fltr_6_tmp = (double *)malloc(rows * cols * sizeof(double));
        SHR *s_kernel6 = (SHR *)malloc(filtersize6 * sizeof(SHR));;
        SHR *s_img_fltr_6 = (SHR *)malloc(inRows * inCols * num_filters6*sizeof(SHR));
        SHR *s_img_fltr_p6 = s_img_fltr_6;
        SHR *s_img_fltr_6_tmp = (SHR *)malloc(inRows * inCols * sizeof(SHR));
        zero = 0.0;
        for(int i=0;i<inRows * inCols * num_filters6;i++){
            uint64_t *shared;
            if(role==SERVER){
                shared=new uint64_t(*(uint64_t*)&zero);
            }
            else{
                shared=NULL;
            }
            s_img_fltr_p6[i]= SHR(shared, bc, party, role);
            delete shared;
        }

        cnt_weight = 0;

        for (int i = 0; i < num_filters6; i++)
        {
            s_img_fltr_p5 = s_img_fltr_5;
            for (int j = 0; j < num_channels6; j++)
            {
                // reading corresponding weights to kernel
                for (int cnt_kernel = 0; cnt_kernel < filtersize6; cnt_kernel++)
                {
                    *(s_kernel6 + cnt_kernel) = s_weights_layer6[cnt_weight + cnt_kernel];
                }

                s_imfilter(s_img_fltr_p5, s_kernel6, s_img_fltr_6_tmp, inRows, inCols, padsize6);
                s_imadd(s_img_fltr_p6, s_img_fltr_6_tmp, inCols, inRows);

                cnt_weight = cnt_weight + filtersize6;
                s_img_fltr_p5 = s_img_fltr_p5 + inRows*inCols;
            }
            s_bias_tmp = s_biases_layer6[i];
            s_PReLU(s_img_fltr_p6, inRows, inCols, s_bias_tmp, prelu_coeff_layer6);
            s_img_fltr_p6 = s_img_fltr_p6 + inRows*inCols;
        }

        free(s_img_fltr_5);
        free(s_img_fltr_5_tmp);
        free(s_kernel5);

        /////////// Layer7
        // Reading weights of 7th layer
        FILE *weights_layer7_ptr;
        weights_layer7_ptr = fopen("weights_layer7.txt", "r");
        if (weights_layer7_ptr == NULL) { printf("Error in the reading weights of 7th layer\n"); };
        double weights_layer7[672];
        SHR s_weights_layer7[672];
        for (int i = 0; i < 672; i++)
        {
            uint64_t *shared;
            if (role == SERVER)
            {
                fscanf(weights_layer7_ptr, "%lf", &weights_layer7[i]);
                shared = new uint64_t(*(uint64_t *)&weights_layer7[i]);
            }
            else
            {
                shared = NULL;
            }
            s_weights_layer7[i] = SHR(shared, bc, party, role);
            delete shared;
        }
        fclose(weights_layer7_ptr);

        // Reading biases of 7th layer
        FILE *biases_layer7_ptr;
        biases_layer7_ptr = fopen("biasess_layer7.txt", "r");
        if (biases_layer7_ptr == NULL) { printf("Error in the reading biases of 7th layer\n"); };
        double biases_layer7[56];
        SHR s_biases_layer7[56];
        for (int i = 0; i < 56; i++)
        {
            uint64_t *shared;
            if (role == SERVER)
            {
                fscanf(biases_layer7_ptr, "%lf", &biases_layer7[i]);
                shared = new uint64_t(*(uint64_t *)&biases_layer7[i]);
            }
            else
            {
                shared = NULL;
            }
            s_biases_layer7[i] = SHR(shared, bc, party, role);
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
        // double *img_fltr_7 = (double *)calloc(rows * cols * num_filters7 , sizeof(double));
        // double *kernel7 = (double *)malloc(filtersize7*sizeof(double));
        // double *img_fltr_p7 = img_fltr_7; // Pointer to img_fltr7
        // double *img_fltr_7_tmp = (double *)malloc(rows * cols * sizeof(double));
        SHR *s_img_fltr_7 = (SHR *)malloc(inRows * inCols * num_filters7 * sizeof(SHR));
        SHR *s_kernel7 = (SHR*)malloc(filtersize7*sizeof(SHR));
        SHR *s_img_fltr_p7 = s_img_fltr_7;
        SHR *s_img_fltr_7_tmp = (SHR*)malloc(inRows * inCols * sizeof(SHR));

        zero = 0.0;
        for(int i=0;i<inRows * inCols * num_filters6;i++){
            uint64_t *shared;
            if(role==SERVER){
                shared=new uint64_t(*(uint64_t*)&zero);
            }
            else{
                shared=NULL;
            }
            *(s_img_fltr_p7+i)= SHR(shared, bc, party, role);
            delete shared;
        }

        cnt_weight = 0;

        for (int i = 0; i < num_filters7; i++)
        {
            s_img_fltr_p6 = s_img_fltr_6;
            for (int j = 0; j < num_channels7; j++)
            {
                // reading corresponding weights to kernel
                for (int cnt_kernel = 0; cnt_kernel < filtersize7; cnt_kernel++)
                {
                    *(s_kernel7 + cnt_kernel) = s_weights_layer7[cnt_weight + cnt_kernel];
                }

                s_imfilter(s_img_fltr_p6, s_kernel7, s_img_fltr_7_tmp, inRows, inCols, padsize7);
                s_imadd(s_img_fltr_p7, s_img_fltr_7_tmp, inCols, inRows);

                cnt_weight = cnt_weight + filtersize7;
                s_img_fltr_p6 = s_img_fltr_p6 + inRows*inCols;
            }
            s_bias_tmp = s_biases_layer7[i];
            s_PReLU(s_img_fltr_p7, inRows, inCols, s_bias_tmp, prelu_coeff_layer7);
            s_img_fltr_p7 = s_img_fltr_p7 + inRows*inCols;
        }

        free(s_img_fltr_6);
        free(s_img_fltr_6_tmp);
        free(s_kernel6);

        /////////// Convolution3 ------------------- Layer 8

        /////////// Layer8
        // Reading weights of 8th layer
        FILE *weights_layer8_ptr;
        weights_layer8_ptr = fopen("weights_layer8.txt", "r");
        if (weights_layer8_ptr == NULL) { printf("Error in the reading weights of 8th layer\n"); };
        double weights_layer8[4536];
        SHR* s_weights_layer8= (SHR*)calloc(4536 , sizeof(SHR));
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

        // Decvolution ==> output is the img_hr
        // double *img_fltr_8 = (double *)calloc((rows*scale) *(cols*scale) * num_filters8 , sizeof(double));
        // double *kernel8 = (double *)malloc(filtersize8*sizeof(double));
        // double *img_fltr_8_tmp = (double *)malloc((rows*scale) *(cols*scale) * sizeof(double));
        SHR *s_img_fltr_8 = (SHR*)calloc(101376 , sizeof(SHR));
	    SHR *s_kernel8 = (SHR*)malloc(filtersize8*sizeof(SHR));
	    SHR *s_img_fltr_8_tmp = (SHR*)calloc(101376 , sizeof(SHR));

        cnt_weight = 0;
        s_img_fltr_p7 = s_img_fltr_7;

        for (int j = 0; j < num_channels8; j++)
        {
            // reading corresponding weights to kernel
            for (int cnt_kernel = 0; cnt_kernel < filtersize8; cnt_kernel++)
            {
                SHR x = *(s_weights_layer8+cnt_weight + cnt_kernel);
			    *(s_kernel8+cnt_kernel) = x;
            }
            cnt_weight = cnt_weight + filtersize8;

            s_deconv(s_img_fltr_p7, s_img_fltr_8_tmp, s_kernel8, inCols, inRows, scale);
		    s_imadd(s_img_fltr_8, s_img_fltr_8_tmp, inCols*scale, inRows*scale);
            
            s_img_fltr_p7 = s_img_fltr_p7 + inRows*inCols;
        }
        SHR *s_outBuf_tmp = (SHR*)malloc(outCols*outRows*sizeof(SHR));
	    SHR *s_outP_tmp = s_outBuf_tmp;
        for (int i=0;i<inRows*scale;i++)
        for (int j = 0;j<inCols*scale; j++)
        {
            int cnt_fnl = i*inCols*scale + j;
            //*(img_hr + cnt_fnl) = *(img_fltr_8 + cnt_fnl) + biases_layer8;
            *(s_outP_tmp + cnt_fnl) = *(s_img_fltr_8 + cnt_fnl) + biases_layer8;
		    *(outP_tmp + cnt_fnl) = s_outP_tmp[cnt_fnl].check();
        }

        /*for ( int i = 0; i < 10; i++)
        {
            printf("%f\n", *(img_hr + i));
        }*/

        free(s_img_fltr_7);
        free(s_img_fltr_7_tmp);
        free(s_kernel7);

        free(s_img_fltr_8);
        free(s_img_fltr_8_tmp);
        free(s_kernel8);
//====== FSRCNN end

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
		fread(inBuf, sizeof(unsigned char), inCols*inRows, inFp);
	    fread(inBuf, sizeof(unsigned char), inCols*inRows / 4, inFp);
        
		inP = inBuf;
		outP = outBuf;

		for (int i = 0; i < inRows / 2; i++)
		for (int j = 0; j < inCols / 2; j++) {

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
	free(inBuf);
	inBuf = NULL;
	free(s_inBuf_tmp);
	free(outBuf);
	outBuf = NULL;
	free(outBuf_tmp);
	outBuf_tmp = NULL;
    //free(s_outBuf_tmp);
}