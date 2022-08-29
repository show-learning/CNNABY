#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <vector>

using namespace std;

ABYParty *party;
BooleanCircuit *bc;
ArithmeticCircuit *ac;
e_role role;

void buildparty(std::string address, std::string p, std::string r, std::string path)
{
	// default argument
	uint16_t port = (uint16_t)std::stoi(p);
	role = (r == "SERVER") ? SERVER : CLIENT;
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

void imfilter(double *img, double *kernel, double *img_fltr, int rows, int cols, int padsize);
void pad_image(double *img, double *img_pad, int rows, int cols, int padsize);
void PReLU(double *img_fltr, int rows, int cols, double bias, double prelu_coeff);
double Max(double a, double b);
double Min(double a, double b);
void imadd(double *img_fltr_crnt, double *img_fltr_prev, int cols, int rows);

// The parameter aᵢ is the coefficient of the negative part which is learnable. PReLU is used to avoid dead features.
// Connecting these five parts forms a eight-layer network as:
//  Conv(5,56,1) — PReLU — Conv(1,12,56) — PReLU — 4 × {Conv(3,12,12) — PReLU} — Conv(1,56,12) — PReLU — DeConv(9,1,12).
// For more information about the design details of L1, please refer to the original paper.

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
			cnt = (i - padsize) * cols + (j - padsize); // counter which shows current pixel in filtered image (central pixel in convolution window)
			sum = 0;
			cnt_krnl = 0; // counter which determines kernel elements
			for (k1 = -padsize; k1 <= padsize; k1++)
				for (k2 = -padsize; k2 <= padsize; k2++)
				{
					cnt_pad = (i + k1) * cols_pad + j + k2; // counter which shows each neighbouring pixel of padded image used for convolution with kernel
					sum = sum + (*(img_pad + cnt_pad)) * (*(kernel + cnt_krnl));
					cnt_krnl++;
				}
			*(img_fltr + cnt) = sum;
		}

	free(img_pad);
	img_pad = NULL;
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
			cnt = (i - padsize) * (cols) + j - padsize;
			double x = *(img + cnt);
			*(img_pad + cnt_pad) = x;
		}
	// Top and Bottom Rows
	for (j = padsize; j < cols_pad - padsize; j++)
		for (k = 0; k < padsize; k++)
		{
			// Top Rows
			cnt_pad = j + k * cols_pad;
			cnt = j - padsize;
			*(img_pad + cnt_pad) = *(img + cnt);
			// Bottom Rows
			cnt_pad = j + (rows_pad - 1 - k) * cols_pad;
			cnt = (j - padsize) + (rows - 1) * cols;
			*(img_pad + cnt_pad) = *(img + cnt);
		}
	// Left and Right Columns
	for (i = padsize; i < rows_pad - padsize; i++)
		for (k = 0; k < padsize; k++)
		{
			// Left Columns
			cnt = (i - padsize) * cols;
			cnt_pad = i * cols_pad + k;
			*(img_pad + cnt_pad) = *(img + cnt);
			// Right Columns
			cnt = (i - padsize) * cols + cols - 1;
			cnt_pad = i * cols_pad + cols_pad - 1 - k;
			*(img_pad + cnt_pad) = *(img + cnt);
		}
	// Corner Pixels
	for (k1 = 0; k1 < padsize; k1++)
		for (k2 = 0; k2 < padsize; k2++)
		{
			// Upper Left Corner
			cnt_pad = k1 * cols_pad + k2;
			*(img_pad + cnt_pad) = *(img);
			// Upper Right Corner
			cnt_pad = k1 * cols_pad + cols_pad - 1 - k2;
			*(img_pad + cnt_pad) = *(img + cols - 1);
			// Lower Left Corner
			cnt_pad = (rows_pad - 1 - k1) * cols_pad + k2;
			*(img_pad + cnt_pad) = *(img + (rows - 1) * cols);
			// Lower Right Corner
			cnt_pad = (rows_pad - 1 - k1) * cols_pad + cols_pad - 1 - k2;
			*(img_pad + cnt_pad) = *(img + (rows - 1) * cols + cols - 1);
		}
}

void PReLU(double *img_fltr, int rows, int cols, double bias, double prelu_coeff)
{
	int cnt = 0;
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
		{
			cnt = i * cols + j;
			*(img_fltr + cnt) = Max(*(img_fltr + cnt) + bias, 0) + prelu_coeff * Min(*(img_fltr + cnt) + bias, 0);
		}
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

void imadd(double *img_fltr_sum, double *img_fltr_crnt, int cols, int rows)
{
	// *img_fltr_crnt ==> pointer to current feature map
	// *img_fltr_sum ==> pointer to the cumulutive feature map

	int cnt = 0;
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
		{
			cnt = i * cols + j;
			*(img_fltr_sum + cnt) = *(img_fltr_sum + cnt) + *(img_fltr_crnt + cnt);
		}
}

int main(int argc, char **argv)
{
	int inCols = 176; // Width of input (downsampled) video
	int inRows = 144; // Height of input (downsampled) video
	// Upsampler parameters
	int scale = 2;
	// Compressed Assault Cube
	int num = 1; // Number of frames to interpolat
	std::string address(argv[1]), port(argv[2]), r(argv[3]), path(argv[4]);
	e_role role = (r == "SERVER") ? SERVER : CLIENT;
	buildparty(address, port, role, path);

	FILE *inFp;
	inFp = fopen("input.yuv", "rb");
	if (inFp == NULL)
	{
		printf("\n We have null pointer \n");
	}

	// To read and write each frame in an unsigned character format
	unsigned char *inBuf = (unsigned char *)malloc(inCols * inRows * sizeof(unsigned char));
	// unsigned char *outBuf = (unsigned char *)malloc(outCols * outRows * sizeof(unsigned char));
	// To work with each pixel in the range of 0~1

	double *inBuf_tmp = (double *)malloc(inCols * inRows * sizeof(double));
	// SHR start
	SHR *s_inBuf_tmp = (SHR *)calloc(inCols * inRows * sizeof(SHR));
	//
	// double *outBuf_tmp = (double *)malloc(outCols * outRows * sizeof(double));

	for (int fcnt = 0; fcnt < num; fcnt++)
	{
		//////// Interpolate each frame using FSRCNN for Y component and simple repitition for U and V components
		// Pointer to obtain value of each tpixel of input frame
		unsigned char *inP = inBuf;
		double *inP_tmp = inBuf_tmp;
		// SHR start
		SHR *s_inP_tmp = s_inBuf_tmp;

		// Y Component
		fread(inBuf, sizeof(unsigned char), inCols * inRows, inFp);
		int i, j;
		for (i = 0; i < inRows; i++)
		{
			for (j = 0; j < inCols; j++)
			{
				int cnt = i * inCols + j;
				int x = *inP++;
				*(inP_tmp + cnt) = (double)(x / 255.0);
				// SHR assign
				uint64_t *shared;
				if (role == CLIENT)
				{
					shared = new uint64_t(*(uint64_t *)&inP_tmp[cnt]);
				}
				else
				{
					shared = NULL;
				}
				s_inP_tmp[cnt] = SHR(shared, bc, party, role);
				delete shared;
			}
		}

		/////////// Convolution1 -------- Layer1
		// Reading weights of first layer
		FILE *weights_layer1_ptr;
		weights_layer1_ptr = fopen("weights_layer1.txt", "r");
		if (weights_layer1_ptr == NULL)
		{
			printf("Error in the reading weights of first layer\n");
		};
		double weights_layer1[1400];
		// SHR start_weight
		SHR s_weights_layer1[1400];
		for (int i = 0; i < 1400; i++)
		{
			uint64_t *shared;
			if (role == SERVER)
			{
				scanf(weights_layer1_ptr, "%lf", &weights_layer1[i]);
				shared = new uint64_t(*(uint64_t *)&weights_layer1[i]);
			}
			else
			{
				shared = NULL;
			}
			s_weights_layer1[i] = SHR(shared, bc, party, role);
			delete shared;
		}
		fclose(weights_layer1_ptr);

		// Reading biases of first layer
		FILE *biases_layer1_ptr;
		biases_layer1_ptr = fopen("biasess_layer1.txt", "r");
		if (biases_layer1_ptr == NULL)
		{
			printf("Error in the reading biases of first layer\n");
		};
		double biases_layer1[56];
		// SHR start_bias
		SHR s_biases_layer5[56];
		for (int i = 0; i < 56; i++)
		{
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
		int filtersize = 25; // 5X5
		int patchsize = 5;
		int padsize = (patchsize - 1) / 2;
		int num_filters = 56;
		double prelu_coeff_layer1 = -0.8986;

		// Convolution
		double *img_fltr_1 = (double *)malloc(inRows * inCols * num_filters * sizeof(double));
		double *kernel = (double *)malloc(filtersize * sizeof(double));
		double *img_fltr_p1 = img_fltr_1; // Pointer to img_fltr1 ==>> Using this way to be able to shift it to access data

		// SHR start
		SHR *s_img_fltr_1 = (SHR *)calloc(rows * cols * num_filters1, sizeof(SHR));
		SHR *s_img_fltr_p1 = s_img_fltr_1;

		int cnt_weight = 0;
		double bias_tmp;
		for (int i = 0; i < num_filters; i++)
		{
			// reading corresponding weights to kernel pointer
			for (int cnt_kernel = 0; cnt_kernel < filtersize; cnt_kernel++)
			{
				*(kernel + cnt_kernel) = s_weights_layer1[cnt_weight + cnt_kernel];
			}
			cnt_weight = cnt_weight + filtersize;
			imfilter(s_inP_tmp, kernel, s_img_fltr_p1, inRows, inCols, padsize);

			bias_tmp = biases_layer1[i];
			PReLU(s_img_fltr_p1, inRows, inCols, bias_tmp, prelu_coeff_layer1);

			s_img_fltr_p1 = s_img_fltr_p1 + inCols * inRows;
		}

		FILE *outFp = fopen("L1_output.txt", "wb");
		for (int i = 0; i < inRows * inCols * num_filters; i++)
		{
			fprintf(outFp, "%lf ", s_img_fltr_1[i].get_d());
		}
		fclose(outFp);
		free(inBuf);
		inBuf = NULL;
		free(inBuf_tmp);
		inBuf_tmp = NULL;
		// free(outBuf);
		// outBuf = NULL;
		// free(outBuf_tmp);
		// outBuf_tmp = NULL;
	}
}