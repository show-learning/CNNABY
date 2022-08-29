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

    int cols = 176; // Width of input (downsampled) video
    int rows = 144; // Height of input (downsampled) video
    std::string address(argv[1]), port(argv[2]), r(argv[3]), path(argv[4]);
    e_role role = (r == "SERVER") ? SERVER : CLIENT;
    buildparty(address, port, role, path);

    FILE *inFp, *outFp;
    inFp = fopen("L2_output.txt", "r");
    int num_filters = 25;
    double *img_fltr_1 = (double *)calloc(rows * cols * num_filters, sizeof(double));
    double *img_fltr_p1 = img_fltr_1;
    // SHR start

    SHR *s_img_fltr_1 = (SHR *)calloc(rows * cols * num_filters * sizeof(SHR));
    SHR *s_img_fltr_p1 = s_img_fltr_1;
    for (int i = 0; i < rows * cols * num_filters; i++)
    {
        uint64_t *shared;
        if (role == SERVER)
        {
            fscanf(inFp, "%lf", (img_fltr_1 + i));
            shared = new uint64_t(*(uint64_t *)&s_img_fltr_1[i]);
        }
        else
        {
            shared = NULL;
        }
        s_img_fltr_1[i] = SHR(shared, bc, party, role);
        delete shared;
    }
    /////////// Layer7
    // Reading weights of 2nd layer
    FILE *weights_layer2_ptr;
    weights_layer2_ptr = fopen("weights_layer2.txt", "r");
    if (weights_layer2_ptr == NULL)
    {
        printf("Error in the reading weights of 2nd layer\n");
    };
    // Note: weights must be saved in a way which that corresponding weights of each channel can be read by pointer concept ==>> for this layer 12X56 matrix is reshaped to (12X56)*1 vector
    double weights_layer2[672];
    // SHR start_weight
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
    if (biases_layer2_ptr == NULL)
    {
        printf("Error in the reading biases of 2nd layer\n");
    };
    double biases_layer2[12];
    // SHR start
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
    int filtersize2 = 1; // 1X1
    int patchsize2 = 1;
    int padsize2 = (patchsize2 - 1) / 2;
    int num_filters2 = 12;
    int num_channels2 = 56;
    double prelu_coeff_layer2 = 0.3236;
    // Convolution
    double *img_fltr_2 = (double *)calloc(rows * cols * num_filters2, sizeof(double)); // use calloc to initialize all variables to zero
    double *img_fltr_2_tmp = (double *)malloc(rows * cols * sizeof(double));
    double *kernel2 = (double *)malloc(filtersize2 * sizeof(double));

    // SHR start
    SHR *s_img_fltr_2 = (SHR *)calloc(rows * cols * num_filters2, sizeof(SHR));
    SHR *s_img_fltr_2_tmp = (SHR *)malloc(rows * cols * sizeof(SHR));
    SHR *s_img_fltr_p2 = s_img_fltr_2;
    ////////////////////////////////////////////////////////////
    int cnt_weight = 0;

    for (int i = 0; i < num_filters2; i++)
    {
        // printf("C");
        s_img_fltr_p1 = s_img_fltr_1; // Return pointer to the first of array which contains feature map of previous layer
        for (int j = 0; j < num_channels2; j++)
        {

            // printf("%d ", j);
            // reading corresponding weights to kernel
            for (int cnt_kernel = 0; cnt_kernel < filtersize2; cnt_kernel++)
            {
                // SHR yet
                *(kernel2 + cnt_kernel) = s_weights_layer2[cnt_weight + cnt_kernel];
            }

            // printf("%d ", j);
            imfilter(s_img_fltr_p1, kernel2, s_img_fltr_2_tmp, rows, cols, padsize2);

            // printf("%d ", j);
            imadd(s_img_fltr_p2, s_img_fltr_2_tmp, cols, rows);

            cnt_weight = cnt_weight + filtersize2;
            s_img_fltr_p1 = s_img_fltr_p1 + rows * cols;
        }
        double bias_tmp;
        bias_tmp = biases_layer2[i];
        // SHR yet
        SHR bias_tmp;

        PReLU(s_img_fltr_p2, rows, cols, bias_tmp, prelu_coeff_layer2);
        s_img_fltr_p2 = s_img_fltr_p2 + rows * cols;
    }

    // printf("pre_done");
    outFp = fopen("L2_output.txt", "wb");
    for (int i = 0; i < rows * cols * num_filters2; i++)
    {
        fprintf(outFp, "%lf ", s_img_fltr_2[i].get_d());
    }
    fclose(outFp);
    // printf("done");
    free(img_fltr_1);
    img_fltr_1 = NULL;
    free(img_fltr_2);
    img_fltr_2 = NULL;
    free(img_fltr_2_tmp);
    img_fltr_2_tmp = NULL;
}