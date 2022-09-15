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
#include <bits/stdc++.h>
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

//use check to get real value
int main(int argc, char **argv)
{
	int inCols = 176; // Width of input (downsampled) video
	int inRows = 144; // Height of input (downsampled) video
	// Upsampler parameters
	int scale = 2;
	// Compressed Assault Cube
	std::string address(argv[1]), port(argv[2]), r(argv[3]), path(argv[4]);
	e_role role = (r == "SERVER") ? SERVER : CLIENT;
	buildparty(address, port, r, path);

	FILE *inFp;
	inFp = fopen("L1_output.txt", "r");
	if (inFp == NULL)
	{
		printf("\n We have null pointer \n");
	}

	// To read and write each frame in an unsigned character format
	// unsigned char *inBuf = (unsigned char *)malloc(inCols * inRows * sizeof(unsigned char));
	// unsigned char *outBuf = (unsigned char *)malloc(outCols * outRows * sizeof(unsigned char));

	// SHR start
	// double *inBuf_tmp = (double *)malloc(inCols * inRows * sizeof(double));
	SHR *s_input = (SHR *)malloc(inCols * inRows *sizeof(SHR));
	double *img = (double *)malloc(inCols * inRows*sizeof(double));
	SHR *s_inP_tmp = s_input;
	for (int i = 0; i < inRows * inCols; i++)
	{
		double x;
		fscanf(inFp, "%lf", &x);
		s_inP_tmp[i] = SHR((uint64_t *)&x, bc, party);
        img[i] = s_inP_tmp[i].check();
	}
	//////// Interpolate each frame using FSRCNN for Y component and simple repitition for U and V components
	// Pointer to obtain value of each tpixel of input frame
	FILE* outFp = fopen("input.txt", "w");
	for (int i = 0; i < inRows * inCols; i++)
	{
		fprintf(outFp, "%lf\n", img[i]);
	}
	fclose(outFp);
	free(s_input);
	
}