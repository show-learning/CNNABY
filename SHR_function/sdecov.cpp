void s_deconv(SHR *img_input, SHR *img_output, SHR *kernel, int cols, int rows, int stride)
{
	int border = 1;
	int fsize = 9;
	int rows_pad = rows + 2 * border;
	int cols_pad = cols + 2 * border;
	SHR *s_img_input_padded = (SHR *)malloc(rows_pad * cols_pad * sizeof(SHR));
	s_pad_image(img_input, s_img_input_padded, rows, cols, border);
	int rows_out_pad = rows_pad * stride;
	int cols_out_pad = cols_pad * stride;
	SHR *s_kernel_modif = (SHR *)malloc(fsize * fsize * sizeof(SHR));
	SHR *s_img_output_tmp = (SHR *)calloc((rows_out_pad + fsize - 1)* (cols_out_pad + fsize - 1), sizeof(SHR));
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
			*(s_kernel_modif + cnt_kernel) = (*(kernel + cnt_kernel))*(*(s_img_input_padded + cnt_img));
			*(s_img_output_tmp + cnt_img_output + k_c) = ((*(s_img_output_tmp + cnt_img_output + k_c)) + (*(s_kernel_modif + cnt_kernel)));
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
		int cnt_img_out_tmp = i_tmp*(cols_out_pad + fsize - 1) + j_tmp; 
        // (cols-pad+fsize-1) is the number of columns in the img_out_tmp
		*(img_output + cnt_img_out) = *(s_img_output_tmp + cnt_img_out_tmp);
	}

	free(s_img_input_padded); s_img_input_padded = NULL;
	free(s_img_output_tmp); s_img_output_tmp = NULL;
	free(s_kernel_modif); s_kernel_modif = NULL;
}