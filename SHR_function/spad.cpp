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