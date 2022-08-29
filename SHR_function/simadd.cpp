void imadd(SHR *img_fltr_sum, SHR *img_fltr_crnt, int cols, int rows)
//void imadd(double *img_fltr_sum, double *img_fltr_crnt, int cols, int rows)
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