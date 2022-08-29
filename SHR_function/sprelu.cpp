SHR Max(SHR a, SHR b)
{
    SHR result=(a>b).mux(a,b);//如果cond成立執行回傳第一個參數,否則回傳第二個參數
    return result;
}

SHR Min(SHR a, SHR b)
{
    SHR result=(a>b).mux(b,a);//如果cond成立執行回傳第一個參數,否則回傳第二個參數
    return result;
}

void PReLU(SHR *img_fltr,int rows, int cols, SHR bias, SHR prelu_coeff)
//void PReLU(double *img_fltr,int rows, int cols, double bias, double prelu_coeff)
{
	int cnt = 0;
	for (int i = 0; i < rows;i++)
	for (int j = 0; j < cols; j++)
	{
		cnt = i*cols + j;
		*(img_fltr + cnt) = Max(*(img_fltr + cnt) + bias, 0) + prelu_coeff * Min(*(img_fltr + cnt) + bias, 0);
	}
}

// double Max(double a, double b)
// {
// 	double c;
// 	c = a > b ? a : b;
// 	return c;
// }



// double Min(double a, double b)
// {
// 	double c;
// 	c = a > b ? b : a;
// 	return c;
// }

