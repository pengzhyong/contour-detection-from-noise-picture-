#include "Gabor.h"
#include <iostream>


Gabor::Gabor()
{
}


Gabor::~Gabor()
{
}

void Gabor::Gabor_kernel(Mat& src, Mat& dst, int size, double theta, double sigma_x, double sigma_y,double x, double y)
{
	setPara(x, y, theta, sigma_x, sigma_y, size);
	double x1,y1;
	double lamda = 3;
	double phi = 0;
	double gama = 0.2;
	double sigma = 0.56*lamda;
	Mat GaborReal(size, size, CV_32FC1);
	Mat GaborImag(size, size, CV_32FC1);
	double F = 0.2;
	double sum = 0;
	for (int i=-size/2;i<=size/2;i++)
	{
		float * pReal = GaborReal.ptr<float>(i + size / 2);
		float* pImag = GaborImag.ptr<float>(i + size / 2);

		for (int j=-size/2;j<=size/2;j++)
		{
			x1 = i*cos(theta) + j*sin(theta);
			y1 = -i*sin(theta) + j*cos(theta);			
			double real = exp(-(x1*x1 + gama*gama*y1*y1) / (2 * sigma*sigma))*cos(2 * 3.14159*x1 / lamda + phi);
			//double real = exp(-0.5*((x1*x1 / (sigmax*sigmax)) + (y1*y1 / (sigmay*sigmay)))) / (2 * 3.14159*sigmax*sigmax)*cos(2*3.14159*F*x1);
			double imag = exp(-(x1*x1 + gama*gama*y1*y1) / (2 * sigma*sigma))*sin(2 * 3.14159*x1 / lamda + phi);

			pReal[j + size / 2] =  real;// sqrt(real*real + imag*imag);
			pImag[j + size / 2] = imag;// sqrt(real*real + imag*imag);
			sum += pReal[j + size / 2];
		}
	}	
	Mat dst1(src.size(), CV_32FC1);
	Mat dst2(src.size(), CV_32FC1);

	filter2D(src, dst1, dst.depth(), GaborReal);
	//return;
	filter2D(src, dst2, dst2.depth(), GaborImag);	
	for (int i = 0; i < dst2.rows; i++)
	{
		float* ptr = dst.ptr<float>(i);
		float* ptr1 = dst1.ptr<float>(i);
		float* ptr2 = dst2.ptr<float>(i);

		for (int j = 0; j < dst2.cols; j++)
		{
			ptr[j] = sqrt(ptr1[j] * ptr1[j] + ptr2[j] * ptr2[j]);
		}
	}

	//normalize(dst, dst, 0, 1, NORM_MINMAX);
	//std::cout << dst << std::endl;
	//GaussianBlur(src, dst, Size(5,5),3,3);
}

void Gabor::setPara(double x, double y, double angle , double sigma_x, double sigma_y, int s)
{
	double x0 = x;
	double y0 = y;
	double theta = angle;
	double sigmax = sigma_x;
	double sigmay = sigma_y;
	int size = s;
}
