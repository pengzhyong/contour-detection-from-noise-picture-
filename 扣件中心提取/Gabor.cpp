#include "Gabor.h"
#include <iostream>


Gabor::Gabor()
{
}


Gabor::~Gabor()
{
}

//void Gabor::Gabor_kernel(Mat& src, Mat& dst, int size, double theta, double sigma_x, double sigma_y,double x, double y)
//{
//	double x1,y1;
//	double lamda = 3;
//	double phi = 0;
//	double gama = 0.1;
//	double sigma = 0.56*lamda;
//	Mat GaborReal(size, size, CV_32FC1);
//	Mat GaborImag(size, size, CV_32FC1);
//	double F = 0.2;
//	double sumReal = 0;
//	double sumImag = 0;
//
//	for (int i=-size/2;i<=size/2;i++)
//	{
//		float * pReal = GaborReal.ptr<float>(i + size / 2);
//		float* pImag = GaborImag.ptr<float>(i + size / 2);
//
//		for (int j=-size/2;j<=size/2;j++)
//		{
//			x1 = i*cos(theta) + j*sin(theta);
//			y1 = -i*sin(theta) + j*cos(theta);			
//			double real = exp(-(x1*x1 + gama*gama*y1*y1) / (2 * sigma*sigma))*cos(2 * 3.14159*x1 / lamda + phi);
//			//double real = exp(-0.5*((x1*x1 / (sigmax*sigmax)) + (y1*y1 / (sigmay*sigmay)))) / (2 * 3.14159*sigmax*sigmax)*cos(2*3.14159*F*x1);
//			double imag = exp(-(x1*x1 + gama*gama*y1*y1) / (2 * sigma*sigma))*sin(2 * 3.14159*x1 / lamda + phi);
//
//			pReal[j + size / 2] =  real;// sqrt(real*real + imag*imag);
//			pImag[j + size / 2] = imag;// sqrt(real*real + imag*imag);
//			sumReal += pReal[j + size / 2];
//			sumImag += pImag[j + size / 2];
//
//			//=====一组12个Gabor滤波器
//			for (int k=1;k<12;k++)
//			{
//				double deltaA = 3.14159 / 12;
//				x1 = i*cos(theta+k*deltaA) + j*sin(theta+ k*deltaA);
//				y1 = -i*sin(theta+ k*deltaA) + j*cos(theta+ k*deltaA);
//				real = exp(-(x1*x1 + gama*gama*y1*y1) / (2 * sigma*sigma))*cos(2 * 3.14159*x1 / lamda + phi);
//				//real = exp(-0.5*((x1*x1 / (sigmax*sigmax)) + (y1*y1 / (sigmay*sigmay)))) / (2 * 3.14159*sigmax*sigmax)*cos(2*3.14159*F*x1);
//				imag = exp(-(x1*x1 + gama*gama*y1*y1) / (2 * sigma*sigma))*sin(2 * 3.14159*x1 / lamda + phi);
//
//				pReal[j + size / 2] += real;// sqrt(real*real + imag*imag);
//				pImag[j + size / 2] += imag;// sqrt(real*real + imag*imag);
//				sumReal += pReal[j + size / 2];
//				sumImag += pImag[j + size / 2];
//			}
//
//
//		}
//	}	
//	GaborReal = GaborReal / sumReal;
//	GaborImag = GaborImag / sumImag;
//	Mat dst1(src.size(), CV_32FC1);
//	Mat dst2(src.size(), CV_32FC1);
//
//	filter2D(src, dst1, dst2.depth(), GaborReal);
//	//std::cout << dst1 << std::endl;
//	//return;
//	filter2D(src, dst2, dst2.depth(), GaborImag);	
//	//normalize(dst1, dst1, 0, 255, NORM_MINMAX);
//	//dst1.convertTo(dst1, CV_8U);
//	//normalize(dst2, dst2, 0, 1, NORM_MINMAX);
//
//	imshow("dst1", dst1);
//	imshow("dst2", dst2);
//
//	for (int i = 0; i < dst2.rows; i++)
//	{
//		float* ptr = dst.ptr<float>(i);
//		float* ptr1 = dst1.ptr<float>(i);
//		float* ptr2 = dst2.ptr<float>(i);
//
//		for (int j = 0; j < dst2.cols; j++)
//		{
//			ptr[j] = sqrt(ptr1[j] * ptr1[j] + ptr2[j] * ptr2[j]);
//		}
//	}
//
//	//normalize(dst, dst, 0, 1, NORM_MINMAX);
//	//std::cout << dst << std::endl;
//	//GaussianBlur(src, dst, Size(5,5),3,3);
//}

void Gabor::Gabor_kernel(Mat& src, Mat& dst, double theta, double lamda1, double gama1, int size, double phi1,  double sigmaRatio)
{
	double x1, y1;
	double lamda = lamda1;
	double phi = phi1;
	double gama = gama1;
	double sigma = sigmaRatio*lamda;
	Mat GaborReal(size, size, CV_32FC1);
	Mat GaborImag(size, size, CV_32FC1);
	double F = 0.2;
	double sumReal = 0;
	double sumImag = 0;

	for (int i = -size / 2; i <= size / 2; i++)
	{
		float * pReal = GaborReal.ptr<float>(i + size / 2);
		float* pImag = GaborImag.ptr<float>(i + size / 2);

		for (int j = -size / 2; j <= size / 2; j++)
		{
			x1 = i*cos(theta) + j*sin(theta);
			y1 = -i*sin(theta) + j*cos(theta);
			double real = exp(-(x1*x1 + gama*gama*y1*y1) / (2 * sigma*sigma))*cos(2 * 3.14159*x1 / lamda + phi);
			//double real = exp(-0.5*((x1*x1 / (sigmax*sigmax)) + (y1*y1 / (sigmay*sigmay)))) / (2 * 3.14159*sigmax*sigmax)*cos(2*3.14159*F*x1);
			double imag = exp(-(x1*x1 + gama*gama*y1*y1) / (2 * sigma*sigma))*sin(2 * 3.14159*x1 / lamda + phi);

			pReal[j + size / 2] = real;// sqrt(real*real + imag*imag);
			pImag[j + size / 2] = imag;// sqrt(real*real + imag*imag);
			sumReal += pReal[j + size / 2];
			sumImag += pImag[j + size / 2];

			//=====一组12个Gabor滤波器
			for (int k = 1; k < 3; k++)
			{
				double deltaA = 3.14159 / 3;
				x1 = i*cos(theta + k*deltaA) + j*sin(theta + k*deltaA);
				y1 = -i*sin(theta + k*deltaA) + j*cos(theta + k*deltaA);
				real = exp(-(x1*x1 + gama*gama*y1*y1) / (2 * sigma*sigma))*cos(2 * 3.14159*x1 / lamda + phi);
				//real = exp(-0.5*((x1*x1 / (sigmax*sigmax)) + (y1*y1 / (sigmay*sigmay)))) / (2 * 3.14159*sigmax*sigmax)*cos(2*3.14159*F*x1);
				imag = exp(-(x1*x1 + gama*gama*y1*y1) / (2 * sigma*sigma))*sin(2 * 3.14159*x1 / lamda + phi);
				//imag = exp(-(x1*x1 + gama*gama*y1*y1) / (2 * sigma*sigma))*cos(2 * 3.14159*x1 / lamda + phi + 3.14159/2);


				//pReal[j + size / 2] += real;// sqrt(real*real + imag*imag);
				//pImag[j + size / 2] += imag;// sqrt(real*real + imag*imag);
				sumReal += pReal[j + size / 2];
				sumImag += pImag[j + size / 2];
			}


		}
	}
	//GaborReal = GaborReal / sumReal;
	//GaborImag = GaborImag / sumImag;
	Mat dst1(src.size(), CV_32FC1);
	Mat dst2(src.size(), CV_32FC1);

	filter2D(src, dst1, dst2.depth(), GaborReal);
	//std::cout << dst1 << std::endl;
	//return;
	filter2D(src, dst2, dst2.depth(), GaborImag);
	//normalize(dst1, dst1, 0, 255, NORM_MINMAX);
	//dst1.convertTo(dst1, CV_8U);
	//normalize(dst2, dst2, 0, 1, NORM_MINMAX);

	imshow("gaborReal", dst1);
	imshow("gaborImag", dst2);

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