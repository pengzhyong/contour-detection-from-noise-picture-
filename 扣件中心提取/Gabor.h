#pragma once
#include<opencv2/opencv.hpp>
using namespace cv;

class Gabor
{
public:
	Gabor();
	~Gabor();
	//void Gabor_kernel(Mat& src, Mat& dst, int size = 5, double theta = 0, double sigma_x = 1, double sigma_y = 3, double x = 0, double y = 0);
	void Gabor_kernel(Mat& src, Mat& dst, double theta = 0, double lamda1 = 3, double gama1 = 0.5, int size = 11, double phi1 = 0, double sigmaRatio = 0.56);
private:
	double x0 = 0;
	double y0 = 0;
	double theta = 0;
	double sigmax = 0.56*5;
	double sigmay =2*sigmax;
	int size = 11;
};

