#pragma once
#include<opencv2/opencv.hpp>
using namespace cv;

class Gabor
{
public:
	Gabor();
	~Gabor();
	void Gabor_kernel(Mat& src, Mat& dst, int size = 5, double theta = 0, double sigma_x = 1, double sigma_y = 3, double x = 0, double y = 0);
private:
	double x0 = 0;
	double y0 = 0;
	double theta = 0;
	double sigmax = 0.56*5;
	double sigmay =2*sigmax;
	int size = 5;
	void setPara(double x, double y, double angle, double sigma_x, double sigma_y, int s);

};

