#include<opencv2/opencv.hpp>
#include<vector>
#include "Gabor.h"
using namespace cv;
using namespace std;
const double PI = 3.14159;
int main()
{
	Mat src = imread("pic/gabor2.png");
	//Mat src = imread("pic/lm.jpg");

	//pyrDown(src, src);
	//pyrDown(src, src);
	imshow("src", src);
	Mat gray;
	cvtColor(src, gray, COLOR_BGR2GRAY);
	GaussianBlur(gray, gray, Size(3, 3), 1);
	Mat gray0;
	gray.copyTo(gray0);

	Mat gabor(gray.size(), CV_32FC1);
	vector<vector<Point> > contours;
	//vector<vector<Vec4i> > hierarchy;
	//Canny(gray, can, 50, 100);
	//=========
	Gabor gabor_ker;
	gabor_ker.Gabor_kernel(gray, gabor,15, PI*0/3);
	gabor.convertTo(gabor, CV_8U);
	imshow("gabor", gabor);
	Mat can;
	gabor.copyTo(can);
	can = can * 255;
	can.convertTo(can, CV_8UC1);
	//imshow("can0", can);
	//normalize(can, can, 0, 255);
	Canny(can, can, 100, 220);
	//imshow("can", can);
	waitKey(0);
}