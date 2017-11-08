#pragma once
#include <opencv2/opencv.hpp>
#include <vector>
using namespace cv;
using namespace std;
class ContourDetecion
{
public:
	ContourDetecion(Mat&);
	~ContourDetecion();
	void GaborDirection(Mat& srcImg, Mat& dstImg, Mat& dirImg, double lamda1 = 3, double gama1 = 0.5, int size = 11, double sigmaRatio = 0.56);
	void GaborGenerator(Mat& gaborReal, Mat& gaborImag, double theta, double lamda1 = 3, double gama1 = 0.5, int size = 11, double sigmaRatio = 0.56);
	void GaborSuppression(Mat& srcImg, Mat& driImg, Mat& dstImg, int supSize = 31);
	void GaborEnergy(Mat& srcImg, vector<Mat >& vGaborE, double lamda1 = 3, double gama1 = 0.5, int size = 11, double sigmaRatio = 0.56);
	void FacilitationEnergy(vector<Mat >& vGaborE, vector<Mat >& faciEnergy, int size = 11);
	void FacProgress(Mat& srcImg, Mat& dstImg);
private:
	void MergeRI(Mat& srcImg, Mat& dstImg, Mat& real, Mat& Imag);
	Mat src;
	Mat dst;
	vector<Mat > vGaborEnergy;
};

