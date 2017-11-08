#include<opencv2/opencv.hpp>
#include<vector>
#include "Gabor.h"
#include "ContourDetecion.h"
using namespace cv;
using namespace std;
const double PI = 3.14159;
void mySobel(Mat& src, Mat& dst, Mat& dstx, Mat& dsty);
void NMS(Mat& srcx, Mat& srcy, Mat& dst, int tl, int th);
void DFS(const Point point, Mat& img, Mat& visited, int th, int* flag);
void LOG(Mat& src, Mat& dst, double sigma1, double sigma2);
void Theta(int, void*);
void Lamda(int, void*);
void Size_gabor(int, void*);
void Gama(int, void*);
void SR(int, void*);
void Phi(int, void*);
void on_Change();

int g_theta = 0;
int g_lamda = 3;
int g_gama10 = 10;
int g_size = 11;
int g_phi = 0;
int g_sigmaRatio100 = 56;
Mat gray;
const cv::String WindowsName("Gabor");
int betaSup100 = 10;
void changeBeta(int, void*);
Mat gaborFilt;// (gray.size(), CV_32F);
Mat suppImg;// (gray.size(), CV_32F);
Mat suppMat;// (gray.size(), CV_32F);
int main()
{
	//Mat src = imread("pic/gabor2.jpg");
	Mat src = imread("pic/faci1.jpg");

	//pyrDown(src, src);
	//pyrDown(src, src);
	imshow("src", src);
	cvtColor(src, gray, COLOR_BGR2GRAY);
	GaussianBlur(gray, gray, Size(3, 3), 1);
	//threshold(gray, gray, 50, 255, THRESH_BINARY_INV);
	//imshow("threshold", gray);

	//====LOG
	/*Mat log(gray.size(), gray.type());
	LOG(gray, log, 1, 2);
	threshold(log, log, 1, 255, THRESH_BINARY);	
	imshow("log", log);*/

	//====MyCanny
	/*Mat mySob;
	Mat Sob;
	Mat mySobx;
	Mat mySoby;
	mySobel(gray, mySob, mySobx, mySoby);
	Sobel(gray, Sob,gray.depth(),1,1);
	imshow("mySobelx", mySobx);
	imshow("mySobely", mySoby);
	imshow("mySobel", mySob);
	imshow("Sobel", Sob);
	Mat myCanny;
	NMS(mySobx, mySoby, myCanny, 50, 100);
	imshow("myCanny", myCanny);*/
	
	//=====Gabor

	Mat srcImg = gray.clone();
	srcImg.convertTo(srcImg, CV_32F);
	ContourDetecion gaborSup(srcImg);

	Mat faciMat(srcImg.size(), CV_32F);
	gaborSup.FacProgress(gray, faciMat);
	

	//Mat gabor(gray.size(), CV_32FC1);
	//Gabor gabor_ker;
	//namedWindow(WindowsName, 0);
	//createTrackbar("theta", WindowsName, &g_theta, 180, Theta);	
	//createTrackbar("lamda", WindowsName, &g_lamda, 51, Lamda);
	//createTrackbar("gama", WindowsName, &g_gama10, 10, Gama);
	//createTrackbar("phi", WindowsName, &g_phi, 180, Phi);
	//createTrackbar("sigmaRatio", WindowsName, &g_sigmaRatio100, 200, SR);
	//createTrackbar("kerner size", WindowsName, &g_size, 50, Size_gabor);
	//gabor_ker.Gabor_kernel(gray, gabor, 0, 2, 0.1, 11);
	//gabor.convertTo(gabor, CV_8U);
	////imwrite("gabor_1st.jpg", gabor);
	////Mat gaborT = gabor.clone();
	////gabor_ker.Gabor_kernel(gabor, gabor, 11, PI/3);// PI / 4);gabor
	////gabor_ker.Gabor_kernel(gabor, gabor, 11, 0);// PI / 4);gabor
	//imshow("gabor", gabor);

	//ContourDetecion contDetec(gray);
	//Mat dirImg(gray.size(),gray.type());
	////Mat gaborFilt(gray.size(), CV_32F);
	//Mat tmp(gray.size(), CV_32F);
	//tmp.copyTo(gaborFilt);
	//tmp.copyTo(suppMat);
	//tmp.copyTo(suppImg);
	//contDetec.GaborDirection(gray, gaborFilt, dirImg,3,0.3,11);
	////normalize(dirImg, dirImg, 0, 1, NORM_MINMAX);
	////imshow("dirImag", dirImg);
	////Mat suppMat(gray.size(), CV_32F);
	//contDetec.GaborSuppression(gaborFilt, dirImg, suppMat, 31);
	////Mat suppImg(gray.size(), CV_32F);
	//double beta = 0.5;
	//suppImg = gaborFilt - beta * suppMat;
	//normalize(suppMat, suppMat, 0, 1, NORM_MINMAX);
	//imshow("suppImg", suppImg);
	//imshow("gaborFilter", gaborFilt);
	//createTrackbar("suppressin cof", "suppImg", &betaSup100, 1000, changeBeta);

	//====Canny
	Mat can;
	faciMat.convertTo(can, CV_8U);
	Canny(can, can, 50, 100);
	imshow("can", can);

	waitKey(0);
}

//考虑像素值与局部均值的差作为算子的相应位置的值
void mySobel(Mat& src, Mat& dst, Mat& dstx, Mat& dsty)
{
	Mat SobKerx = (Mat_<int>(3, 3) << 1, 0, -1, -2, 0, -2, 1, 0, -1);
	Mat SobKery = (Mat_<int>(3, 3) << 1, 2, 1, 0, 0, 0, -1, -2, -1);
	std::cout << SobKerx << endl;
	std::cout << SobKery << endl;
	src.copyTo(dst);
	src.copyTo(dstx);
	src.copyTo(dsty);

	for (int i=1;i<src.rows-1;i++)
	{
		uchar* ptr = src.ptr<uchar>(i);
		uchar* ptrDst = dst.ptr<uchar>(i);
		uchar* ptrDstx = dstx.ptr<uchar>(i);
		uchar* ptrDsty = dsty.ptr<uchar>(i);

		for (int j=1;j<src.cols-1;j++)
		{
			//求均值
			int sum = 0;
			for (int m=-1;m<2;m++)
			{
				for (int n=-1;n<2;n++)
				{
					sum += (ptr + m)[j + n];
				}
			}
			int average = sum / 9;
			//滤波
			int rx = 0;
			int ry = 0;
			rx = SobKerx.at<int>(0, 0)*((ptr - 1)[j - 1] - average) + SobKerx.at<int>(0, 1)*((ptr - 1)[j] - average) + SobKerx.at<int>(0, 2)*((ptr - 1)[j + 1] - average)
				+ SobKerx.at<int>(1, 0)*((ptr)[j - 1] - average) + SobKerx.at<int>(1, 1)*((ptr)[j + 1] - average) + SobKerx.at<int>(1, 2)*((ptr)[j + 1] - average)
				+ SobKerx.at<int>(2, 0)*((ptr + 1)[j - 1] - average) + SobKerx.at<int>(2, 1)*((ptr + 1)[j] - average) + SobKerx.at<int>(2, 2)*((ptr + 1)[j + 1] - average);
			ry = SobKery.at<int>(0, 0)*((ptr - 1)[j - 1] - average) + SobKery.at<int>(0, 1)*((ptr - 1)[j] - average) + SobKery.at<int>(0, 2)*((ptr - 1)[j + 1] - average)
				+ SobKery.at<int>(1, 0)*((ptr)[j - 1] - average) + SobKery.at<int>(1, 1)*((ptr)[j] - average) + SobKery.at<int>(1, 2)*((ptr)[j + 1] - average)
				+ SobKery.at<int>(2, 0)*((ptr + 1)[j - 1] - average) + SobKery.at<int>(2, 1)*((ptr + 1)[j] - average) + SobKery.at<int>(2, 2)*((ptr + 1)[j + 1] - average);
			ptrDst[j] = sqrt(rx*rx + ry*ry);
			ptrDstx[j] = abs(rx);
			ptrDsty[j] = abs(ry);
		}
	}

}

//非极大值抑制和滞后阈值连接
void NMS(Mat& srcx, Mat& srcy, Mat& dst, int tl, int th)
{
	const double pi = 3.14159;
	Mat mag(srcx.size(), CV_32F);// gray.type());
	Mat pha(srcx.size(), CV_32F);// gray.type());
	std:cout << srcx.size() << ", " << srcy.size() << endl;

	srcx.convertTo(srcx, CV_32F);
	srcy.convertTo(srcy, CV_32F);

	cartToPolar(srcx, srcy, mag, pha);

	dst = mag.clone();
	for (int r = 1; r < srcx.rows - 1; r++)
	{
		for (int c = 1; c < srcx.cols - 1; c++)
		{
			float angle = pha.at<float>(r, c);
			if ((angle >= 1.875*pi || angle < 0.125*pi) || (angle >= 0.875 * pi && angle < 1.125*pi))//第一个条件要用||,处在圆首尾交接处
			{
				uchar v0 = mag.at<float>(r, c);
				uchar v1 = mag.at<float>(r, c - 1);
				uchar v2 = mag.at<float>(r, c + 1);
				if (v0 < v1 || v0 < v2)
				{
					dst.at<float>(r, c) = 0;
				}
			}
			if ((angle >= 0.125*pi && angle < 0.375*pi) || (angle >= 1.125 * pi && angle < 1.375*pi))
			{
				uchar v0 = mag.at<float>(r, c);
				uchar v1 = mag.at<float>(r + 1, c - 1);
				uchar v2 = mag.at<float>(r - 1, c + 1);
				if (v0 < v1 || v0 < v2)
				{
					dst.at<float>(r, c) = 0;
				}
			}
			if ((angle >= 0.375*pi && angle < 0.625*pi) || (angle >= 1.375 * pi && angle < 1.625*pi))
			{
				uchar v0 = mag.at<float>(r, c);
				uchar v1 = mag.at<float>(r - 1, c);
				uchar v2 = mag.at<float>(r + 1, c);
				if (v0 < v1 || v0 < v2)
				{
					dst.at<float>(r, c) = 0;
				}
			}
			if ((angle >= 0.625*pi && angle < 0.875*pi) || (angle >= 1.625 * pi && angle < 1.875*pi))
			{
				uchar v0 = mag.at<float>(r, c);
				uchar v1 = mag.at<float>(r - 1, c + 1);
				uchar v2 = mag.at<float>(r + 1, c - 1);
				if (v0 < v1 || v0 < v2)
				{
					dst.at<float>(r, c) = 0;
				}
			}
		}
	}

	dst.convertTo(dst, CV_8U);

	//threshold(dst0, dst0, th - 1, 255, THRESH_);
	//threshold(dst, dst, tl - 1, 255, THRESH_TOZERO);

	for (int r = 0; r < dst.rows; r++)
	{
		for (int c = 0; c < dst.cols; c++)
		{
			if (dst.at<uchar>(r, c) >= th)
			{
				dst.at<uchar>(r, c) = 255;
			}
			if (dst.at<uchar>(r, c) < tl)
			{
				dst.at<uchar>(r, c) = 0;
			}
		}
	}
	//threshold(dst, dst, 20, 255, THRESH_BINARY);
	imshow("dst0-1", dst);
	Mat visitedMat(dst.size(), dst.type(), Scalar(0));

	for (int i = 0; i < dst.rows; i++)
	{
		for (int j = 0; j < dst.cols; j++)
		{
			if (dst.at<uchar>(i, j) >= th && visitedMat.at<uchar>(i, j) == 0)
			{

				int flag = 0;
				DFS(Point(i, j), dst, visitedMat, th, &flag);
			}
		}
	}
	threshold(dst, dst, th - 1, 255, THRESH_BINARY);
	//imshow("dst1", dst);
}

//深度优先搜索，用于滞后阈值连接
void DFS(const Point point, Mat& img, Mat& visited, int th, int* flag)
{
	if (point.x<0 || point.x>img.rows - 1 ||
		point.y<0 || point.y>img.cols - 1)
		return;
	if (visited.at<uchar>(point.x, point.y) == 1 || img.at<uchar>(point.x, point.y) == 0)
	{
		return;
	}
	visited.at<uchar>(point.x, point.y) = 1;
	img.at<uchar>(point.x, point.y) = 255;
	for (int i = -1; i < 2; i++)
	{
		for (int j = -1; j < 2; j++)
		{
			if (i != 0 || j != 0)
			{
				DFS(Point(point.x + i, point.y + j), img, visited, th, flag);
			}
		}
	}

}

//高斯差分
void LOG(Mat& src, Mat& dst, double sigma1, double sigma2)
{
	Mat tmp1;
	Mat tmp2;
	//src.convertTo(src, CV_8U);
	GaussianBlur(src, tmp1, Size(5,5), sigma1, sigma1);
	GaussianBlur(src, tmp2, Size(5, 5), sigma2, sigma2);
	imshow("tmp1", tmp1);
	imshow("tmp2", tmp2);

	//dst = tmp1 - tmp2;
	for (int i=0;i<src.rows;i++)
	{
		uchar* ptr1 = tmp1.ptr<uchar>(i);
		uchar* ptr2 = tmp2.ptr<uchar>(i);
		uchar* ptrdst = dst.ptr<uchar>(i);
		for (int j=0;j<src.cols;j++)
		{
			ptrdst[j] = abs(ptr1[j] - ptr2[j]);
		}
	}
}

//滑动条回掉函数
void on_Change()
{
	Mat gabor(gray.size(), CV_32FC1);
	Gabor gabor_ker;
	double theta1 = g_theta / 180.0 * PI;
	if (g_lamda % 2 == 0)
		g_lamda += 1;
	double lamda1 = g_lamda;
	double gama1 = g_gama10 / 10.0;
	if (g_size % 2 == 0)
		g_size += 1;
	int size1 = g_size;
	double phi1 = g_phi / 180.0 * PI;
	double sr = g_sigmaRatio100 / 100.0;
	gabor_ker.Gabor_kernel(gray, gabor, theta1, lamda1, gama1, size1, phi1, sr);
	gabor.convertTo(gabor, CV_8U);
	imshow("Gabor Result", gabor);
}
void Theta(int, void*)
{
	on_Change();
}
void Lamda(int, void*)
{
	on_Change();
}
void Gama(int, void*)
{
	on_Change();
}
void Size_gabor(int, void*)
{
	on_Change();
}
void Phi(int, void*)
{
	on_Change();
}
void SR(int, void*)
{
	on_Change();
}

void changeBeta(int, void*)
{
	suppImg = gaborFilt - betaSup100/100.0 * suppMat;
	imshow("suppImg", suppImg);

}
