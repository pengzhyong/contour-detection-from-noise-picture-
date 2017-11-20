#include "PostProgress.h"
#include <opencv2/opencv.hpp>
using namespace cv;
void main()
{
	Mat src = imread("pic/lm.jpg");

	//pyrDown(src, src);
	//pyrDown(src, src);
	//imshow("src", src);

	Mat gray;
	cvtColor(src, gray, COLOR_BGR2GRAY);
	GaussianBlur(gray, gray, Size(3, 3), 1);

	Mat can;
	Canny(gray, can, 30, 80);
	imshow("can", can);
	//waitKey(0); 

	//测试图
	//{
	//	Mat testImg(gray.size(), gray.type(), Scalar(0));
	//	line(testImg, Point(50, 30), Point(50, 150), Scalar(255));
	//	line(testImg, Point(50, 150), Point(60, 100), Scalar(255));
	//	line(testImg, Point(10, 70), Point(70, 120), Scalar(255));
	//	line(testImg, Point(10, 10), Point(11, 10), Scalar(255));
	//	line(testImg, Point(11, 10), Point(11, 11), Scalar(255));
	//	line(testImg, Point(11, 11), Point(12, 11), Scalar(255));
	//	//line(testImg, Point(15, 10), Point(15, 11), Scalar(255));
	//	//line(testImg, Point(15, 11), Point(16, 18), Scalar(255));
	//	circle(testImg, Point(50, 70), 20, Scalar(255));
	//	imshow("testImg", testImg);
	//}

	
	//多次迭代去毛刺
	PostProgress postPro;
	vector<Point> vEndPoints;
	postPro.FindEndPoint(can, vEndPoints);
	vector<Point> vCrossPoints;
	for(int t=0;t<20;t++)
	{
		std::cout << "t= " << t << std::endl;
		vEndPoints.clear();
		vCrossPoints.clear();
		postPro.FindEndPoint(can, vEndPoints);
		postPro.FindCrossPoint(can, vCrossPoints);
		postPro.RemoveBurr(can, vEndPoints, vCrossPoints, 10 + 0.5 * t);
		imshow("testImg_while", can);
		waitKey(100);
	}

	//标记端点和交叉点
	{
		if(can.channels()==1)
			cvtColor(can, can, COLOR_GRAY2BGR);
		for (int i = 0; i < vCrossPoints.size(); i++)
		{
			//circle(can, vCrossPoints[i], 3, Scalar(0,0,255),-1);
			can.at<Vec3b>(vCrossPoints[i].y, vCrossPoints[i].x)[0] = 0;
			can.at<Vec3b>(vCrossPoints[i].y, vCrossPoints[i].x)[1] = 0;
			can.at<Vec3b>(vCrossPoints[i].y, vCrossPoints[i].x)[2] = 255;
		}

		for (int i = 0; i < vEndPoints.size(); i++)
		{
			//circle(can, vEndPoints[i], 3, Scalar(255,0,0));
			//circle(can, vEndPoints[i], 1, Scalar(255));
			can.at<Vec3b>(vEndPoints[i].y, vEndPoints[i].x)[0] = 255;

			can.at<Vec3b>(vEndPoints[i].y, vEndPoints[i].x)[1] = 0;
			can.at<Vec3b>(vEndPoints[i].y, vEndPoints[i].x)[2] = 0;
		}
		namedWindow("End&Cross", 0);
		imshow("End&Cross", can);
		waitKey(0);
		if (can.channels() == 3)
			cvtColor(can, can, COLOR_BGR2GRAY);
	}
}