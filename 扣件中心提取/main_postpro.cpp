#include "PostProgress.h"
#include <vector>
#include <opencv2/opencv.hpp>
using namespace cv;
using namespace std;
void main()
{
	//Mat src = imread("pic/gabor2.jpg");
	Mat src = imread("pic/kj2.jpg");

	//pyrDown(src, src);
	//pyrDown(src, src);
	imshow("src", src);

	Mat gray;
	cvtColor(src, gray, COLOR_BGR2GRAY);
	GaussianBlur(gray, gray, Size(3, 3), 1);

	PostProgress postPro;
	Mat dst;
	Mat dstx;
	Mat dsty;

	postPro.mySobel(gray, dst, dstx, dsty);
	//imshow("afterSobel", dst);
	//imshow("afterSobel_x", dstx);
	//imshow("afterSobel_y", dsty);
	//waitKey(0);
	Mat nms;
	postPro.NMS(dstx, dsty, nms);
	//imshow("contour",nms);
	//waitKey(0);


	Mat can;
	vector<vector<Point> > contour;
	//vector<vector<Vec4i> > hierarchy;
	Canny(gray, can, 50, 100);
	Mat ker = (Mat_<uchar>(3, 3) << 1, 1, 1, 1, 1, 1, 1, 1, 1);
	std::cout << ker << std::endl;
	//morphologyEx(can, can, MORPH_CLOSE, ker);
	//imshow("can", can);
	Mat temp1;
	temp1 = can.clone();
	findContours(temp1, contour, RETR_CCOMP, CHAIN_APPROX_NONE);
	Mat hiera(can.size(), can.type(), Scalar(0));
	drawContours(hiera, contour, -1, Scalar(255), 1);

	imshow("beforeSort_contours", hiera);
	imwrite("before sorted contours.jpg", hiera);
	vector<vector<Point> > outCont;

	postPro.SortContours(contour, outCont);
	Mat hiera1(can.size(), can.type(), Scalar(0));
	drawContours(hiera1, contour, -1, Scalar(255), 1);
	imshow("after sorted contours", hiera1);
	waitKey(0);

	//对轮廓进行多边形逼近，去噪
	
	vector<vector<Point> > newContour;
	double curRadio = 10;//总点数/线段数
	for (auto iter1 = contour.begin();iter1!=contour.end();iter1++)
	{
		double xdata[3000] = { 0 };
		double ydata[3000] = { 0 };
		int pNums = 0;
		for (auto iter2 = iter1->begin(); iter2 < iter1->end(); iter2++)
		{
			xdata[pNums] = iter2->x;
			ydata[pNums] = iter2->y;
			pNums++;
		}
		int resArr[300] = { 0 };
		int ind = 0;
		double th = 10;//直线拟合距离阈值
		postPro.RamerFunc(xdata, ydata, 0, pNums - 1, resArr, ind, th);
		if (pNums/ind > curRadio && pNums > 30)
		{
			newContour.push_back(*iter1);
		}
	}
	Mat newHiera(can.size(), can.type(), Scalar(0));

	drawContours(newHiera, newContour, -1, Scalar(255), 1);
	imshow("newContours", newHiera);
	waitKey(0);
	//cv::RetrievalModes;
}