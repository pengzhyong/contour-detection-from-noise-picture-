#include<opencv2/opencv.hpp>
#include<vector>
#include "contours.h"
using namespace cv;
using namespace std;
int main_quzao()
{
	Mat src = imread("pic/kj.jpg");
	pyrDown(src, src);
	pyrDown(src, src);
	//imshow("src", src);
	Mat gray;
	cvtColor(src, gray, COLOR_BGR2GRAY);
	imshow("gray", gray);

	//threshold(gray, gray, 100, 255, THRESH_BINARY_INV);
	GaussianBlur(gray, gray, Size(3, 3), 1);

	//Mat temp;
	//GaussianBlur(gray, temp, Size(3, 3), 1);
	//gray = gray - temp;
	//threshold(gray, gray, 1, 255, THRESH_BINARY);

	//imshow("gray_edge", gray);
	//equalizeHist(gray, gray);

	//adaptiveThreshold(gray, gray, 255, ADAPTIVE_THRESH_MEAN_C, THRESH_BINARY_INV, 5, 5);
	Mat can;
	vector<vector<Point> > contour;
	//vector<vector<Vec4i> > hierarchy;
	Canny(gray, can, 30, 100);
	//imshow("can", can);
	Mat temp1;
	temp1 = can.clone();
	findContours(temp1, contour, RETR_CCOMP, CHAIN_APPROX_SIMPLE);
	Mat hiera(can.size(), can.type(), Scalar(0));
	drawContours(hiera, contour, -1, Scalar(255), 1);
	imshow("contours", hiera);
	vector<vector<Point> > new_contours;


	contours cont(contour, 50, 0.12);
	cont.cutShortContours();
	Mat new_hiera(can.size(), can.type(), Scalar(0));
	drawContours(new_hiera, cont.conts, -1, Scalar(255), 1);
	imshow("newcontours", new_hiera);
	//cont.sortContours();
	//cout << cont.conts[10] << endl;
	cont.filterContours();
	vector<vector<Point> > cc;
	cc.push_back(cont.conts[1]);
	Mat new_hiera1(can.size(), can.type(), Scalar(0));
	drawContours(new_hiera1, cont.conts, -1, Scalar(255), 1);
	imshow("newcontours1", new_hiera1);

	//int lenThresh = 50;
	//for (int i=1;i<contours.size();i++)
	//{
	//	//std:cout << contours[i].size() << endl;
	//	if (contours[i].size()>lenThresh)
	//	{
	//		new_contours.push_back(contours[i]);
	//	}
	//}
	//vector<vector<Point> > kk;
	//kk.push_back(new_contours[1]);

	//Mat new_hiera(can.size(), can.type(), Scalar(0));
	//drawContours(new_hiera, kk, -1, Scalar(255), 1);
	////imshow("newcontours", new_hiera);
	//cout << new_contours[1] << endl;
	//circle(new_hiera, Point(new_contours[1][0].x, new_contours[1][0].y), 2, Scalar(255), -1);
	//circle(new_hiera, Point(new_contours[1][1].x, new_contours[1][1].y), 2, Scalar(255), -1);
	//circle(new_hiera, Point(new_contours[1][2].x, new_contours[1][2].y), 3, Scalar(255), -1);
	//circle(new_hiera, Point(new_contours[1][3].x, new_contours[1][2].y), 3, Scalar(255), -1);
	//imshow("newcontours", new_hiera);


	waitKey(0);
	return 1;
}
void sortContours(vector<vector<Point> >& contours, vector<vector<Point> >& new_contours)
{

}