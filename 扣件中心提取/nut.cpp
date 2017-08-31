#include<opencv2/opencv.hpp>
#include<vector>
using namespace cv;
using namespace std;
int main()
{
	Mat src = imread("pic/lm.jpg");
	cvtColor(src, src, COLOR_BGR2GRAY);
	//GaussianBlur(src, src, Size(3, 3),3);
	imshow("src", src);
	Mat sob;
	Mat lap;
	//Sobel(src, sob, src.depth(), 1, 0);
	Laplacian(src, lap, src.depth(),3);
	normalize(lap, lap, 0, 255, NORM_MINMAX);

	convertScaleAbs(lap, lap);
	imshow("sob", lap);
	src = lap;
	vector<Point2f> coners;
	int maxConers = 100;
	goodFeaturesToTrack(src, coners, maxConers, 0.01, 5);
	for (int i=0;i<coners.size();i++)
	{
		circle(src, coners[i], 2, Scalar(255), -1, 8, 0);
	}
	//Mat blurImg;
	////GaussianBlur(src, blurImg, Size(3,3), 1);
	////src = src - blurImg;
	////Canny(src, src, 100, 200);
	//Mat dst;
	//Mat dstNorm;
	//dst = Mat::zeros(src.size(), CV_32FC1);
	//cornerHarris(src, dst, 3, 7, 0.01,BORDER_DEFAULT);
	//normalize(dst, dstNorm, 0, 255, NORM_MINMAX,CV_32FC1,Mat());
	//convertScaleAbs(dstNorm, dstNorm);
	//for (int j = 0; j < dst.rows; j++)
	//{
	//	for (int i = 0; i < dst.cols; i++)
	//	{
	//		if ((int)dstNorm.at<uchar>(j, i) > 30)
	//		{
	//			circle(src, Point(i, j), 5, Scalar(0), 1, 8, 0);
	//			//circle(src, Point(i, j), 5, Scalar(255, 0, 0), -1, 8, 0);
	//		}
	//	}
	//}
	//imshow("dst", dstNorm);
	//imshow("dst0", src);
	imshow("dst", src);
	waitKey(0);
	return 1;
}