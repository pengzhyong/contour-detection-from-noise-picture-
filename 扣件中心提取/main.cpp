#include<opencv2/opencv.hpp>
#include<vector>
using namespace cv;
using namespace std;
int main()
{
	Mat src = imread("pic/c10.jpg");
	pyrDown(src, src);
	//pyrDown(src, src);
	imshow("src", src);
	Mat gray;
	cvtColor(src, gray, COLOR_BGR2GRAY);
	//threshold(gray, gray, 100, 255, THRESH_BINARY_INV);
	//GaussianBlur(gray, gray, Size(3, 3), 3);

	//equalizeHist(gray, gray);

	//adaptiveThreshold(gray, gray, 255, ADAPTIVE_THRESH_MEAN_C, THRESH_BINARY_INV, 5, 5);
	Mat can;
	Canny(gray, can, 50, 100);

	imshow("canny", can);
	//imshow("thresh", gray);
	int sumThresh = 255 * 12;
	for (int i=5;i<gray.rows-5;i++)
	{
		uchar* ptr = can.ptr(i);
		for (int j=5;j<gray.cols-5;j++)
		{
			int sum = 0;
			for (int k=-5;k<6;k++)
			{
				for (int m=-5;m<6;m++)
				{
					sum += *(ptr + k * gray.step + (j + m)*gray.channels());
				}
				
			}
			if (sum>sumThresh)
			{
				*(ptr + j) = 0;
			}
		}
	}
	imshow("sumThresh", can);
	vector<Vec3f> circles;
	HoughCircles(can, circles, HOUGH_GRADIENT, 2, 1,100,160);
	for (size_t i = 0; i < circles.size(); i++)
	{
		Point center(cvRound(circles[i][0]), cvRound(circles[i][1]));
		int radius = cvRound(circles[i][2]);
		// draw the circle center
		circle(gray, center, 3, Scalar(255,255,255), -1, 8, 0);
		// draw the circle outline 
		circle(gray, center, radius, Scalar(255,255,255), 3, 8, 0);
		std::cout << "radius: " << radius << endl;
	}
	HoughCircles(can, circles, HOUGH_GRADIENT, 2, 1, 100, 40, 20, 30);
	for (size_t i = 0; i < circles.size(); i++)
	{
		Point center(cvRound(circles[i][0]), cvRound(circles[i][1]));
		int radius = cvRound(circles[i][2]);
		// draw the circle center
		circle(gray, center, 3, Scalar(255,255,255), -1, 8, 0);
		// draw the circle outline 
		circle(gray, center, radius, Scalar(255,255,255), 3, 8, 0);
		std::cout << "radius: " << radius << endl;
	}
	
	imshow("HoughCircle", gray);
	//Mat close_kenerl = (Mat_<uchar>(3, 3) <<
	//	1, 1, 1,
	//	1, 1, 1,
	//	1, 1, 1);
	//morphologyEx(gray, gray, MORPH_ERODE, close_kenerl);
	//morphologyEx(gray, gray, MORPH_ERODE, close_kenerl);
	//morphologyEx(gray, gray, MORPH_CLOSE, close_kenerl);

	//imshow("MORPH", gray);

	//GaussianBlur(src, src, Size(3, 3), 5);
	//Mat cann;
	//Canny(src, cann, 100, 200);
	//imshow("cann", cann);
	//Mat can;
	//Mat close_kenerl = (Mat_<uchar>(5, 5) <<
	//	1, 1, 1, 1, 1,
	//	1, 1, 1, 1, 1,
	//	1, 1, 1, 1, 1,
	//	1, 1, 1, 1, 1,
	//	1, 1, 1, 1, 1);
	//morphologyEx(cann, can, MORPH_CLOSE, close_kenerl);
	//imshow("can_close", can);
	//imshow("cha", cann-(can - cann));

	//vector<vector<Point> > contours;
	//vector<Vec4i> hierarchy;
	//findContours(can, contours, hierarchy, CV_RETR_CCOMP, CHAIN_APPROX_NONE);
	//int num_cont = 0;
	//int max_len = 0;
	//for (int i = 0; i < contours.size(); i++)
	//{
	//	num_cont += contours[i].size();
	//	if (contours.size() > max_len)
	//	{
	//		max_len = contours.size();
	//	}
	//}
	//int avr_len = num_cont / contours.size();
	//vector<vector<Point> > new_contours;
	//for (int i = 0; i < contours.size(); i++)
	//{
	//	if (contours[i].size() > 5 * avr_len)
	//	{
	//		new_contours.push_back(contours[i]);
	//	}
	//}
	//Mat cont(can.size(), can.type(), Scalar(0));
	//Scalar sca(255, 255, 255);
	//for (int idx = 0; idx < new_contours.size(); idx++)//idx>=0;idx=hierarchy[idx][0])
	//{
	//	//drawContours(cont, new_contours, idx, sca);
	//	drawContours(cont, new_contours, idx, sca);

	//}
	//imshow("new_contours", cont);
	waitKey(0);
	return 1;
}