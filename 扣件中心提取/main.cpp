#include<opencv2/opencv.hpp>
#include<vector>
using namespace cv;
using namespace std;
int main()
{
	Mat src = imread("pic/c11.jpg");
	pyrDown(src, src);
	//pyrDown(src, src);
	imshow("src", src);
	Mat gray;
	cvtColor(src, gray, COLOR_BGR2GRAY);
	//threshold(gray, gray, 100, 255, THRESH_BINARY_INV);
	//GaussianBlur(gray, gray, Size(3, 3), 1);

	//equalizeHist(gray, gray);

	//adaptiveThreshold(gray, gray, 255, ADAPTIVE_THRESH_MEAN_C, THRESH_BINARY_INV, 5, 5);
	Mat can;
	vector<vector<Point> > contours;
	//vector<vector<Vec4i> > hierarchy;
	Canny(gray, can, 50, 100);
	
	imshow("canny", can);
	Mat temp;
	temp = can.clone();
	findContours(temp,contours, RETR_CCOMP, CHAIN_APPROX_SIMPLE);
	Mat hiera(can.size(), can.type(),Scalar(0));
	drawContours(hiera, contours, -1, Scalar(255), 1);
	imshow("contours", hiera);
	vector<vector<Point> > new_contours;
	int longThresh = 50;
	double radio = 1.1;
	for (int i=0;i<contours.size();i++)
	{
		int counts = 1;
		
		if (contours[i].size() < longThresh)
			continue;
		cout << "contours.size: " << contours[i].size() << endl;
		for (int j=1;j<contours[i].size()-2;j++)
		{
			Point orient;
			Point prevorient;
			Point nextorient;
			orient.x = contours[i][j + 1].x - contours[i][j].x;
			orient.y = contours[i][j + 1].y - contours[i][j].y;
			prevorient.x = contours[i][j].x - contours[i][j - 1].x;
			prevorient.y = contours[i][j].y - contours[i][j - 1].y;
			nextorient.x = contours[i][j+2].x - contours[i][j+1].x;
			nextorient.y = contours[i][j+2].y - contours[i][j+1].y;
			if (((orient.x==prevorient.x && orient.y==prevorient.y)||
				(orient.x==prevorient.x && abs(orient.y-prevorient.y)==1)||
				(orient.y == prevorient.y && abs(orient.x - prevorient.x) == 1)) &&
				((orient.x == nextorient.x && orient.y == nextorient.y) ||
				(orient.x == nextorient.x && abs(orient.y - nextorient.y) == 1) ||
					(orient.y == nextorient.y && abs(orient.x - nextorient.x) == 1)))
			{
				continue;
			}
			counts++;//轮廓上拐点数
		}
		double AreaRadio = contourArea(contours[i])*1.0 / contours[i].size();
		cout << "AreaRadio:" << contourArea(contours[i])*1.0 / contours[i].size() << endl;
		cout << "radio: " << contours[i].size()*1.0 / counts << endl;
		if ((contours[i].size()*1.0 / counts) > radio && AreaRadio < 0.4)//轮廓总长度与拐点数的比值
		{
			new_contours.push_back(contours[i]);
			
		}
	}
	cout << "new_contours size: " << new_contours.size() << endl;
	Mat hiera1(can.size(), can.type(), Scalar(0));
	drawContours(hiera1, new_contours, -1, Scalar(255), 1);
	imshow("newcontours", hiera1);

	can = hiera1;

	//imshow("thresh", gray);
	//======去除密度大的斑块
	//int sumThresh = 255 * 12;
	//for (int i=5;i<gray.rows-5;i++)
	//{
	//	uchar* ptr = can.ptr(i);
	//	for (int j=5;j<gray.cols-5;j++)
	//	{
	//		int sum = 0;
	//		for (int k=-5;k<6;k++)
	//		{
	//			for (int m=-5;m<6;m++)
	//			{
	//				sum += *(ptr + k * gray.step + (j + m)*gray.channels());
	//			}
	//			
	//		}
	//		if (sum>sumThresh)
	//		{
	//			*(ptr + j) = 0;
	//		}
	//	}
	//}
	//imshow("sumThresh", can);
	//imwrite("prof4.jpg",can);
	vector<Vec3f> circles;
	HoughCircles(can, circles, HOUGH_GRADIENT, 2, 1, 100, 30, 55, 70);
	//cout << "cirlces[0]: " << circles[0] << endl;

	for (size_t i = 0; i < 1 && i<circles.size(); i++)
	{
		Point center(cvRound(circles[i][0]), cvRound(circles[i][1]));
		int radius = cvRound(circles[i][2]);
		// draw the circle center
		//circle(gray, center, 3, Scalar(255,255,255), -1, 8, 0);
		// draw the circle outline 
		circle(gray, center, radius, Scalar(255,255,255), 1, 8, 0);
		std::cout << "radius: " << radius << endl;
	}
	if (circles.size()<1)
	{
		cout << "unable to detect the big circle " << endl;
		return 0;
	}
	int bigCricleLfx = circles[0][0] - (int)circles[0][2];
	int bigCricleLfy = circles[0][1] - (int)circles[0][2];

	Mat can1(can, Rect(circles[0][0]-(int)circles[0][2], circles[0][1]-(int)circles[0][2],2* (int)circles[0][2],2*(int)circles[0][2]));
	vector<Vec3f> circles1;

	HoughCircles(can1, circles1, HOUGH_GRADIENT, 2, 1, 100, 40, 20, 55);
	//cout << "cirlces[1]: " << circles1[0] << endl;

	for (size_t i = 0; i < 1 &&i < circles1.size(); i++)
	{
		Point center(cvRound(circles1[i][0])+bigCricleLfx, cvRound(circles1[i][1])+bigCricleLfy);
		int radius = cvRound(circles1[i][2]);
		// draw the circle center
		circle(gray, center, 3, Scalar(255,255,255), -1, 8, 0);
		// draw the circle outline 
		circle(gray, center, radius, Scalar(200,200,200), 1, 8, 0);
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


//#include<iostream>  
//#include<opencv2/opencv.hpp>  
//#include<vector>  
//
//using namespace std;
//using namespace cv;
//
//int g_CannyThred = 180, g_CannyP = 0, g_CannySize = 0, g_HoughThred = 15, g_HoughThick = 9;
//int g_Blue = 255, g_Green = 255, g_Red = 0;
//int g_nWay = 0;
//int g_nHoughLineMax = 10, g_nHoughLineMin = 50;
//
//int g_nDp = 0;
//int g_nMinDist = 5;
//int g_nMinRadius = 0, g_nMaxRadius = 0;
//
//int main()
//{
//	//此路径为系统盘里的图片位置路径
//	Mat srcImage = imread("pic/c9.jpg");
//
//	/*Mat srcImage = imread("C:\\Users\\Yao%20Yao\\AppData\\Local\\Temp\\1.jpg");*/
//	imshow("【原图】", srcImage);
//
//	Mat grayImage;
//	cvtColor(srcImage, grayImage, CV_BGR2GRAY);
//	GaussianBlur(grayImage, grayImage, Size(9, 9), 2, 2);
//
//	Mat cannyImage;
//	vector<Vec3f> circles;
//
//	namedWindow("【滚动条窗口】", 0);
//	createTrackbar("dp", "【滚动条窗口】", &g_nDp, 100, 0);
//	createTrackbar("MinDist", "【滚动条窗口】", &g_nMinDist, 100, 0);
//	createTrackbar("CannyThred", "【滚动条窗口】", &g_CannyThred, 300, 0);
//	createTrackbar("HoughThred", "【滚动条窗口】", &g_HoughThred, 255, 0);
//	createTrackbar("Blue", "【滚动条窗口】", &g_Blue, 255, 0);
//	createTrackbar("Green", "【滚动条窗口】", &g_Green, 255, 0);
//	createTrackbar("Red", "【滚动条窗口】", &g_Red, 255, 0);
//	createTrackbar("Bgr/Gray", "【滚动条窗口】", &g_nWay, 1, 0);
//	createTrackbar("HoughThred", "【滚动条窗口】", &g_HoughThred, 200, 0);
//	createTrackbar("MinRadius", "【滚动条窗口】", &g_nMinRadius, 100, 0);
//	createTrackbar("MaxRadius", "【滚动条窗口】", &g_nMaxRadius, 100, 0);
//	createTrackbar("HoughThick", "【滚动条窗口】", &g_HoughThick, 100, 0);
//
//	char key;
//	Mat dstImage;
//	while (1)
//	{
//		HoughCircles(grayImage, circles, CV_HOUGH_GRADIENT, (double)((g_nDp + 1.5)), (double)g_nMinDist + 1
//			, (double)g_CannyThred + 1, g_HoughThred + 1);// , g_nMinRadius, g_nMaxRadius);
//
//		/*HoughCircles(grayImage, circles, CV_HOUGH_GRADIENT, 1.5, 10, 200, 100, 0, 0);*/
//
//		//显示线段  
//		for (size_t i = 0; i < circles.size(); i++)
//		{
//			if (g_nWay)
//				grayImage.copyTo(dstImage);
//			else
//				srcImage.copyTo(dstImage);
//
//			circle(dstImage, Point(cvRound(circles[i][0]), cvRound(circles[i][1])), cvRound(circles[i][2])
//				, Scalar(g_Blue, g_Green, g_Red), g_HoughThick);
//		}
//
//		imshow("【处理后】", dstImage);
//
//		key = waitKey(1);
//		if (key == 27)
//			break;
//	}
//
//	return 0;
//}