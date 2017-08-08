#include<opencv2/opencv.hpp>
#include<vector>
using namespace cv;
using namespace std;
void DoG(Mat& img1, Mat& dst, double sigma1, double step = 0.1 );
int _main()
{

	//Mat src = imread("pic/c10.jpg");
	//pyrDown(src, src);
	////pyrDown(src, src);
	//imshow("src", src);
	//Mat gray;
	//cvtColor(src, gray, COLOR_BGR2GRAY);
	////threshold(gray, gray, 100, 255, THRESH_BINARY_INV);
	////GaussianBlur(gray, gray, Size(3, 3), 1);

	////equalizeHist(gray, gray);

	////adaptiveThreshold(gray, gray, 255, ADAPTIVE_THRESH_MEAN_C, THRESH_BINARY_INV, 5, 5);
	//Mat can;
	//vector<vector<Point> > contours;
	////vector<vector<Vec4i> > hierarchy;
	//Canny(gray, can, 50, 100);

	//imshow("canny", can);
	//Mat temp;
	//temp = can.clone();
	//findContours(temp, contours, RETR_CCOMP, CHAIN_APPROX_SIMPLE);
	//Mat hiera(can.size(), can.type(), Scalar(0));
	//drawContours(hiera, contours, -1, Scalar(255), 1);
	//imshow("contours", hiera);
	//vector<vector<Point> > new_contours;
	//int longThresh = 50;
	//double radio = 1.1;
	//for (int i = 0; i < contours.size(); i++)
	//{
	//	int counts = 1;

	//	if (contours[i].size() < longThresh)
	//		continue;
	//	cout << "contours.size: " << contours[i].size() << endl;
	//	for (int j = 1; j < contours[i].size() - 2; j++)
	//	{
	//		Point orient;
	//		Point prevorient;
	//		Point nextorient;
	//		orient.x = contours[i][j + 1].x - contours[i][j].x;
	//		orient.y = contours[i][j + 1].y - contours[i][j].y;
	//		prevorient.x = contours[i][j].x - contours[i][j - 1].x;
	//		prevorient.y = contours[i][j].y - contours[i][j - 1].y;
	//		nextorient.x = contours[i][j + 2].x - contours[i][j + 1].x;
	//		nextorient.y = contours[i][j + 2].y - contours[i][j + 1].y;
	//		if (((orient.x == prevorient.x && orient.y == prevorient.y) ||
	//			(orient.x == prevorient.x && abs(orient.y - prevorient.y) == 1) ||
	//			(orient.y == prevorient.y && abs(orient.x - prevorient.x) == 1)) &&
	//			((orient.x == nextorient.x && orient.y == nextorient.y) ||
	//			(orient.x == nextorient.x && abs(orient.y - nextorient.y) == 1) ||
	//				(orient.y == nextorient.y && abs(orient.x - nextorient.x) == 1)))
	//		{
	//			continue;
	//		}
	//		counts++;//轮廓上拐点数
	//	}
	//	double AreaRadio = contourArea(contours[i])*1.0 / contours[i].size();
	//	cout << "AreaRadio:" << contourArea(contours[i])*1.0 / contours[i].size() << endl;
	//	cout << "radio: " << contours[i].size()*1.0 / counts << endl;
	//	if ((contours[i].size()*1.0 / counts) > radio && AreaRadio < 0.4)//轮廓总长度与拐点数的比值
	//	{
	//		new_contours.push_back(contours[i]);

	//	}
	//}
	//cout << "new_contours size: " << new_contours.size() << endl;
	//Mat hiera1(can.size(), can.type(), Scalar(0));
	//drawContours(hiera1, new_contours, -1, Scalar(255), 1);
	//imshow("newcontours", hiera1);



	//Mat src = imread("pic/c9.jpg");
	////pyrDown(src, src);
	//pyrDown(src, src);
	//imshow("src", src);
	//Mat gray;
	//cvtColor(src, gray, COLOR_BGR2GRAY);
	//equalizeHist(gray, gray);
	//Mat can;
	//Canny(gray, can, 100, 150);
	//imshow("can", can);
	Mat can = imread("pic/c7.jpg");
	Mat templateImg = imread("pic/template.jpg");
	imshow("tempalte", templateImg);

	/*Mat templateImg(Size(150, 150), can.type(), Scalar(0));
	circle(templateImg, Point(75, 75), 62, Scalar(255));
	circle(templateImg, Point(75, 75), 27, Scalar(255));
	imshow("templateImage", templateImg);*/
	Mat results(Size(can.rows - templateImg.rows + 1, can.cols - templateImg.cols + 1), CV_32FC1);
	matchTemplate(can, templateImg, results, TM_SQDIFF);
	imshow("results", results);
	//cout << results << endl;
	double max = results.at<float>(0, 0);
	Point minLoc;
	for (int i=0;i<results.rows;i++)
	{
		for (int j = 0; j < results.cols;j++)
		{
			if (results.at<float>(i,j)<max)
			{
				max = results.at<float>(i, j);
				minLoc.x = i;
				minLoc.y = j;
			}
		}
	}
// 	int radius1 = 62;
// 	int radius2 = 27;
	cout << "minLoc: " << minLoc.x << " " << minLoc.y << endl;
	circle(can, Point(minLoc.y+templateImg.rows/2-5 , minLoc.x +templateImg.cols/2-27), 50, Scalar(255,255,255));
	//circle(can, Point(minLoc.x +templateImg.rows/2-5, minLoc.y +templateImg.cols/2-27), 27, Scalar(255,255,255));
	circle(templateImg, Point(templateImg.rows / 2-5,templateImg.cols / 2-27), 50, Scalar(255, 255, 255));
	imshow("templateImg", templateImg);

	imshow("final", can);




	//Mat img_1;
	//Mat img_2;
	//gray.copyTo(img_1);
	//gray.copyTo(img_2);
	//Mat dst;
	//GaussianBlur(img_1, img_1, Size(5, 5), 0.5);
	//GaussianBlur(img_2, img_2, Size(5, 5),0.8);
	//dst = img_1 -img_2;
	//threshold(dst, dst, 10, 255, THRESH_BINARY);
	////normalize(dst, dst, 0, 255);
	//GaussianBlur(gray, gray, Size(3, 3), 0.5);
	//imshow("dog", dst);
	/*Mat dog[5];
	dog[0] = gray;
	for (int i=1;i<5;i++)
	{
		DoG(dog[i-1], dog[i], 1, 1.0);

	}
	imshow("dog[i]", dog[4]);
	waitKey(0);*/


		//////threshold(gray, gray, 100, 255, THRESH_BINARY_INV);
	//////GaussianBlur(gray, gray, Size(3, 3), 1);
	//////equalizeHist(gray, gray);
	//////adaptiveThreshold(gray, gray, 255, ADAPTIVE_THRESH_MEAN_C, THRESH_BINARY_INV, 5, 5);
	////Mat can;
	////vector<vector<Point> > contours;
	//////vector<vector<Vec4i> > hierarchy;
	////Canny(gray, can, 50, 100);
	////imshow("canny", can);
	////Mat temp;
	////temp = can.clone();
	////findContours(temp, contours, RETR_CCOMP, CHAIN_APPROX_SIMPLE);
	////Mat hiera(can.size(), can.type(), Scalar(0));
	////drawContours(hiera, contours, -1, Scalar(255), 1);
	////imshow("contours", hiera);
	////vector<vector<Point> > new_contours;
	////int longThresh = 50;
	////double radio = 1.1;
	////for (int i = 0; i < contours.size(); i++)
	////{
	////	int counts = 1;
	////	if (contours[i].size() < longThresh)
	////		continue;
	////	cout << "contours.size: " << contours[i].size() << endl;
	////	for (int j = 1; j < contours[i].size() - 2; j++)
	////	{
	////		Point orient;
	////		Point prevorient;
	////		Point nextorient;
	////		orient.x = contours[i][j + 1].x - contours[i][j].x;
	////		orient.y = contours[i][j + 1].y - contours[i][j].y;
	////		prevorient.x = contours[i][j].x - contours[i][j - 1].x;
	////		prevorient.y = contours[i][j].y - contours[i][j - 1].y;
	////		nextorient.x = contours[i][j + 2].x - contours[i][j + 1].x;
	////		nextorient.y = contours[i][j + 2].y - contours[i][j + 1].y;
	////		if (((orient.x == prevorient.x && orient.y == prevorient.y) ||
	////			(orient.x == prevorient.x && abs(orient.y - prevorient.y) == 1) ||
	////			(orient.y == prevorient.y && abs(orient.x - prevorient.x) == 1)) &&
	////			((orient.x == nextorient.x && orient.y == nextorient.y) ||
	////			(orient.x == nextorient.x && abs(orient.y - nextorient.y) == 1) ||
	////				(orient.y == nextorient.y && abs(orient.x - nextorient.x) == 1)))
	////		{
	////			continue;
	////		}
	////		counts++;//轮廓上拐点数
	////	}
	////	double AreaRadio = contourArea(contours[i])*1.0 / contours[i].size();
	////	cout << "AreaRadio:" << contourArea(contours[i])*1.0 / contours[i].size() << endl;
	////	cout << "radio: " << contours[i].size()*1.0 / counts << endl;
	////	if ((contours[i].size()*1.0 / counts) > radio && AreaRadio < 0.4)//轮廓总长度与拐点数的比值
	////	{
	////		new_contours.push_back(contours[i]);
	////	}
	////}
	////cout << "new_contours size: " << new_contours.size() << endl;
	////Mat hiera1(can.size(), can.type(), Scalar(0));
	////drawContours(hiera1, new_contours, -1, Scalar(255), 1);
	////imshow("newcontours", hiera1);
	////can = hiera1;
	////vector<Vec3f> circles;
	////HoughCircles(can, circles, HOUGH_GRADIENT, 2, 1, 100, 30, 55, 70);
	//////cout << "cirlces[0]: " << circles[0] << endl;
	////for (size_t i = 0; i < 1 && i < circles.size(); i++)
	////{
	////	Point center(cvRound(circles[i][0]), cvRound(circles[i][1]));
	////	int radius = cvRound(circles[i][2]);
	////	// draw the circle center
	////	//circle(gray, center, 3, Scalar(255,255,255), -1, 8, 0);
	////	// draw the circle outline 
	////	circle(gray, center, radius, Scalar(255, 255, 255), 1, 8, 0);
	////	std::cout << "radius: " << radius << endl;
	////}
	////if (circles.size() < 1)
	////{
	////	cout << "unable to detect the big circle " << endl;
	////	return 0;
	////}
	////int bigCricleLfx = circles[0][0] - (int)circles[0][2];
	////int bigCricleLfy = circles[0][1] - (int)circles[0][2];
	////Mat can1(can, Rect(circles[0][0] - (int)circles[0][2], circles[0][1] - (int)circles[0][2], 2 * (int)circles[0][2], 2 * (int)circles[0][2]));
	////vector<Vec3f> circles1;
	////HoughCircles(can1, circles1, HOUGH_GRADIENT, 2, 1, 100, 40, 20, 55);
	//////cout << "cirlces[1]: " << circles1[0] << endl;
	////for (size_t i = 0; i < 1 && i < circles1.size(); i++)
	////{
	////	Point center(cvRound(circles1[i][0]) + bigCricleLfx, cvRound(circles1[i][1]) + bigCricleLfy);
	////	int radius = cvRound(circles1[i][2]);
	////	// draw the circle center
	////	circle(gray, center, 3, Scalar(255, 255, 255), -1, 8, 0);
	////	// draw the circle outline 
	////	circle(gray, center, radius, Scalar(200, 200, 200), 1, 8, 0);
	////	std::cout << "radius: " << radius << endl;
	////}
	////imshow("HoughCircle", gray);
	waitKey(0);
	return 1;
}
void DoG(Mat& img1, Mat& dst, double sigma1, double step)
{
	Mat img_1;
	Mat img_2;
	GaussianBlur(img1, img_1, Size(3, 3), sigma1);
	GaussianBlur(img1, img_2, Size(3, 3), sigma1 + step);
	dst = img_1 - img_2;
	threshold(dst, dst, 1, 255,THRESH_BINARY);
}