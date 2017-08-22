#include<opencv2/opencv.hpp>
#include<vector>
#define WINDOW_NAME "DoG for edge detect"
using namespace cv;
using namespace std;
void Dog(Mat& src, Mat& dst);
void Sigma1(int, void*);
void Sigma2(int, void*);

Mat g_gray;
Mat g_edge;
double g_sigma1 = 1;
double g_sigma2 = 1;
int g_s1 = 10;
int g_s2 = 20;
int main_dog()
{
	Mat src = imread("pic/c10.jpg");
	//pyrDown(src, src);
	//pyrDown(src, src);
	imshow("src", src);
	cvtColor(src, g_gray, COLOR_BGR2GRAY);
	namedWindow(WINDOW_NAME, 1);
	createTrackbar("sigma1", WINDOW_NAME, &g_s1, 10, Sigma1);
	createTrackbar("sigma2", WINDOW_NAME, &g_s2, 10, Sigma2);
	Sigma1(g_s1,0);
	Sigma1(g_s2, 0);

	//Dog(g_gray, g_edge);
	//imshow("edge", g_edge);
	waitKey(0);
	return 1;
}
void Dog(Mat& src, Mat& dst)
{
	Mat img = src;
	Mat blur1, blur2;
	GaussianBlur(img, blur1, Size(3, 3),g_sigma1);
	//imshow("sigma1", blur1);

	GaussianBlur(img, blur2, Size(3, 3), g_sigma2);
	//imshow("sigma2", blur2);

	dst = blur1 - blur2;
	threshold(dst, dst, 1, 255, CV_THRESH_BINARY);
	imshow(WINDOW_NAME, dst);
}

void Sigma1(int, void*)
{
	g_sigma1 =  g_s1 / 10.0;
	Dog(g_gray, g_edge);
}
void Sigma2(int, void*)
{
	g_sigma2 = g_s2 /10.0;
	Dog(g_gray, g_edge);
}
