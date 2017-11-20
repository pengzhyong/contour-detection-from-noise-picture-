#pragma once
#include <opencv2/opencv.hpp>
#include <vector>
using namespace std;
using namespace cv;
class PostProgress
{
public:
	PostProgress();
	~PostProgress();
	//�Ǽ���ֵ����
	void NMS(Mat& srcx, Mat& srcy, Mat& dst);
	//Sobel����
	void MySobel(Mat& src, Mat& dst, Mat& dstx, Mat& dsty);
	//���������ж���αƽ���ͨ������������������Ŀ�ı�ֵ���ж������Ƿ����������������������ĵ���������
	void RamerFunc(double* xdata, double* ydata, int beginpoint, int endpoint, int* ResultArray, int &Index, double th);
	//����������������һ��Ӧ�ڶ���αƽ�֮ǰ���У��������Ϊ�������򣬿��ܻᵼ�±ƽ����ִ���Ľ��
	void SortContours(vector<vector<Point> >& inCont, vector<vector<Point> >& outCont);
	//�ҳ������еĽ����,���ؽ����������ţ�������ģʽ�У�ֻ��T��ģʽ�Ľṹ���ǽ���㣬����8��T��ģʽ
	void FindCross(const vector<Point >& inCnt, vector<int >& crossIdx);
	
	//================���º�������ͼ�������ֱ�ӽ���
	//�ҵ�ͼ���еĶ˵㣬Ϊȥë����׼��
	void FindEndPoint(Mat& srcImg, vector<Point>& vEndPoints);
	void FindCrossPoint(Mat& srcImg, vector<Point>& vCrossPoints);
	bool IsCrossPont(const Mat& srcImg, Point pt);
	void RemoveBurr(Mat& srcImg, vector<Point>& vEndPoint, vector<Point>& vCrossPoint, int lenThresh = 20);
private:
	Point findEndPoints(vector<Point>&, int& );
	bool isNeigbour(Point & a, Point & b, int isFourNeig);
	bool isNeigbour(const Point & a, const Point & b);
	int distance(Point&, Point&);
	void splitContours(vector<Point >& srcConts, vector<vector<Point > >& dstCnt);//�����ж��������ߵ��������зֽ⣬ʹ��ÿһ��������û�з�֧;���õݹ鷽ʽ
	void removeRepeat(vector<Point>&);
	//�ҳ�û�н����������������ӽ�����һ������㿪ʼ���٣�ֱ���ߵ��˵������һ�������
	void getSingleCnt(int beginPointIdx, vector<Point >& inCnt, int* isVisited, vector<Point >& outCnt);
};

