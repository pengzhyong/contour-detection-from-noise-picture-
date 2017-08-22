#pragma once
#include<opencv2/opencv.hpp>
#include<vector>
using namespace cv;
using namespace std;
class contours
{
public:
	contours() {};
	contours(vector<vector<Point> >& c, int cutLen, double radioT);
	void cutShortContours();
	void sortContours();
	void filterContours();
	vector<vector<Point> > conts;
private:

	int lenThresh = 100;
	double radioThresh = 0.2;
	Point findEndPoints(vector<Point>&);
	bool isNeigbour(Point&, Point&);
	int distance(Point&, Point&);
};
