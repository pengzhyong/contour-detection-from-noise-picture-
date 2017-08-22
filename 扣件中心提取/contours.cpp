#include<opencv2/opencv.hpp>
#include<vector>
#include "contours.h"
using namespace cv;
using namespace std;
contours::contours(vector<vector<Point> >& c, int cutLen, double radioT)
{
	conts = c;
	lenThresh = cutLen;
	radioThresh = radioT;
};
void contours::cutShortContours()
{
	vector<vector<Point> > temp;
	for (int i = 1; i < conts.size(); i++)
	{
		if (conts[i].size() > lenThresh)
		{
			temp.push_back(conts[i]);
		}
	}
	conts = temp;
}

void contours::sortContours()
{
	for (int i=0;i<conts.size();i++)
	{
		//remove repeat points
		vector<Point> temp0;
		//cout << conts[i] << endl;
		temp0.push_back(conts[i][0]);
		for (int j=1;j<conts[i].size();j++)
		{
			bool flag = 1;
			for (int k=0;k<temp0.size();k++)
			{
				if (conts[i][j].x==temp0[k].x && conts[i][j].y == temp0[k].y)
				{
					flag = 0;
					break;
				}
			}
			if (flag==1)
			{
				temp0.push_back(conts[i][j]);
			}
		}
		conts[i] = temp0;
		int* isVisited = new int[conts[i].size()];
		isVisited[0] = 1;
		Point beginP = findEndPoints(conts[i]);
		vector<Point> temp;
		temp.push_back(beginP);
		bool flag = 1;
		while (flag)
		{
			flag = 0;
			for (int j = 0; j < conts[i].size(); j++)
			{
				if (isVisited[j] != 1 && isNeigbour(temp[temp.size()-1], conts[i][j]))
				{
					isVisited[j] = 1;
					temp.push_back(conts[i][j]);
					flag = 1;
					break;
				}
			}
		}		
		conts[i] = temp0;
		delete[] isVisited;
	}
}

void contours::filterContours()
{
	vector<vector<Point> > temp;
	for (int i = 0; i < conts.size();i++)
	{
		int* distOriten = new int[conts[i].size()];//距离增加为1，减小为0
		distOriten[0] = 1;
		for (int j=1;j<conts[i].size();j++)
		{
			if (distance(conts[i][j], conts[i][0]) >= distance(conts[i][j - 1], conts[i][0]))
			{
				distOriten[j] = 1;
			}
			else
				distOriten[j] = 0;
		}
		int backNums = 0;
		for (int i=1;i<conts[i].size();i++)
		{
			if (distOriten[i]!=distOriten[i-1])
			{
				backNums++;
			}
		}
		double rt = backNums*1.0 / conts.size();
		if (rt<radioThresh)
		{
			temp.push_back(conts[i]);
		}
		delete[] distOriten;
	}
	conts = temp;
}

Point contours::findEndPoints(vector<Point>& vp)
{
	int minx = vp[0].x;
	int miny = vp[0].y;

	for (int i=1;i<vp.size();i++)
	{
		if (vp[i].x<minx)
		{
			minx = vp[i].x;
			miny = vp[i].y;
		}
	}
	return Point(minx, miny);
}

bool contours::isNeigbour(Point& a, Point& b)
{
	if (a.x==b.x && a.y==b.y)
	{
		return 0;
	}
	if (abs(a.x-b.x)<2 && abs(a.y-b.y)<2)
	{
		return 1;
	}
	return 0;
}

int contours::distance(Point& a, Point& b)
{
	return abs(a.x - b.x) + abs(a.y - b.y);
}
