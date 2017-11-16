#include "ContourDetecion.h"
#include <vector>
#include <iostream>

ContourDetecion::ContourDetecion(Mat& srcImg)
{
	src = srcImg;
}

ContourDetecion::~ContourDetecion()
{
}

void ContourDetecion::GaborDirection(Mat& srcImg, Mat& dstImg, Mat& dirImg, double lamda1, double gama1, int size, double sigmaRatio)
{
	//srcImg.convertTo(srcImg, CV_32FC1);
	dirImg.convertTo(dirImg, CV_32FC1);

	Mat gaborFilter(srcImg.size(),CV_32F);
	//Mat gaborReal;
	//Mat gaborImag;
	std::vector<Mat> vgaborReal;
	std::vector<Mat> vgaborImag;
	int directions = 3;

	double deltaA = 3.14159 / directions;

	
	for (int i = 0; i < srcImg.rows; i++)
	{
		float* ptr = gaborFilter.ptr<float>(i);
		float* ptrDir = dirImg.ptr<float>(i);

		for (int j = 0; j < srcImg.cols; j++)
		{
			ptr[j] = 0;
			ptrDir[j] = 0;
		}
	}

	Mat tmp1(srcImg.size(), CV_8U);
	Mat tmp2(srcImg.size(), CV_8U);
	for (int k = 0; k < directions; k++)
	{
		Mat gaborReal(size, size, CV_32FC1);
		Mat gaborImag(size, size, CV_32FC1);
		GaborGenerator(gaborReal, gaborImag, k*deltaA, lamda1, gama1, size, sigmaRatio);
		vgaborReal.push_back(gaborReal);
		vgaborImag.push_back(gaborImag);
		filter2D(srcImg, tmp1, srcImg.depth(), vgaborReal[k]);
		filter2D(srcImg, tmp2, srcImg.depth(), vgaborImag[k]);
		for (int i=0;i<srcImg.rows;i++)
		{
			float* ptr = gaborFilter.ptr<float>(i);
			uchar* ptr1 = tmp1.ptr<uchar>(i);
			uchar* ptr2 = tmp2.ptr<uchar>(i);
			float* ptrDir = dirImg.ptr<float>(i);
			for (int j=0;j<srcImg.cols;j++)
			{
				double tmpSqrt = ptr2[j];// sqrt(ptr1[j] * ptr1[j] + ptr2[j] * ptr2[j]);
				if (ptr[j] < tmpSqrt)
				{
					ptr[j] = tmpSqrt;
					ptrDir[j] = k*deltaA;
				}
			}
		}						
		/*imshow("tmp1", tmp1);
		imshow("tmp2", tmp2);
		imshow("gaborFilter", gaborFilter);
		waitKey(0);*/
		
	}
	normalize(gaborFilter, gaborFilter, 0, 1, NORM_MINMAX);
	dstImg = gaborFilter;

	//for (int i = size / 2; i < src.rows - size / 2; i++)
	//{
	//	uchar* ptrSrc = srcImg.ptr<uchar>(i);
	//	for (int j = size / 2; j < src.cols - size / 2; j++)
	//	{
	//		double maxConv = 0;
	//		int maxIdx = 0;
	//		//求出16个方向中使得卷积值最大的方向,即最优响应方向
	//		for (int k = 0; k < 8; k++)
	//		{
	//			double tmpSumReal = 0;
	//			double tmpSumImag = 0;
	//			float* ptrReal = vgaborReal[k].ptr<float>(int(size / 2));
	//			float* ptrImag = vgaborImag[k].ptr<float>(int(size / 2));

	//			for (int m = -size / 2; m <= size / 2; m++)
	//			{
	//				for (int n = -size / 2; n <= size / 2; n++)
	//				{
	//					//tmpSumReal += (ptrSrc + m)[j + n] * (vgaborReal[k].at<float>((int)size/2 + m, (int)size / 2 + n));
	//					//tmpSumImag += (ptrSrc + m)[j + n] * (vgaborImag[k].at<float>((int)size / 2 + m, (int)size / 2 + n));
	//					tmpSumReal += (ptrSrc + m)[j + n] * ((ptrReal + m)[int(size / 2)] + n);
	//					tmpSumImag += (ptrSrc + m)[j + n] * ((ptrImag + m)[int(size / 2)] + n);
	//				}
	//			}
	//			double tmp = tmpSumImag;// sqrt(tmpSumReal * tmpSumReal + tmpSumImag * tmpSumImag);
	//			if (maxConv < tmpSumImag)
	//			{
	//				maxConv = tmpSumImag;
	//				maxIdx = k;
	//			}
	//		}
	//		gaborFilter.at<uchar>(i, j) = maxConv;
	//		dirImg.at<uchar>(i, j) = maxIdx;
	//	}		
	//}
	//gaborFilter.convertTo(gaborFilter, CV_8U);

	//imshow("gaborFilter", gaborFilter);
	//waitKey(0);
}

void ContourDetecion::GaborGenerator(Mat& gaborReal, Mat& gaborImag, double theta, double lamda1, double gama1, int size, double sigmaRatio)
{
	double x1, y1;
	double lamda = lamda1;
	double gama = gama1;
	double sigma = sigmaRatio*lamda;	

	//Mat gaborDerect(size, size, CV_32FC1);
	double sumReal = 0;
	double sumImag = 0;

	for (int i = -size / 2; i <= size / 2; i++)
	{
		float * pReal = gaborReal.ptr<float>(i + size / 2);
		float* pImag = gaborImag.ptr<float>(i + size / 2);
		//float* pDerect = gaborImag.ptr<float>(i + size / 2);

		for (int j = -size / 2; j <= size / 2; j++)
		{
			x1 = i*cos(theta) + j*sin(theta);
			y1 = -i*sin(theta) + j*cos(theta);
			double real = exp(-(x1*x1 + gama*gama*y1*y1) / (2 * sigma*sigma))*cos(2 * 3.14159*x1 / lamda);
			double imag = exp(-(x1*x1 + gama*gama*y1*y1) / (2 * sigma*sigma))*sin(2 * 3.14159*x1 / lamda);

			pReal[j + size / 2] = real;
			pImag[j + size / 2] = imag;		
			sumReal += real;
			sumImag += imag;

		}
	}
	//gaborReal = gaborReal / sumReal;
	//gaborImag = gaborImag / sumImag;

	//for test
	//Mat tmp1(src.size(), CV_8U);
	//Mat tmp2(src.size(), CV_8U);
	//Mat tmp(src.size(), CV_8U);

	//filter2D(src, tmp2, src.depth(), gaborImag);
	//filter2D(src, tmp1, src.depth(), gaborReal);
	//for (int i = 0; i < tmp1.rows; i++)
	//{
	//	uchar* ptr1 = tmp1.ptr<uchar>(i);
	//	uchar* ptr2 = tmp2.ptr<uchar>(i);
	//	uchar* ptr = tmp.ptr<uchar>(i);

	//	for (int j = 0; j < tmp1.cols; j++)
	//	{
	//		ptr[j] = sqrt(ptr1[j] * ptr1[j] + ptr2[j] * ptr2[j]);
	//	}
	//}
	//imshow("tmp test for gabor gene", tmp);
	//normalize(tmp1, tmp1, 0, 1, NORM_MINMAX);
	//imshow("tmp11", tmp1);
	//imshow("tmp22", tmp2);
	//////namedWindow("gaborReal", WINDOW_NORMAL);
	//////namedWindow("gaborImag", WINDOW_NORMAL);

	/////*imshow("GaborReal", gaborReal);
	////imshow("GaborImag", gaborImag);*/


	//waitKey(0);

}

void ContourDetecion::GaborSuppression(Mat & srcImg, Mat& dirImg, Mat & dstImg, int supSize)
{
	if (supSize % 2 == 0)
		supSize += 1;
	for (int i=0;i<srcImg.rows;i++)
	{
		float* ptrDst = dstImg.ptr<float>(i);

		for (int j=0;j<srcImg.cols;j++)
		{
			ptrDst[j] = 0;
		}
	}	
	for (int i = supSize / 2; i < srcImg.rows - supSize / 2; i++)
	{		
		float* ptrDir = dirImg.ptr<float>(i);
		float* ptrSrc = srcImg.ptr<float>(i);
		float* ptrDst = dstImg.ptr<float>(i);


		for (int j = supSize / 2; j < srcImg.cols - supSize / 2; j++)
		{			
			
			for (int m = -supSize / 2; m <= supSize / 2; m++)
			{
				for (int n = -supSize; n <= supSize / 2; n++)
				{
					if (abs(n / (m*1.0)) > (CV_PI * 3 / 8.0) || sqrt(m*m+n*n) < supSize/6)//非侧向抑制区
					{
						continue;
					}
					double coefDist = exp(-(m*m + n*n) / (2 * 1.5*1.5));
					double deltaAng = min(abs((ptrDir + m)[j + n] - ptrDir[j]), (float)CV_PI - abs((ptrDir + m)[j + n] - ptrDir[j]));
					double coefAng = exp(-(deltaAng*deltaAng) / (2 * supSize*supSize));
					ptrDst[j] += (ptrSrc + m)[j + n] * coefAng*coefDist;
				}
			}
				
		}
	}
	imshow("srcImg in supp", srcImg);
	normalize(dstImg, dstImg, 0, 1, NORM_MINMAX);
	imshow("supp", dstImg);
	waitKey(0);
}

Mat ContourDetecion::GaborEnergy(Mat & srcImg, Mat& dirImg, vector<Mat >& vGaborE, double lamda1, double gama1, int size, double sigmaRatio)
{
	int directions = 12;
	double deltaA = 3.14159 / directions;
	std::vector<Mat> vgaborReal;
	std::vector<Mat> vgaborImag;
	Mat tmp1(srcImg.size(), CV_8U);
	Mat tmp2(srcImg.size(), CV_8U);

	Mat maxGaborE(srcImg.size(), CV_32F);
	double maxE = 0;
	for (int k = 0; k < directions; k++)
	{
		Mat tmpRet(srcImg.size(), CV_32FC1);
		Mat gaborReal(size, size, CV_32FC1);
		Mat gaborImag(size, size, CV_32FC1);
		GaborGenerator(gaborReal, gaborImag, k*deltaA, lamda1, gama1, size, sigmaRatio);		
		filter2D(srcImg, tmp1, srcImg.depth(), gaborReal);
		filter2D(srcImg, tmp2, srcImg.depth(), gaborImag);

		//for test
		//tmp1 = tmp1 / 255.0;
		//tmp2 = tmp2 / 255.0;
		//imshow("tmpReal", tmp1);
		//imshow("tmpImag", tmp2);
		//waitKey(0);

		for (int i = 0; i<srcImg.rows; i++)
		{
			float* ptr1 = tmp1.ptr<float>(i);
			float* ptr2 = tmp2.ptr<float>(i);
			float* ptrRet = tmpRet.ptr<float>(i);
			float* ptrMax = maxGaborE.ptr<float>(i);
			float* ptrDir = dirImg.ptr<float>(i);
			for (int j = 0; j<srcImg.cols; j++)
			{
				double tmpSqrt = sqrt(ptr1[j] * ptr1[j] + ptr2[j] * ptr2[j]);
				ptrRet[j] = tmpSqrt;
				if (maxE < tmpSqrt)
				{
					maxE = tmpSqrt;
					ptrMax[j] = tmpSqrt;
					ptrDir[j] = i;
				}
			}
		}	
		//normalize(tmpRet, tmpRet, 0, 1, NORM_MINMAX);

		vGaborE.push_back(tmpRet);
		/*imshow("vGaborE[i]", tmpRet);
		waitKey(0);*/
	}	
	return maxGaborE;
}

void ContourDetecion::FacilitationEnergy(vector<Mat>& vGaborE, vector<Mat>& faciEnergy, Mat& maxGabor, Mat& dirImg, int size)
{
	int Re =10;
	double psi = CV_PI / 6;
	double sigma_c = 0.17;
	double sigma_d = 7; 




	//尝试用数组保存Matrix,加快速度
	const int rows = src.rows;
	const int cols = src.cols;
	vector<float** > vGaborEarr;

	for (int i = 0; i < 12; i++)
	{
		float** tmpArr = new float*[rows];
		for (int i = 0; i < rows; i++)
		{
			tmpArr[i] = new float[cols];
		}
		for (int r = 0; r < src.rows; r++)
		{
			for (int c = 0; c < src.cols; c++)
			{
				tmpArr[r][c] = vGaborE[i].at<float>(r, c);
			}
		}
		vGaborEarr.push_back(tmpArr);
// 		for (int i = 0; i < rows; i++)
// 		{
// 			delete[] tmpArr[i];
// 		}
// 		delete[] tmpArr;
	}

	//======
	//Mat MatF(src.size(), CV_32FC1, Scalar(0));
	//
	//for (int m = Re; m < src.rows - Re; m++)
	//{
	//	float* ptrF = MatF.ptr<float>(m);

	//	for (int n = Re; n < src.cols - Re; n++)
	//	{
	//		double supValue = 0;
	//		double alpha = dirImg.at<float>(m, n) * CV_PI / 12;
	//		for (int k = -Re; k < Re; k++)
	//		{

	//			for (int l = -Re; l < Re; l++)
	//			{
	//				double tmpDist = sqrt(k*k + l*l);
	//				double wd = exp(-tmpDist*tmpDist / (2 * sigma_d*sigma_d));					
	//				double beta = 0;
	//				double gama = atan(l / (k*1.0));
	//				if (gama < 0)
	//					gama += CV_PI;
	//				if (tmpDist > Re || abs(gama - alpha) > psi)
	//				{
	//					cofDistMat[Re + k][Re + l] = 0;
	//					cofAngMat[Re + k][Re + l] = 0;
	//					idxGaborE[Re + k][Re + l] = 0;
	//					continue;
	//				}
	//				tmp = 2 * gama - alpha;
	//				if (tmp < 0)
	//					beta = tmp + CV_PI;
	//				if (tmp >= 0 && tmp < CV_PI)
	//					beta = tmp;
	//				if (tmp >= CV_PI)
	//					beta = tmp - CV_PI;
	//				int betaIdx = (beta + 0.5*CV_PI / 12) / (CV_PI / 12);
	//				betaIdx = betaIdx % 12;
	//				double curv = 0;
	//				if (tmp >= 0 && tmp < CV_PI)
	//				{
	//					curv = 2 * sin(abs((beta - alpha) / 2.0)) / tmpDist;
	//				}
	//				else
	//				{
	//					curv = 2 * cos(abs((beta - alpha) / 2.0)) / tmpDist;
	//				}
	//				double wc = exp(-curv*curv / (2 * sigma_c*sigma_c));
	//				supValue += wd * wc * vGaborEarr[idxGaborE[Re + k][Re + l]][Re + k][Re + l];

	//			}
	//		}
	//		ptrF[n] = supValue;

	//		//for test
	//		//std::cout << "supValue:" << supValue << std::endl;
	//	}
	//}

	for (int i = 0;i < 12;i++)//alpha
	{
		Mat MatF(src.size(), CV_32FC1, Scalar(0));
		double alpha = i*CV_PI / 12;

		//======
		double tmpDist;
		double beta = 0;
		double gama;
		double tmp;
		double curv;
		double wc;
		double wd;

		Mat cofAng(20, 20, CV_32F);
		double cofDistMat[20][20];
		double cofAngMat[20][20];
		int idxGaborE[20][20];
		for (int k = -Re; k < Re; k++)
		{
			for (int l = -Re; l < Re; l++)
			{
				tmpDist = sqrt(k*k + l*l);
				wd = exp(-tmpDist*tmpDist / (2 * sigma_d*sigma_d));
				cofDistMat[Re + k][Re + l] = wd;
				beta = 0;
				gama = atan(-k / (l*1.0));//坐标系容易混乱，此处把坐标系转为常规坐标系
				if (gama < 0)
					gama += CV_PI;
				if (tmpDist > Re || abs(gama - alpha) > psi)
				{
					cofDistMat[Re + k][Re + l] = 0;
					cofAngMat[Re + k][Re + l] = 0;
					idxGaborE[Re + k][Re + l] = 0;
					cofAng.at<float>(Re + k, Re + l) = 0;
					continue;
				}
				tmp = 2 * gama - alpha;
				if (tmp < 0)
					beta = tmp + CV_PI;
				if (tmp >= 0 && tmp < CV_PI)
					beta = tmp;
				if (tmp >= CV_PI)
					beta = tmp - CV_PI;
				int betaIdx = (beta + 0.5*CV_PI / 12) / (CV_PI / 12);
				betaIdx = betaIdx % 12;

				if (tmp >= 0 && tmp < CV_PI)
				{
					curv = 2 * sin(abs((beta - alpha) / 2.0)) / tmpDist;
				}
				else
				{
					curv = 2 * cos(abs((beta - alpha) / 2.0)) / tmpDist;
				}
				wc = exp(-curv*curv / (2 * sigma_c*sigma_c));
				cofAngMat[Re + k][Re + l] = wc;
				idxGaborE[Re + k][Re + l] = betaIdx;
				cofAng.at<float>(Re+k, Re+l) = wc;
			}
		}
		/*std::cout << cofAng << std::endl;
		normalize(cofAng, cofAng, 0, 1, NORM_MINMAX);
		imshow("cofang", cofAng);
		waitKey(0);*/
		for (int m=Re;m<src.rows-Re;m++)
		{
			float* ptrF = MatF.ptr<float>(m);
			
			for (int n=Re;n<src.cols-Re;n++)
			{
				double supValue = 0;
				
				for (int k = -Re; k<Re ; k++)
				{
					
					for (int l=-Re;l<Re;l++)
					{
						//float* ptrIdx = vGaborE[idxGaborE[Re + k][Re + l]].ptr<float>(m + k);
						supValue += cofDistMat[Re + k][Re + l] * cofAngMat[Re + k][Re + l] * vGaborE[idxGaborE[Re + k][Re + l]].at<float>(m + k, n + l);
						//supValue += cofDistMat[Re + k][Re + l] * cofAngMat[Re + k][Re + l] * vGaborEarr[idxGaborE[Re + k][Re + l]][m+ k][n+ l];

					}
				}
				ptrF[n] = supValue;
				
				//for test
				//std::cout << "supValue:" << supValue << std::endl;
			}
		}
		/*normalize(MatF, MatF, 0, 1, NORM_MINMAX);
		imshow("MatF", MatF);
		waitKey(0);*/
		faciEnergy.push_back(MatF);
	}
}

void ContourDetecion::FacProgress(Mat & srcImg, Mat & dstImg)
{
	srcImg.convertTo(srcImg, CV_32F);
	vector<Mat > vGaborE;
	vector<Mat > vFaciEnergy;
	vector<Mat > vDst;
	double cofFaci = 0.5;
	Mat tmpSrc = srcImg.clone();
	Mat maxGaborE;
	Mat dirImg(src.size(), CV_32F);
	maxGaborE = GaborEnergy(tmpSrc,dirImg,vGaborE);	
	FacilitationEnergy(vGaborE, vFaciEnergy, maxGaborE,dirImg);
	for (int i = 0; i < 12; i++)
	{
// 		imshow("vFaciE", vFaciEnergy[i]);
// 		waitKey(0);
	}
	for (int i=0;i<12;i++)
	{
		vDst.push_back(vGaborE[i] + cofFaci * vFaciEnergy[i]);
	}
	for (int i=0;i<srcImg.rows;i++)
	{
		float* ptrDst = dstImg.ptr<float>(i);
		for (int j = 0; j < srcImg.cols; j++)
		{
			double max = 0;
			for (int v=0;v<12;v++)
			{
				double tmp = vDst[v].at<float>(i, j);
				if (max < tmp)
				{
					max = tmp;
				}
			}
			ptrDst[j] = max;
		}
	}
	normalize(dstImg, dstImg, 0, 1, NORM_MINMAX);
	imshow("dstImg", dstImg);
	waitKey(0);
}

void ContourDetecion::MergeRI(Mat& srcImg, Mat& dstImg, Mat& gaborReal, Mat& gaborImag)
{
	Mat dst1(src.size(), CV_32FC1);
	Mat dst2(src.size(), CV_32FC1);

	filter2D(srcImg, dst1, dst2.depth(), gaborReal);
	filter2D(srcImg, dst2, dst2.depth(), gaborImag);

	imshow("gaborReal", dst1);
	imshow("gaborImag", dst2);

	for (int i = 0; i < dst2.rows; i++)
	{
		float* ptr = dstImg.ptr<float>(i);
		float* ptr1 = dst1.ptr<float>(i);
		float* ptr2 = dst2.ptr<float>(i);

		for (int j = 0; j < dst2.cols; j++)
		{
			ptr[j] = sqrt(ptr1[j] * ptr1[j] + ptr2[j] * ptr2[j]);
		}
	}
}