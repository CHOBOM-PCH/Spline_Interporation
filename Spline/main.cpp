#include "Spline.h"
#include <iostream>




int main()
{
	//DWORD dwStart = timeGetTime();
	int time=1000;
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> x_m;
	std::vector<double> y_m;
	std::vector<double>* x_z;
	std::vector<double>* y_z;

	std::vector<double> x_s;
	std::vector<double> y_s;
	std::vector<double>* x_v;
	std::vector<double>* y_v;

	//Mat plot_img(512,512,CV_8UC3, Scalar(0,0,0));

	int z = 0;
	for (int i = 0; i <= 5;i++){
		z = 10 + 100*i;
		x.push_back(z);
	}
	y.push_back(120);
	y.push_back(400);
	y.push_back(50);
	y.push_back(220);
	y.push_back(120);
	y.push_back(130);

	x_z = &x_m;
	y_z = &y_m;

	mono_spline(time,x,y,x_z,y_z);

	x_v = &x_s;
	y_v = &y_s;
	cubic_spline(x,y,x_v,y_v);

	for (int k = 0;k<x.size();k++){//입력된 값에 의한 그래프
		printf("입력된 값 x = %lf , y = %lf \n" ,x[k],y[k]);
	//	if (k == 0);

	//	else
	//		line(plot_img,Point(x[k-1],y[k-1]),Point(x[k],y[k]),Scalar(0,255,255),2,8,0);
	}

	for (int k =0;k<x_m.size();k++){//mono spline후 그래프
		printf("mono계산된 값 x[%d] = %lf , y[%d] = %lf \n", k, x_m[k], k, y_m[k]);
		//if (k == 0);

	//	else
	//		line(plot_img,Point(x_m[k-1],y_m[k-1]),Point(x_m[k],y_m[k]),Scalar(255,255.0),2,8,0);
	}
	for (int k =0;k<x_s.size()-1;k++){//cubic spline후 그래프
		printf("cubic 계산된 값 x[%d] = %lf , y[%d] = %lf \n", k, x_s[k], k, y_s[k]);
		//if (k == 0);

	//	else
	//		line(plot_img,Point(x_s[k-1],y_s[k-1]),Point(x_s[k],y_s[k]),Scalar(255,100,255),2,8,0);
	}
	//imshow("PLOT",plot_img);

	//DWORD dwEnd = timeGetTime();
	//printf("end %d ms \n", dwEnd-dwStart);

	//waitKey(0);
	system("pause");
	return 0;
}