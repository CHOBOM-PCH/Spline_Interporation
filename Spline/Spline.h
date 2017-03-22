#include <math.h>
#include <vector>
//#include <iostream>
#define eps 1*(10^(-10))

bool mono_spline(int time_limit,//monotone ���ö��� ������//��ũ ����� ���� �ð� ����
	const std::vector<double> x_mrc,//x ��ǥ��
	const std::vector<double> y_mrc,//y ��ǥ��
	std::vector<double>* destX,
	std::vector<double>* destY);

bool cubic_spline(std::vector<double>x_meries,//cubic ���ö��� ������
	std::vector<double>y_meries,
	std::vector<double>* destX,
	std::vector<double>* destY);

