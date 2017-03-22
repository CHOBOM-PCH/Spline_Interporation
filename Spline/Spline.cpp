#include "Spline.h"



double h00(double t)
{
	return 2*t*t*t - 3*t*t +1;
}
double h10(double t)
{
	return t*(1-t)*(1-t);
}
double h01(double t)
{
	return t*t*(3-2*t);
}
double h11(double t)
{
	return t*t*(t-1);
}

bool mono_spline(int time_limit,//monotone ���ö��� ������//��ũ ����� ���� �ð� ����
	const std::vector<double> x_mrc,//x ��ǥ��
	const std::vector<double> y_mrc,//y ��ǥ��
	std::vector<double>* destX,
	std::vector<double>* destY)
{
	// 0-based index ���.
	int n = (int)x_mrc.size();//x��ǥ�� ���� n���� x��ǥ�� �������� x_(n-1)�̵�
	int k = 0;//n�� ����?
	double *m = new double[n];//���� m
	m[0] = (y_mrc[1] - y_mrc[0])/(x_mrc[1] - x_mrc[0]);//ó�� x�� ���õ� ����
	m[n-1] = (y_mrc[n-1] - y_mrc[n-2])/(x_mrc[n-1]-x_mrc[n-2]);//������ x����

	for (k = 1; k<n-1; k++){//��x��ǥ�� ���� y������ ����
		m[k] = (y_mrc[k] - y_mrc[k-1])/(2*(x_mrc[k] - x_mrc[k-1])) + (y_mrc[k+1] - y_mrc[k])/(2*(x_mrc[k+1]-x_mrc[k]));
	}
	for (k = 0; k<n-1; k++){
		double delta_k = (y_mrc[k+1]-y_mrc[k])/(x_mrc[k+1]-x_mrc[k]);

		if (fabs(delta_k) <= eps){
			m[k] = m[k+1] = 0;
		}
		else{
			double ak = m[k]/delta_k;
			double bk = m[k+1]/delta_k;
			if(ak*ak + bk*bk > 9){//������Ʈ �Ǻ�: 9���� ũ�� ��ũk�� ���� 1���� �۴ٴ¸�
				m[k] = 3/(sqrt(ak*ak+bk*bk))*ak*delta_k;//m_k�� m_k+1�� ���Ӱ� ���
				m[k+1] = 3/(sqrt(ak*ak+bk*bk))*bk*delta_k;//�̴� ���� a�� bũ�⸦ ������ 3�� ���� ũ��� ����
			}
		}
	}

	for (k = 0; k<n-1; k++){//cubic ������
		double cur_x = (double)((int)(0.5 + x_mrc[k]));//x_lower current
		double next_x = (double)((int)(x_mrc[k+1]));//x_upper
		double cur_y = y_mrc[k];
		double next_y = y_mrc[k+1];
		double h = next_x - cur_x;
		double x = 0;
		double inc = (next_x - cur_x)*0.1;//���� x���� ���� x�����̿��� 10�����ؼ� ����������
		for (x = cur_x; x<next_x; x+=inc){
			if (x > time_limit) break;
			double t = (x-cur_x)/h;
			if (destX != NULL){
				destX->push_back(x);//destX���ڿ� ���� x�� �������
			}
			double y = cur_y*h00(t) + h*m[k]*h10(t) + next_y*h01(t) + h*m[k+1]*h11(t);
			destY->push_back(y);
		}
	}

	delete m;

	return true;
}



bool cubic_spline(std::vector<double>x_meries,//cubic ���ö��� ������
	std::vector<double>y_meries,
	std::vector<double> *destX,
	std::vector<double>* destY)
{   
	int n = x_meries.size() -1;
	// Step 1.
	double *h = new double[n+1];
	double *alpha = new double[n+1];
	int i = 0;
	for(i = 0; i<=n-1; i++){
		h[i] = (x_meries)[i+1] - (x_meries)[i];
	}

	// Step 2.
	for(i = 1; i<=n-1;i++){
		alpha[i]= 3*((y_meries)[i+1]-(y_meries)[i])/h[i]-3*((y_meries)[i]-(y_meries)[i-1])/h[i-1];
	}

	// Step 3.
	double *l = new double[n+1];
	double *u = new double[n+1];
	double *z = new double[n+1];
	double *c = new double[n+1];
	double *b = new double[n+1];
	double *d = new double[n+1];

	l[0] = 1; u[0] = 0; z[0] = 0;

	// Step 4.
	for(i = 1; i<=n-1; i++){
		l[i] = 2*((x_meries)[i+1] - (x_meries)[i-1]) - h[i-1]*u[i-1];
		u[i] = h[i]/l[i];
		z[i] = (alpha[i] - h[i-1]*z[i-1]) / l[i];
	}

	// Step 5.
	l[n] = 1;     z[n] = 0;     c[n] = 0;

	// Step 6.
	for(i = n-1; i>=0; i--){
		c[i] = z[i] - u[i]*c[i+1];
		b[i] = ((y_meries)[i+1] - (y_meries)[i])/h[i] - h[i]*(c[i+1] + 2*c[i])/3;
		d[i] = (c[i+1] - c[i]) / (3*h[i]);
	}

	for(i = 0; i<=n-2;i++){
		double x = (x_meries)[i];
		double inc = ((x_meries)[i+1] - (x_meries)[i])*0.1;//1/10������ ���
		for(; x < (x_meries)[i+1]; x+=inc){
			double x_offset = x - (x_meries)[i];
			double Sx = (y_meries)[i] + b[i]*x_offset + c[i]*x_offset*x_offset + d[i]*x_offset*x_offset*x_offset;
			if(destX != NULL){
				destX->push_back(x);
			}
			if(destY != NULL){
				destY->push_back(Sx);
			}
		}           
	}
	delete [] h;
	delete [] alpha;
	delete [] l;
	delete [] u;
	delete [] z;
	delete [] c;
	delete [] b;
	delete [] d;


	return true;
}