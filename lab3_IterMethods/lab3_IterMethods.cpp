// lab3_IterMethods.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include "MatrixOperation.h"
#include "VectorsOperation.h"
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;

const double epsylon = 0.0001;

double tauCalc(vector<vector<double>> a)
{
	double y = 0.9;
	return y*(2 / EuclideanNorm(MultMatrix(Transpose(a),a)));
}

double qCalc(vector<double> xPrev, vector<double> xCur, vector<double> xNext)
{
	return FirstVectorNorm(SubtractionVector(xNext, xCur)) / FirstVectorNorm(SubtractionVector(xCur,xPrev));
}

double Pogreshnost(vector<double> xCur, vector<double> xNext, double q)
{
	return  q*FirstVectorNorm(SubtractionVector(xNext, xCur))/(1-q);
}

vector<double> vectorNevyazki(vector<vector<double>> a, vector<double> x, vector<double> b)
{
	vector<double> n(b.size());

	for (int i = 0; i < a.size(); ++i)
	{
		for (int j = 0; j < a[i].size(); ++j)
			n[i] += a[i][j] * x[j];
		n[i] -= b[i];
	}

	return n;
}

vector<double> sqrtMethod(vector<vector<double>> a, vector<double> b)
{
	vector<vector<double>> S(a.size(), (vector<double>(a.size())));
	vector<vector<double>> d(a.size(),(vector<double>(b.size())));


	for(int i=0;i<S.size();++i)
	{
		double sum = 0;

		for (int k = 0; k < i; ++k)
			sum += S[k][i] * S[k][i]*d[k][k];

		d[i][i] = copysign(1.0, a[i][i] - sum);
		S[i][i] = sqrt(a[i][i] - sum);

		double temp = 1 / (S[i][i] * d[i][i]);

		for(int j=i+1;j<S.size();++j)
		{
			sum = 0;

			for (int k = 0; k < i; ++k)
				sum += S[k][i] * S[k][j]*d[k][k];

			S[i][j] = (a[i][j] - sum) *temp;
		}

	}

	vector<double> y(b.size());
	vector<double> x(b.size());

	y[0] = b[0] / (S[0][0]*d[0][0]);

	for (int i = 1; i < b.size(); ++i)
	{
		double sum = 0;
		for (int j = 0; j < i; ++j)
			sum += S[j][i] * y[j]*d[j][j];

		y[i] = (b[i] - sum) / (S[i][i]*d[i][i]);
	}

	x[b.size() - 1] = y[b.size() - 1] / S[b.size() - 1][b.size() - 1];

	for (int i = x.size() - 2; i > -1; --i)
	{
		double sum = 0;
		for (int j = i; j<x.size(); ++j)
			sum += S[i][j] * x[j];

		x[i] = (y[i] - sum) / S[i][i];
	}

	return x;
}

void simpleIteration(vector<vector<double>> a, vector<double> b,vector<double> SolvedX)
{
	vector<double> xPrev(4,0);
	vector<double> xCur=b;
	vector<double> xNext(4);
	vector<double> nev;
	double nevNorm, tau=tauCalc(a), q=0, pogr=1;
	double pogrNorm = 1;

	for(int k=1;abs(pogrNorm)-epsylon>0;++k)
	{
		for(int i=0;i<a.size();++i)
		{
			double sum = 0;
			for (int j = 0; j < a[i].size(); ++j)
				sum += a[i][j] * xCur[j];
			xNext[i] = xCur[i] + tau * (b[i] - sum);
		}


		nev = vectorNevyazki(a, xCur, b);
		nevNorm = FirstVectorNorm(nev);
		q =(k==1)?FirstVectorNorm(SubtractionVector(xCur,xNext)):qCalc(xPrev, xCur, xNext);
		pogr = Pogreshnost(xCur, xNext,q);
		pogrNorm = FirstVectorNorm(SubtractionVector(xNext, SolvedX));
		cout <<k<<" | " << setprecision(4) << tau <<setprecision(3) <<" | " << q << " | " <<setw(10) <<setprecision(7) << nevNorm << " | " << pogrNorm<<" | " << pogr << endl;

		xPrev = xCur;
		xCur = xNext;

		
	}
}

void fast_gardient(vector<vector<double>> A, vector<double> b,vector<double> X_q) //	Ìåòîä ñêîðåéøåãî ñïóñêà
{

	cout << "Ìåòîä íàèñêîðåéøåãî ñïóñêà" << endl;
	cout << setw(4) << "Iter|" << setw(8) << "tau|" << setw(8) << "q|" << setw(12)
		<< "Íîðìà r|" << setw(12) << "Íîðìà ïîãð|" << setw(12) << "Îöåíêà ïîãð|"
		<< setw(8) << "X[0]" << setw(8) << "X[1]" << setw(8) << "X[2]" << setw(8) << "X[3]" << endl;


	vector<double> X(A.size()); // X[k+1]
	vector<double> X_pred(A.size()); // X[k]
	vector<double> X_prpr(A.size()); //X[k-1]

	X_pred = b;
	X_prpr = b;
	
	int iter = 0;

	double t, q, pg, norm_pg, tau;

	vector<double> r(A.size()); // âåêòîð íåâÿçêè


	do {


		r = vectorNevyazki(A, X_pred, b);
		tau = ScalarMult(r, r) / ScalarMult(MultMatrixVector(A, r), r);


		X = AdditionVector(MultVectorNum(SubtractionVector(b, MultMatrixVector(A, X_pred)), tau), X_pred);


		q = (iter == 0) ? FirstVectorNorm(SubtractionVector(X_pred, X)) : qCalc(X_prpr, X_pred, X);
		pg = Pogreshnost(X_pred, X, q);
		norm_pg = FirstVectorNorm(SubtractionVector(X, X_q));

		X_prpr = X_pred;
		X_pred = X;

		cout << setw(4) << iter + 1 << "|" << fixed << setprecision(4) << tau << " | " << setprecision(3) << q
			<< " |" << setprecision(7) << setw(10) << FirstVectorNorm(r) << " |" << setw(10) << norm_pg << " |" << setw(10) << pg << " | ";

		for (int i = 0; i < b.size(); i++)
		{
			cout << setw(8) << fixed << setprecision(5) << X[i];
		}
		cout << endl;
		iter++;

	} while (FirstVectorNorm(SubtractionVector(X, X_q)) / FirstVectorNorm(X) > epsylon);
}

void conjugate_gardient(vector<vector<double>> A, vector<double> b, vector<double> X_q) //ìåòîä ñîïðÿæ¸ííûõ ãðàäèåíòîâ
{

	cout << "Ìåòîä ñîïðÿæåííûé ãðàäèåòîâ" << endl;
	cout << setw(4) << "Iter|" << setw(8) << "tau|" << setw(8) << "q|" << setw(12) << "Íîðìà r|"
		<< setw(12) << "Íîðìà ïîãð|" << setw(12) << "Îöåíêà ïîãð|"
		<< setw(8) << "X[0]" << setw(8) << "X[1]" << setw(8) << "X[2]" << setw(8) << "X[3]" << endl;

	vector<double> X(A.size()); // X[k]
	vector<double> X_now(A.size()); //X[k+1]
	vector<double> X_pred(A.size()); // X[k-1]

	vector<double> r(A.size()); // âåêòîð íåâÿçêè

	X_pred = b;
	X_now = b;
	X = b;

	int iter = 0;
	double tau, tau_old, alpha = 1, alpha_old = 1, pg, norm_pg, q;

	r = vectorNevyazki(A, X_pred, b);
	tau = ScalarMult(r, r) / ScalarMult(MultMatrixVector(A, r), r);

	do
	{
		if (iter > 0) {

			r = vectorNevyazki(A, X, b);
			tau = ScalarMult(r, r) / ScalarMult(MultMatrixVector(A, r), r);
			alpha = 1 / (1 - (tau / tau_old) * (1 / alpha_old) * ScalarMult(r, r) /
				ScalarMult(vectorNevyazki(A, X_pred, b), vectorNevyazki(A, X_pred, b)));
		}


		X_now = SubtractionVector(AdditionVector(MultVectorNum(X, alpha), MultVectorNum(X_pred, 1 - alpha)), MultVectorNum(r, alpha * tau));


		q = (iter == 0) ? FirstVectorNorm(SubtractionVector(X, X_now)) : qCalc(X_pred, X, X_now);
		pg = Pogreshnost(X, X_now, q);
		norm_pg = FirstVectorNorm(SubtractionVector(X_now, X_q));


		cout << setw(4) << iter + 1 << "|" << fixed << setprecision(4) << tau << " | " << setprecision(3) << q
			<< " |" << setprecision(7) << setw(10) << FirstVectorNorm(r) << " |" << setw(10) << norm_pg << " |" << setw(10) << pg << " | ";


		for (int i = 0; i < b.size(); i++)
		{
			cout << setw(8) << fixed << setprecision(5) << X_now[i];
		}
		cout << fixed << setprecision(14) << " | " << alpha << endl;


		alpha_old = alpha;
		tau_old = tau;
		iter++;
		X_pred = X;
		X = X_now;

	} while (FirstVectorNorm(SubtractionVector(X_now, X_q)) / FirstVectorNorm(X_now) > epsylon);
}

void SORNext(vector<vector<double>> a, vector<double> b, double w, vector<double>& xCur, vector<double>& xNext)
{
	double sum;
	for(int i=0;i<b.size();++i)
	{
		sum = 0;

		for (int j = 0; j < i; ++j)
			sum += a[i][j] * xNext[j];
		for (int j = i + 1; j < b.size(); ++j)
			sum += a[i][j] * xCur[j];

		xNext[i] = xCur[i] + w * ((b[i] - sum) / a[i][i] - xCur[i]);
	}
}

double findOptimalW(vector<vector<double>> a, vector<double> b)
{
	double eps = 0.01, w, wOpt = 0.1,n;
	vector<double> nev(4);
	int itr, minItr = INT_MAX;
	for(w=0.1;w<2;w+=0.1)
	{
		itr = 0;
		n = 1;
		vector<double> xCur(4);
		vector<double> xNext(4);

		while(n>eps)
		{
			itr++;

			double sum;

			SORNext(a, b, w, xCur, xNext);

			nev = vectorNevyazki(a, xCur, b);
			n = FirstVectorNorm(nev);

			for (int i = 0; i < b.size(); ++i)
				xCur[i] = xNext[i];
		}

		if(itr<minItr)
		{
			wOpt = w;
			minItr = itr;
		}

		cout << "w= " << setprecision(2) << w << " Itr= " << itr << endl;
	}

	cout << "w* = " << setprecision(4) << wOpt << " ItrMin = " << minItr << endl;
	return wOpt;
}

void SOR(vector<vector<double>> a, vector<double> b, vector<double> solvedX)
{
	double w = findOptimalW(a, b);
	vector<double> xPrev(4);
	vector<double> xCur = b;
	vector<double> xNext(4);
	vector<double> nev;
	double nevNorm;
	double pogrNorm=1;

	double q, pogr=1;
	for(int itr=1;abs(pogrNorm)-epsylon>0;++itr)
	{
		nev = vectorNevyazki(a, xCur, b);
		SORNext(a, b, w, xCur, xNext);
		q = (itr == 1) ? FirstVectorNorm(SubtractionVector(xCur, xNext)) : qCalc(xPrev, xCur, xNext);
		pogr = Pogreshnost(xCur, xNext, q);
		nevNorm = FirstVectorNorm(nev);
		pogrNorm = FirstVectorNorm(SubtractionVector(xNext, solvedX));

		cout << itr << " | " << setprecision(3) << w <<setprecision(3)<< " | " << q << " | " << setprecision(7) << nevNorm << " | " << pogrNorm<<" | "<<pogr << endl;

		xPrev = xCur;
		xCur = xNext;
	}

}

int main()
{
	vector<double> b(4);
	b[0] = 1, b[1] = 2, b[2] = 3, b[3] = 4;

	ifstream fin("input.txt");

	vector<vector<double>> a(4, (vector<double>(4)));

	for (int i = 0; i < a.size(); ++i)
		for (int j = 0; j < a[i].size(); ++j)
			fin >> a[i][j];

	fin.close();

	vector<vector<double>> RevA = ReverseMatrix(a);
	double numberOb = EuclideanNorm(MultMatrix(Transpose(a), a)) * EuclideanNorm(MultMatrix(Transpose(RevA), RevA));

	cout << log(1 / epsylon) / 2 * numberOb << endl;
	cout << log(1 / epsylon) / 4 * sqrt(numberOb) << endl;
	cout << log(2 / epsylon) / 2 * sqrt(numberOb) << endl;

	for (int i = 0; i < a.size(); ++i)
	{
		for (int j = 0; j < a[i].size(); ++j)
			cout << a[i][j] << " ";
		cout << endl;
	}
	cout<<endl;

	
	vector<double> x = sqrtMethod(a, b);
	simpleIteration(a, b,x);
	cout << endl;
	fast_gardient(a, b, x);
	cout << endl;
	SOR(a, b, x);
	cout << endl;
	conjugate_gardient(a, b,x);

	for (int i = 0; i < x.size(); ++i)
		cout << x[i] << endl;

	

}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
