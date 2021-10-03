// lab3_IterMethods.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include "MatrixOperation.h"
#include "VectorsOperation.h"
#include <iomanip>
#include <fstream>

using namespace std;

const double epsylon = 0.0001;

double tauCalc(vector<vector<double>> a)
{
	double y = 0.9;
	return y*(2 / EuclideanNorm(MultMatrix(Transpose(a),a)));
}

double qCalc(vector<double> xPrev, vector<double> xCur, vector<double> xNext)
{
	return FirstVectorNorm(SubtractionVector(xNext, xCur)) / FirstVectorNorm(SubtractionVector(xCur,xPrev)); //Как тебе F#-way?
}

double Pogreshnost(vector<double> xCur, vector<double> xNext)
{
	return  FirstVectorNorm(SubtractionVector(xNext, xCur)); //Мне очень нравится композиция функций в ФП
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



void simpleIteration(vector<vector<double>> a, vector<double> b)
{
	vector<double> xPrev(4,0);
	vector<double> xCur=b;
	vector<double> xNext(4);
	vector<double> nev;
	double nevNorm=1, tau=tauCalc(a), q=0, pogr=0;

	for(int i=1;nevNorm>epsylon;++i)
	{
		for(int j=0;j<a.size();++j)
		{
			double sum = 0;
			for (int k = 0; k < a[j].size(); ++k)
				sum += a[j][k] * xCur[k];
			xNext[j] = xCur[j] + tau * (b[j] - sum);
		}

		nev = vectorNevyazki(a, xCur, b);
		nevNorm = FirstVectorNorm(nev);
		q =(i==0)?FirstVectorNorm(xNext):qCalc(xPrev, xCur, xNext);
		pogr = Pogreshnost(xCur, xNext);
		cout << setprecision(4) << tau << " | " << q << " | " << setprecision(8) << nevNorm << " | " << pogr  << endl;

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

	for (int i = 0; i < a.size(); ++i)
	{
		for (int j = 0; j < a[i].size(); ++j)
			cout << a[i][j] << " ";
		cout << endl;
	}

	simpleIteration(a, b);

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
