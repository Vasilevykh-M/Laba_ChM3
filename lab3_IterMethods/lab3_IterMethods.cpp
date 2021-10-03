// lab3_IterMethods.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include "MatrixOperation.h"
#include "VectorsOperation.h"

using namespace std;

double tau(vector<vector<double>> a)
{
	double y = 0.9;
	return 2 * y / EuclideanNorm(a);
}

double qCalc(vector<double> xPrev, vector<double> xCur, vector<double> xNext)
{
	return ThirdVectorNorm(SubtractionVector(xNext, xCur)) / ThirdVectorNorm(SubtractionVector(xCur,xPrev)); //Как тебе F#-way?
}

double Pogreshnost(vector<double> xCur, vector<double> xNext)
{
	vector<double> x = SubtractionVector(xNext, xCur);
	return min(ThirdVectorNorm(x), min(SecondVectorNorm(x), FirstVectorNorm(x))); //Мне очень нравится композиция функций в ФП
}

vector<double> vectorNevyazki(vector<vector<double>> a, vector<double> x, vector<double> b)
{
	vector<double> n(b.size());

	for (int i = 0; i < a.size(); ++i)
	{
		for (int j = 0; j < a[i].size(); ++j)
			n[i] += a[i][j] * x[i];
		n[i] -= b[i];
	}

	return n;
}

int main()
{

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
