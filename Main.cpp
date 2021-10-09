#include <algorithm>
#include<vector>
#include <iomanip>
#include<iostream>
#include <fstream>

using namespace std;

const double eps = 0.0001;

std::vector<double> SubtractionVector(std::vector<double> a, std::vector<double> b)
{
    std::vector<double> res;
    for (int i = 0; i < a.size(); ++i)
        res.push_back(a[i] - b[i]);

    return res;
}

double FirstVectorNorm(std::vector<double> a)
{
    return abs(*std::max_element(a.begin(), a.end(), [](double x, double y) {return abs(x) < abs(y); }));
}

double SecondVectorNorm(std::vector<double> a)
{
    double sum = 0;

    for (int i = 0; i < a.size(); ++i)
        sum += std::abs(a[i]);

    return sum;
}

double ScalarMult(std::vector<double> a, std::vector<double> b)
{
    double sum = 0;

    for (int i = 0; i < a.size(); ++i)
        sum += a[i] * b[i];

    return sum;
}

double ThirdVectorNorm(std::vector<double> a)
{
    return ScalarMult(a, a);
}

std::vector<std::vector<double>> Transpose(std::vector<std::vector<double>> a)
{
    std::vector<std::vector<double>> res(a.size(), (std::vector<double>(a[0].size())));

    for (int i = 0; i < a.size(); ++i)
        for (int j = 0; j < a[i].size(); ++j)
            res[j][i] = a[i][j];

    return res;
}

std::vector<std::vector<double>> MultMatrix(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b)
{
    std::vector<std::vector<double>> res(a.size(), (std::vector<double>(a[0].size())));

    for (int i = 0; i < a.size(); ++i)
        for (int j = 0; j < a.size(); ++j)
        {
            double s = 0;

            for (int k = 0; k < a.size(); ++k)
                s += a[i][k] * b[k][j];

            res[i][j] = s;
        }

    return res;
}

std::vector<double> MultMatrixVector(std::vector<std::vector<double>> a, std::vector<double> b)
{
    std::vector<double> result(b.size());

    for (int i = 0; i < a.size(); ++i)
    {
        double s = 0;
        for (int j = 0; j < b.size(); ++j)
            s += a[i][j] * b[j];

        result[i] = s;
    }

    return result;
}

void search_max(int& i, int& j, std::vector<std::vector<double>> A)
{
    double max = LONG_MIN;

    for (int x = 0; x < A.size(); x++)
    {
        for (int y = 0; y < A[x].size(); y++)
        {
            if (abs(A[x][y]) > max && x != y)
            {
                max = abs(A[x][y]);
                i = x;
                j = y;
            }
        }
    }
}

double summ(std::vector<std::vector<double>> A)
{
    double sum = 0;
    for (int i = 0; i < A.size(); i++)
    {
        for (int j = 0; j < A.size(); j++)
        {
            if (i != j)
            {
                sum += pow(A[i][j], 2.0);
            }
        }
    }
    return sum;
}

double EuclideanNorm(std::vector<std::vector<double>> A)
{
    int x, y;

    std::vector<std::vector<double>> T(A.size(), std::vector<double>(A.size()));
    std::vector<std::vector<double>> TT;

    std::vector<std::vector<double>> AA;


    for (int i = 0; i < T.size(); i++)
    {
        for (int j = 0; j < T.size(); j++)
        {
            T[i][j] = 0;
            if (i == j) T[i][j] = 1;
        }
    }

    double alpha;

    while (abs(summ(A)) > eps) {

        for (int i = 0; i < T.size(); i++)
        {
            for (int j = 0; j < T.size(); j++)
            {
                T[i][j] = 0;
                if (i == j) T[i][j] = 1;
            }
        }

        search_max(x, y, A);
        double d = A[x][y];

        if (abs(A[x][x] - A[y][y]) < eps) {
            alpha = atan(1);
        }
        else
        {
            alpha = abs(atan(2 * A[x][y] / (A[x][x] - A[y][y]))) / 2;
        }

        T[x][x] = cos(alpha);
        T[x][y] = -sin(alpha);
        T[y][x] = sin(alpha);
        T[y][y] = cos(alpha);

        TT = Transpose(T);

        AA = MultMatrix(TT, A);

        A = MultMatrix(AA, T);
    }

    double min = LONG_MAX;
    double max = LONG_MIN;

    for (int i = 0; i < A.size(); i++)
    {
        if (A[i][i] > max)
            max = A[i][i];

        if (A[i][i] < min)
            min = A[i][i];
    }

    return sqrt(max);
}

std::vector<std::vector<double>> SubtractionMatrix(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b)
{
    std::vector<std::vector<double>> res(a.size(), (std::vector<double>(a.size())));

    for (int i = 0; i < res.size(); ++i)
        for (int j = 0; j < res.size(); ++j)
            res[i][j] = a[i][j] - b[i][j];

    return res;
}

double cubic_norm(std::vector<std::vector<double>> mass)
{
    double max = LONG_MIN;
    double summ = 0;

    for (int i = 0; i < mass.size(); i++)
    {
        for (int j = 0; j < mass[i].size(); j++)
        {
            summ += abs(mass[i][j]);
        }

        if (summ - max > eps)
            max = summ;

        summ = 0;
    }

    return max;
}

double octahedral_norm(std::vector<std::vector<double>> mass)
{
    double max = LONG_MIN;
    double summ = 0;

    for (int i = 0; i < mass[0].size(); i++)
    {
        for (int j = 0; j < mass.size(); j++)
        {
            summ += abs(mass[j][i]);
        }
        if (summ - max > eps) max = summ;
        summ = 0;
    }
    return max;
}

std::vector<double> AdditionVector(std::vector<double> A, std::vector<double> B) // я не помню есть ли это в твоей бибе или нет
{
    std::vector<double> res;
    for (int i = 0; i < A.size(); ++i)
        res.push_back(A[i] + B[i]);
    return res;
}

std::vector<double> MultVectorNum(std::vector<double> A, double B) // я не помню есть ли это в твоей бибе или нет
{
    std::vector<double> res;
    for (int i = 0; i < A.size(); ++i)
        res.push_back(A[i] * B);
    return res;
}

double tauCalc(vector<vector<double>> a)
{
    double y = 0.9;
    return y * (2 / EuclideanNorm(MultMatrix(Transpose(a), a)));
}

double qCalc(vector<double> xPrev, vector<double> xCur, vector<double> xNext)
{
    return FirstVectorNorm(SubtractionVector(xNext, xCur)) / FirstVectorNorm(SubtractionVector(xCur, xPrev)); //Как тебе F#-way?
}

double Pogreshnost(vector<double> xCur, vector<double> xNext, double q)
{
    return  q * FirstVectorNorm(SubtractionVector(xNext, xCur)) / (1 - q); //Мне очень нравится композиция функций в ФП
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
    vector<vector<double>> d(a.size(), (vector<double>(b.size())));


    for (int i = 0; i < S.size(); ++i)
    {
        double sum = 0;

        for (int k = 0; k < i; ++k)
            sum += S[k][i] * S[k][i] * d[k][k];

        d[i][i] = copysign(1.0, a[i][i] - sum);
        S[i][i] = sqrt(a[i][i] - sum);

        double temp = 1 / (S[i][i] * d[i][i]);

        for (int j = i + 1; j < S.size(); ++j)
        {
            sum = 0;

            for (int k = 0; k < i; ++k)
                sum += S[k][i] * S[k][j] * d[k][k];

            S[i][j] = (a[i][j] - sum) * temp;
        }

    }

    vector<double> y(b.size());
    vector<double> x(b.size());

    y[0] = b[0] / (S[0][0] * d[0][0]);

    for (int i = 1; i < b.size(); ++i)
    {
        double sum = 0;
        for (int j = 0; j < i; ++j)
            sum += S[j][i] * y[j] * d[j][j];

        y[i] = (b[i] - sum) / (S[i][i] * d[i][i]);
    }

    x[b.size() - 1] = y[b.size() - 1] / S[b.size() - 1][b.size() - 1];

    for (int i = x.size() - 2; i > -1; --i)
    {
        double sum = 0;
        for (int j = i; j < x.size(); ++j)
            sum += S[i][j] * x[j];

        x[i] = (y[i] - sum) / S[i][i];
    }

    return x;
}


void fast_gardient(vector<vector<double>> A, vector<double> b) //	Метод скорейшего спуска
{

    cout << "Метод наискорейшего спуска" << endl;
    cout << setw(4) << "Iter|" << setw(8) << "tau|" << setw(8) << "q|" << setw(12) 
        << "Норма r|" << setw(12) << "Норма погр|" << setw(12) << "Оценка погр|"
        << setw(8) << "X[0]" << setw(8) << "X[1]" << setw(8) << "X[2]" << setw(8) << "X[3]" << endl;


    vector<double> X(A.size()); // X[k+1]
    vector<double> X_pred(A.size()); // X[k]
    vector<double> X_prpr(A.size()); //X[k-1]
    vector<double> X_q(A.size());//X решение слау(X*)

    X_pred = b;
    X_prpr = b;

    X_q = sqrtMethod(A, b);

    int iter = 0;

    double t, q, pg, norm_pg, tau;

    vector<double> r(A.size()); // вектор невязки


    do {


        r = vectorNevyazki(A, X_pred, b);
        tau = ScalarMult(r, r) / ScalarMult(MultMatrixVector(A, r), r);


        X = AdditionVector(MultVectorNum(SubtractionVector(b,MultMatrixVector(A,X_pred)), tau), X_pred);


        q = (iter == 0) ? FirstVectorNorm(SubtractionVector(X_pred, X)) : qCalc(X_prpr, X_pred, X);
        pg = Pogreshnost(X_pred, X, q);
        norm_pg = FirstVectorNorm(SubtractionVector(X, X_q));

        X_prpr = X_pred;
        X_pred = X;

        cout << setw(4) << iter + 1 << "|" << fixed << setprecision(4) << tau << " | " << setprecision(3) << q
            <<" |" << setprecision(7) << setw(10) << FirstVectorNorm(r) << " |" << setw(10) << norm_pg << " |" << setw(10) << pg << " | ";

        for (int i = 0; i < 4; i++)
        {
            cout << setw(8) << fixed << setprecision(5) << X[i];
        }
        cout << endl;
        iter++;

    } while (FirstVectorNorm(SubtractionVector(X,X_q)) / FirstVectorNorm(X) > eps);
}

void conjugate_gardient(vector<vector<double>> A, vector<double> b) //метод сопряжённых градиентов
{

    cout << "Метод сопряженный градиетов" << endl;
    cout << setw(4) << "Iter|" << setw(8) << "tau|"  << setw(8) << "q|" <<setw(12)<<"Норма r|"
        <<setw(12)<<"Норма погр|"<<setw(12)<<"Оценка погр|" 
        <<setw(8)<<"X[0]"<< setw(8) << "X[1]" << setw(8) << "X[2]" << setw(8) << "X[3]" << endl;

    vector<double> X(A.size()); // X[k]
    vector<double> X_now(A.size()); //X[k+1]
    vector<double> X_pred(A.size()); // X[k-1]
    vector<double> X_q(A.size()); // X(*)

    vector<double> r(A.size()); // вектор невязки

    X_pred = b;
    X_now = b;
    X = b;
    X_q = sqrtMethod(A, b);

    int iter = 0;
    double tau, tau_old, alpha = 1, alpha_old = 1, pg,norm_pg, q;

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


        q = (iter == 0)? FirstVectorNorm(SubtractionVector(X, X_now)) : qCalc(X_pred, X, X_now);
        pg = Pogreshnost(X, X_now, q);
        norm_pg = FirstVectorNorm(SubtractionVector(X_now, X_q));

        
        cout <<setw(4)<<iter+1<<"|" << fixed << setprecision(4) << tau << " | " << setprecision(3) << q
            << " |" <<setprecision(7)<<setw(10) << FirstVectorNorm(r) << " |" << setw(10) << norm_pg << " |" << setw(10) << pg << " | ";


        for (int i = 0; i < 4; i++)
        {
            cout <<setw(8)<<fixed<< setprecision(5)<< X_now[i];
        }
        cout <<fixed<< setprecision(14)<<" | "<<alpha<< endl;


        alpha_old = alpha;
        tau_old = tau;
        iter++;
        X_pred = X;
        X = X_now;

    } while (FirstVectorNorm(SubtractionVector(X_now, X_q)) / FirstVectorNorm(X_now) > eps);
}



void simpleIteration(vector<vector<double>> a, vector<double> b)
{
    vector<double> xPrev(4, 0);
    vector<double> xCur = b;
    vector<double> xNext(4);
    vector<double> nev;
    double nevNorm = 1, tau = tauCalc(a), q = 0, pogr = 0;

    for (int i = 1; nevNorm > eps; ++i)
    {
        for (int j = 0; j < a.size(); ++j)
        {
            double sum = 0;
            for (int k = 0; k < a[j].size(); ++k)
                sum += a[j][k] * xCur[k];
            xNext[j] = xCur[j] + tau * (b[j] - sum);
        }

        nev = vectorNevyazki(a, xCur, b);
        nevNorm = FirstVectorNorm(nev);
        q = (i == 0) ? FirstVectorNorm(xNext) : qCalc(xPrev, xCur, xNext);
        pogr = Pogreshnost(xCur, xNext,q);
        cout << setprecision(4) << tau << " | " << q << " | " << setprecision(8) << nevNorm << " | " << pogr << endl;

        xPrev = xCur;
        xCur = xNext;
    }
}

int main()
{
    setlocale(LC_ALL, "rus");
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

    conjugate_gardient(a, b);

    fast_gardient(a, b);

}