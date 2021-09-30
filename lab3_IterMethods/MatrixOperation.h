#pragma once
#include <vector>

std::vector<std::vector<double>> Transpose(std::vector<std::vector<double>>);
std::vector<std::vector<double>> MultMatrix(std::vector<std::vector<double>>, std::vector<std::vector<double>>);
std::vector<double> MultMatrixVector(std::vector<std::vector<double>>, std::vector<double>);
double EuclideanNorm(std::vector<std::vector<double>>);
