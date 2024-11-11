#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <algorithm>

using namespace std;

double f(double x) {
    return exp(x) - 1.0 / x;
}

double f_prime(double x) {
    return exp(x) + 1 / (x * x);
}

double f_double_prime(double x) {
    return exp(x) - 2 / (x * x * x);
}

vector<vector<double>> divided_differences(const vector<double>& x, const vector<double>& y) {
    int n = x.size();
    vector<vector<double>> a(n, vector<double>(n));

    for (int i = 0; i < n; ++i) {
        a[i][0] = y[i];
    }

    for (int j = 1; j < n; ++j) {
        for (int i = 0; i < n - j; ++i) {
            a[i][j] = (a[i + 1][j - 1] - a[i][j - 1]) / (x[i + j] - x[i]);
        }
    }

    return a;
}

double newton_derivative(const vector<double>& x, const vector<vector<double>>& a, double t, int order) {
    int n = x.size();
    double result = 0;

    if (order == 1) {
        for (int i = 1; i < n; ++i) {
            double sum = 0;
            for (int j = 0; j < i; ++j) {
                double prod = 1;
                for (int k = 0; k < i; ++k) {
                    if (k != j) {
                        prod *= (t - x[k]);
                    }
                }
                sum += prod;
            }
            result += a[0][i] * sum;
        }
    }
    else if (order == 2) {
        for (int i = 2; i < n; ++i) {
            double sum = 0;
            for (int j = 0; j < i; ++j) {
                for (int k = 0; k < i; ++k) {
                    if (k != j) {
                        double prod = 1;
                        for (int l = 0; l < i; ++l) {
                            if (l != j && l != k) {
                                prod *= (t - x[l]);
                            }
                        }
                        sum += prod;
                    }
                }
            }
            result += a[0][i] * sum;
        }
    }

    return result;
}

vector<int> nearest_points(const vector<double>& arr, double value, int n) {
    vector<pair<double, int>> diff(arr.size());
    for (int i = 0; i < arr.size(); ++i) {
        diff[i] = { abs(arr[i] - value), i };
    }
    sort(diff.begin(), diff.end());

    vector<int> result(n);
    for (int i = 0; i < n; ++i) {
        result[i] = diff[i].second;
    }
    return result;
}

int main() {
    setlocale(LC_ALL, "");

    int n;
    cout << "Enter the polynomial degree: ";
    cin >> n;

    vector<double> x(11), y(11), newX(21), newY_prime(21), newY_double_prime(21);

    for (int i = 0; i <= 10; ++i) {
        x[i] = 0.2 + 0.5 * i;
        y[i] = f(x[i]);
    }

    for (int j = 0; j < 21; ++j) {
        newX[j] = 0.2 + 0.25 * j;

        auto indices = nearest_points(x, newX[j], n + 1);
        vector<double> subX(n + 1), subY(n + 1);

        for (int i = 0; i <= n; ++i) {
            subX[i] = x[indices[i]];
            subY[i] = y[indices[i]];
        }

        auto subA = divided_differences(subX, subY);

        newY_prime[j] = newton_derivative(subX, subA, newX[j], 1);
        newY_double_prime[j] = newton_derivative(subX, subA, newX[j], 2);
    }

    printf("%5s %27s %16s %27s %17s\n", "X", "First derivative", "Error", "Second derivative", "Error");

    for (int j = 0; j < 21; ++j) {
        cout << std::defaultfloat
            << setw(0) << "X[" << j << "] = " << newX[j]
            << setw(17) << newY_prime[j]
            << setw(23) << std::scientific << abs(newY_prime[j] - f_prime(newX[j]))
            << setw(23) << std::defaultfloat << newY_double_prime[j]
            << setw(23) << std::scientific << abs(newY_double_prime[j] - f_double_prime(newX[j]))
            << "\n";
    }

    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout << "Press Enter to exit...";
    std::cin.get();


    return 0;
}
