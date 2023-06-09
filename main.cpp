#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>

// 補助関数
std::vector<double> A_array(double n) {
    double A0 = 1 + std::pow(n, 2) / 4.0 + std::pow(n, 4) / 64.0;
    double A1 = -1.5 * (n - std::pow(n, 3) / 8.0 - std::pow(n, 5) / 64.0);
    double A2 = 15.0 / 16.0 * (std::pow(n, 2) - std::pow(n, 4) / 4.0);
    double A3 = -35.0 / 48.0 * (std::pow(n, 3) - 5.0 / 16.0 * std::pow(n, 5));
    double A4 = 315.0 / 512.0 * std::pow(n, 4);
    double A5 = -693.0 / 1280.0 * std::pow(n, 5);
	 std::cout << "A0: " << A0 << std::endl;
	 std::cout << "A1: " << A1 << std::endl;
	 std::cout << "A2: " << A2 << std::endl;
	 std::cout << "A3: " << A3 << std::endl;
	 std::cout << "A4: " << A4 << std::endl;
	 std::cout << "A5: " << A5 << std::endl;
	 return {A0, A1, A2, A3, A4, A5};
}

std::vector<double> alpha_array(double n) {
    double a0 = std::nan(""); // dummy
    double a1 = 0.5 * n - 2.0 / 3.0 * std::pow(n, 2) + 5.0 / 16.0 * std::pow(n, 3) + 41.0 / 180.0 * std::pow(n, 4) - 127.0 / 288.0 * std::pow(n, 5);
    double a2 = 13.0 / 48.0 * std::pow(n, 2) - 3.0 / 5.0 * std::pow(n, 3) + 557.0 / 1440.0 * std::pow(n, 4) + 281.0 / 630.0 * std::pow(n, 5);
    double a3 = 61.0 / 240.0 * std::pow(n, 3) - 103.0 / 140.0 * std::pow(n, 4) + 15061.0 / 26880.0 * std::pow(n, 5);
    double a4 = 49561.0 / 161280.0 * std::pow(n, 4) - 179.0 / 168.0 * std::pow(n, 5);
    double a5 = 34729.0 / 80640.0 * std::pow(n, 5);
	 std::cout << "a0: " << a0 << std::endl;
	 std::cout << "a1: " << a1 << std::endl;
	 std::cout << "a2: " << a2 << std::endl;
	 std::cout << "a3: " << a3 << std::endl;
	 std::cout << "a4: " << a4 << std::endl;
	 std::cout << "a5: " << a5 << std::endl;
    return {a0, a1, a2, a3, a4, a5};
}

std::pair<double, double> calc_xy(double phi_deg, double lambda_deg, double phi0_deg, double lambda0_deg) {
    // 緯度経度・平面直角座標系原点をラジアンに直す
    double phi_rad = M_PI * phi_deg / 180.0;
    double lambda_rad = M_PI * lambda_deg / 180.0;
    double phi0_rad = M_PI * phi0_deg / 180.0;
    double lambda0_rad = M_PI * lambda0_deg / 180.0;
	 std::cout<< "phi_rad: " << phi_rad << std::endl;
	 std::cout<< "lambda_rad: " << lambda_rad << std::endl;
	 std::cout<< "phi0_rad: " << phi0_rad << std::endl;
	 std::cout<< "lambda0_rad: " << lambda0_rad << std::endl;
    // 定数 (a,F: 世界測地系-測地基準系1980（GRS80）楕円体)
    double m0 = 0.9999;
    double a = 6378137.0;
    double F = 298.257222101;

	 // sin_values ベクターを作成
	 std::vector<double> sin_values;
	 for (double i : {1, 2, 3, 4, 5}) {
		  sin_values.push_back(std::sin(2 * phi0_rad * i));
	 }


    // (1) n, A_i, alpha_iの計算
    double n = 1.0 / (2 * F - 1);
    std::vector<double> A_arr = A_array(n);
    std::vector<double> alpha_arr = alpha_array(n);

    // (2), S, Aの計算
    double A_ = (m0 * a) / (1.0 + n) * A_arr[0]; // [m]
	 double S_ = (m0 * a) / (1.0 + n) * (A_arr[0] * phi0_rad + std::inner_product(A_arr.begin() + 1, A_arr.end(), sin_values.begin(), 0.0)); // [m]
	 std::cout << "A_: " << A_ << std::endl;
	 std::cout <<"S_:"  <<	S_<< std::endl;

    // (3) lambda_c, lambda_sの計算
    double lambda_c = std::cos(lambda_rad - lambda0_rad);
    double lambda_s = std::sin(lambda_rad - lambda0_rad);
	 std::cout << "lambda_c: " << lambda_c << std::endl;
	 std::cout <<"lambda_s:" <<	lambda_s<<std::endl;

    // (4) t, t_の計算
    double t = std::sinh(std::atanh(std::sin(phi_rad)) - ((2 * std::sqrt(n)) / (1 + n)) * std::atanh(((2 * std::sqrt(n)) / (1 + n)) * std::sin(phi_rad)));
    double t_ = std::sqrt(1 + t * t);
	 std::cout << "t: " << t << std::endl;
	 std::cout <<"t_:" <<	t_<<std::endl;


    // (5) xi', eta'の計算
    double xi2 = std::atan(t / lambda_c); // [rad]
    double eta2 = std::atanh(lambda_s / t_);
	 std::cout << "xi2: " << xi2 << std::endl;
	 std::cout <<"eta2:" <<	eta2<<std::endl;

	 // cos_sin_values および cosh_sinh_values ベクターを作成
	 std::vector<double> cos_sin_values;
	 std::vector<double> cosh_sinh_values;
	 for (double i : {1, 2, 3, 4, 5}) {
		  cos_sin_values.push_back(std::cos(2 * xi2 * i));
		  cosh_sinh_values.push_back(std::cosh(2 * eta2 * i));
	 }


    // (6) x, yの計算
	 double x = A_ * (xi2 + std::inner_product(alpha_arr.begin() + 1, alpha_arr.end(),
															 std::vector<double>{1, 2, 3, 4, 5}.begin(), 0.0, std::plus<double>(),
															 [&](double alpha, double i) {
																  return alpha * std::sin(2 * xi2 * i) * std::cosh(2 * eta2 * i);
															 })) - S_; // [m]
	 std::cout << "x: " << x << std::endl;

	 double y = A_ * (eta2 + std::inner_product(alpha_arr.begin() + 1, alpha_arr.end(),
															  std::vector<double>{1, 2, 3, 4, 5}.begin(), 0.0, std::plus<double>(),
															  [&](double alpha, double i) {
																	return alpha * std::cos(2 * xi2 * i) * std::sinh(2 * eta2 * i);
															  }));

	 std::cout << "x: " << x << std::endl;
	 std::cout <<"y: " <<y<<std::endl;
    return std::make_pair(x, y); // [m]
}

int main() {
    double phi_deg = 36.103774791666666;      // 緯度
    double lambda_deg = 140.08785504166664;  // 経度
    double phi0_deg = 36.;     // 平面直角座標系の原点緯度
    double lambda0_deg = 139+50./60; // 原点経度


	 std::cout << "phi_deg: " << phi_deg << ", lambda_deg: " << lambda_deg << std::endl;
	 std::cout << "phi0_deg: " << phi0_deg << ", lambda0_deg: " << lambda0_deg << std::endl;
    std::pair<double, double> xy = calc_xy(phi_deg, lambda_deg, phi0_deg, lambda0_deg);
    std::cout << "x: " << xy.first << ", y: " << xy.second << std::endl;

    return 0;
  }
