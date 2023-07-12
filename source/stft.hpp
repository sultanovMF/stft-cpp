#pragma once

#include <eigen3/Eigen/Core>
#include <eigen3/unsupported/Eigen/FFT>
#include <complex>

typedef double (*wfunc_t)(int n,int N);
typedef Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> CM;

extern wfunc_t window;

CM stft(const std::vector<double> &signal,int nperseg = 256,int overlap = 128);
std::vector<double>  istft(const CM& spectrum,int nperseg = 256,int overlap = 128);
std::vector<double> get_time_samples(int frames_count, int overlap, double fs);
std::vector<double> get_frequency_samples(int nperseg,  double fs);