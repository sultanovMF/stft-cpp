#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>

#include <eigen3/unsupported/Eigen/FFT>
#include <complex>
#include <iostream>
#include <vector>

#include "stft.hpp"

double hann_window(int n,int N){
	if(n<0||n>N)
		return 0;
	double s = sin(M_PI*n/(N-1));
	return s*s;
}
double fullspec_window(int n,int N){
	return 1;
}
wfunc_t window = hann_window;

CM stft(const std::vector<double> &signal,int nperseg,int overlap){
	
	//
    std::vector<double> padded;
	int window_step = nperseg - overlap;
	std::cout << "STFT"<<std::endl;
    int rpad = window_step * std::ceil(1.0*signal.size()/window_step) - signal.size();
    padded.reserve(rpad+signal.size()+window_step);
    for(int i=0;i<window_step;i++)
        padded.push_back(0);
    for(double s : signal)
        padded.push_back(s);
    for(int i=0;i<rpad;i++)
        padded.push_back(0);


	Eigen:Eigen::FFT<double>fft;
	fft.SetFlag(fft.HalfSpectrum);
	int T = ((padded.size()-nperseg)/window_step)+1;
	CM stftm(nperseg/2+1,T);
	for(int i=0;i<T;i++){
		std::vector<double> f(nperseg);
		for(int j=0;j<nperseg;j++){
			f[j] = window(j,nperseg) * padded[j+i*window_step];
		}
		std::vector<std::complex<double>> frq(nperseg/2+1,0);
		fft.fwd(frq,f);
		for(int j=0;j<nperseg/2+1;j++){
			stftm(j,i) = frq[j];
		}


	}

	std::cout << "STFT OK"<<std::endl;
	return stftm;
}


std::vector<double>  istft(const CM& spectrum,int nperseg,int overlap){
	std::cout << "ISTFT"<<std::endl;
	Eigen::FFT<double> fft;
	fft.SetFlag(fft.HalfSpectrum);
	int window_step = nperseg - overlap;
	int T = spectrum.cols();
	int length = (T)*window_step+nperseg;
	std::vector<double> signal(length,0);
	std::vector<double> norm(length,0);
	for(int i=0;i<T;i++){
		std::vector<std::complex<double>> frq(nperseg/2+1,0);
		std::vector<double> f(nperseg);
		for(int j = 0; j < nperseg/2+1;j++)
			frq[j] = spectrum(j,i);
		fft.inv(f,frq);

		for(int j = 0; j < nperseg;j++)
			if(window(j,nperseg)>1e-5){
				signal[j+window_step*i] += f[j]*window(j,nperseg);
                norm[j+window_step*i]+=window(j,nperseg)*window(j,nperseg);

			}

	}
	for(int i=0;i<length;i++)
        if(norm[i]>1e-10)
            signal[i]/=norm[i];


	return std::vector<double>(signal.begin()+window_step,signal.end());	
}

std::vector<double> get_time_samples(int frames_count, int overlap, double fs) {
    std::vector<double> time_samples (frames_count);

    for (int i = 0; i < frames_count; ++i) {
        time_samples[i] = i * overlap / fs;
    }

    return time_samples;
}
std::vector<double> get_frequency_samples(int freq_count,  double fs) {
    std::vector<double> frequency_samples (freq_count);

    for (int i = 0; i < freq_count; ++i) {
        frequency_samples[i] = i * fs / freq_count;
    }

    return frequency_samples;
}

PYBIND11_MODULE(koval_stft, m) {
  m.doc() = "Short time fourier transform";

  m.def("stft", &stft);
  m.def("istft", &istft);
  m.def("get_time_samples", &get_time_samples);
  m.def("get_frequency_samples", &get_frequency_samples);
}
