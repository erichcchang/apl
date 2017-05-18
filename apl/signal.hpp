#ifndef SIGNAL
#define SIGNAL

#include <cassert>
#include <complex>
#include <iostream>

namespace apl {

	class signal {

	public:
		std::complex<double>* data;
		int first;
		int last;
		unsigned int length;

		signal(int firstIndex, int lastIndex);
		signal(const std::complex<double>* array, int firstIndex, int lastIndex);
		signal(int firstIndex, int lastIndex, std::complex<double>* array);
		signal(const signal& that);
		signal(signal&& that);
		signal& operator=(const signal& that);
		signal& operator=(signal&& that);
		inline std::complex<double> operator[](int index) const { int temp = index - first; assert(temp >= 0 && temp < (int)length); return data[temp]; }
		friend std::ostream& operator<<(std::ostream& os, const signal& that);

		bool operator==(const signal& that) const; // add support for zero padded signals
		bool operator!=(const signal& that) const;
		signal& operator+=(const signal& that);
		signal& operator+=(const std::complex<double>& constant);
		signal operator+(const signal& that) const;
		signal operator+(const std::complex<double>& constant) const;
		friend signal operator+(const std::complex<double>& constant, const signal& that);
		signal& operator-=(const signal& that);
		signal& operator-=(const std::complex<double>& constant);
		signal operator-(const signal& that) const;
		signal operator-(const std::complex<double>& constant) const;
		signal operator-() const;
		friend signal operator-(const std::complex<double>& constant, const signal& that);
		signal& operator*=(const signal& that);
		signal& operator*=(const std::complex<double>& constant);
		signal operator*(const signal& that) const;
		signal operator*(const std::complex<double>& constant) const;
		friend signal operator*(const std::complex<double>& constant, const signal& that);
		signal& operator/=(const signal& that);
		signal& operator/=(const std::complex<double>& constant);
		signal operator/(const signal& that) const;
		signal operator/(const std::complex<double>& constant) const;
		friend signal operator/(const std::complex<double>& constant, const signal& that);
		signal& operator^=(const std::complex<double>& constant);
		signal operator^(const std::complex<double>& constant) const;
		friend signal operator^(const std::complex<double>& constant, const signal& that);

		std::complex<double> signal::sum() const;
		double lpnorm(double p) const;
		signal cconj() const;
		signal xcorr(const signal& that) const;
		signal conv(const signal& that) const;
		// circconv
		// fftconv
		signal dft();
		signal idft();
		~signal();

	private: 
		int idftmode;

		signal dft(const signal& x);

	};

};

#endif SIGNAL