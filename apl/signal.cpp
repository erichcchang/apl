#define _USE_MATH_DEFINES

#include <ostream>
#include <iostream>
#include "signal.hpp"

using namespace apl;
using namespace std;

signal::signal(int firstIndex, int lastIndex) {
	first = firstIndex;
	last = lastIndex;
	length = lastIndex - firstIndex + 1;
	assert(length);
	data = new complex<double>[length];
	for (int i = 0; i < (signed)length; ++i) {
		data[i] = first + i;
	}
}

signal::signal(const complex<double>* array, int firstIndex, int lastIndex) {
	first = firstIndex;
	last = lastIndex;
	length = lastIndex - firstIndex + 1;
	assert(length);
	data = new complex<double>[length];
	for (unsigned int i = 0; i < length; ++i) {
		data[i] = array[i];
	}
}

signal::signal(int firstIndex, int lastIndex, complex<double>* array) {
	first = firstIndex;
	last = lastIndex;
	length = lastIndex - firstIndex + 1;
	assert(length);
	data = array;
}

signal::signal(const signal& that) {
	first = that.first;
	last = that.last;
	length = that.length;
	data = new complex<double>[length];
	for (unsigned int i = 0; i < length; ++i) {
		data[i] = that.data[i];
	}
}

signal::signal(signal&& that) {
	first = that.first;
	last = that.last;
	length = that.length;
	data = that.data;
	that.data = nullptr;
}

signal& signal::operator=(const signal& that) {
	if (this != &that) {
		first = that.first;
		last = that.last;
		length = that.length;
		delete[] data;
		data = new complex<double>[length];
		for (unsigned int i = 0; i < length; ++i) {
			data[i] = that.data[i];
		}
	}
	return *this;
}

signal& signal::operator=(signal&& that) {
	if (this != &that) {
		first = that.first;
		last = that.last;
		length = that.length;
		delete[] data;
		data = that.data;
		that.data = nullptr;
	}
	return *this;
}

ostream& apl::operator<<(ostream& os, const signal& that) {
	os << "[ ";
	for (unsigned int i = 0; i < that.length - 1; ++i) {
		os << "(" << (int)i + that.first << ", " << that.data[i] << "), ";
	}
	os << "(" << that.last << ", " << that.data[that.length - 1] << ") ]";
	return os;
}

bool signal::operator==(const signal& that) const {
	if (length != that.length || first != that.first) {
		return false;
	}
	for (unsigned int i = 0; i < length; ++i) {
		if (abs(data[i] - that.data[i]) < 10*(-DBL_MAX)) {
			return false;
		}
	}
	return true;
}

bool signal::operator!=(const signal& that) const {
	return !(*this == that);
}

signal& signal::operator+=(const signal& that) {
	*this = *this + that;
	return *this;
}

signal& signal::operator+=(const std::complex<double>& constant) {
	for (unsigned int i = 0; i < length; ++i) {
		data[i] += constant;
	}
	return *this;
}

signal signal::operator+(const signal& that) const {
	if (last < that.first) {
		unsigned int ylength = that.last - first + 1;
		unsigned int offset = (unsigned int)(that.first - first);
		complex<double>* ydata = new complex<double>[ylength];
		for (unsigned int i = 0; i < length; ++i) {
			ydata[i] = data[i];
		}
		for (unsigned int i = length; i < offset; ++i) {
			ydata[i] = 0;
		}
		for (unsigned int i = 0; i < that.length; ++i) {
			ydata[i + offset] = that.data[i];
		}
		return signal(first, that.last, ydata);
	}
	else if (first > that.last) {
		unsigned int ylength = last - that.first + 1;
		unsigned int offset = (unsigned int)(first - that.first);
		complex<double>* ydata = new complex<double>[ylength];
		for (unsigned int i = 0; i < that.length; ++i) {
			ydata[i] = that.data[i];
		}
		for (unsigned int i = that.length; i < offset; ++i) {
			ydata[i] = 0;
		}
		for (unsigned int i = 0; i < length; ++i) {
			ydata[i + offset] = data[i];
		}
		return signal(that.first, last, ydata);
	}
	else if (first < that.first) {
		unsigned int offset = (unsigned int)(that.first - first);
		if (last < that.last) {
			complex<double>* ydata = new complex<double>[that.last - first + 1];
			for (unsigned int i = 0; i < offset; ++i) {
				ydata[i] = data[i];
			}
			for (unsigned int i = offset; i < length; ++i) {
				ydata[i] = data[i] + that.data[i - offset];
			}
			for (unsigned int i = last - that.first + 1; i < that.length; ++i) {
				ydata[i + offset] = that.data[i];
			}
			return signal(first, that.last, ydata);
		}
		else {
			complex<double>* ydata = new complex<double>[length];
			for (unsigned int i = 0; i < offset; ++i) {
				ydata[i] = data[i];
			}
			unsigned int endoffset = that.length + offset;
			for (unsigned int i = offset; i < endoffset; ++i) {
				ydata[i] = data[i] + that.data[i - offset];
			}
			for (unsigned int i = endoffset; i < length; ++i) {
				ydata[i] = data[i];
			}
			return signal(first, last, ydata);
		}
	}
	else {
		unsigned int offset = (unsigned int)(first - that.first);
		if (that.last < last) {
			complex<double>* ydata = new complex<double>[last - that.first + 1];
			for (unsigned int i = 0; i < offset; ++i) {
				ydata[i] = that.data[i];
			}
			for (unsigned int i = offset; i < that.length; ++i) {
				ydata[i] = that.data[i] + data[i - offset];
			}
			for (unsigned int i = that.last - first + 1; i < length; ++i) {
				ydata[i + offset] = data[i];
			}
			return signal(that.first, last, ydata);
		}
		else {
			complex<double>* ydata = new complex<double>[that.length];
			for (unsigned int i = 0; i < offset; ++i) {
				ydata[i] = that.data[i];
			}
			unsigned int endoffset = length + offset;
			for (unsigned int i = offset; i < endoffset; ++i) {
				ydata[i] = that.data[i] + data[i - offset];
			}
			for (unsigned int i = endoffset; i < that.length; ++i) {
				ydata[i] = that.data[i];
			}
			return signal(that.first, that.last, ydata);
		}
	}
}

signal signal::operator+(const std::complex<double>& constant) const {
	complex<double>* ydata = new complex<double>[length];
	for (unsigned int i = 0; i < length; ++i) {
		ydata[i] = data[i] + constant;
	}
	return signal(first, last, ydata);
}

signal apl::operator+(const std::complex<double>& constant, const signal& that) {
	return that + constant;
}

signal& signal::operator-=(const signal& that) {
	*this = *this - that;
	return *this;
}

signal& signal::operator-=(const std::complex<double>& constant) {
	for (unsigned int i = 0; i < length; ++i) {
		data[i] -= constant;
	}
	return *this;
}

signal signal::operator-(const signal& that) const {
	if (last < that.first) {
		unsigned int ylength = that.last - first + 1;
		unsigned int offset = (unsigned int)(that.first - first);
		complex<double>* ydata = new complex<double>[ylength];
		for (unsigned int i = 0; i < length; ++i) {
			ydata[i] = data[i];
		}
		for (unsigned int i = length; i < offset; ++i) {
			ydata[i] = 0;
		}
		for (unsigned int i = 0; i < that.length; ++i) {
			ydata[i + offset] = -that.data[i];
		}
		return signal(first, that.last, ydata);
	}
	else if (first > that.last) {
		unsigned int ylength = last - that.first + 1;
		unsigned int offset = (unsigned int)(first - that.first);
		complex<double>* ydata = new complex<double>[ylength];
		for (unsigned int i = 0; i < that.length; ++i) {
			ydata[i] = -that.data[i];
		}
		for (unsigned int i = that.length; i < offset; ++i) {
			ydata[i] = 0;
		}
		for (unsigned int i = 0; i < length; ++i) {
			ydata[i + offset] = data[i];
		}
		return signal(that.first, last, ydata);
	}
	else if (first < that.first) {
		unsigned int offset = (unsigned int)(that.first - first);
		if (last < that.last) {
			complex<double>* ydata = new complex<double>[that.last - first + 1];
			for (unsigned int i = 0; i < offset; ++i) {
				ydata[i] = data[i];
			}
			for (unsigned int i = offset; i < length; ++i) {
				ydata[i] = data[i] - that.data[i - offset];
			}
			for (unsigned int i = last - that.first + 1; i < that.length; ++i) {
				ydata[i + offset] = -that.data[i];
			}
			return signal(first, that.last, ydata);
		}
		else {
			complex<double>* ydata = new complex<double>[length];
			for (unsigned int i = 0; i < offset; ++i) {
				ydata[i] = data[i];
			}
			unsigned int endoffset = that.length + offset;
			for (unsigned int i = offset; i < endoffset; ++i) {
				ydata[i] = data[i] - that.data[i - offset];
			}
			for (unsigned int i = that.length; i < length; ++i) {
				ydata[i + offset] = data[i];
			}
			return signal(first, last, ydata);
		}
	}
	else {
		unsigned int offset = (unsigned int)(first - that.first);
		if (that.last < last) {
			complex<double>* ydata = new complex<double>[last - that.first + 1];
			for (unsigned int i = 0; i < offset; ++i) {
				ydata[i] = -that.data[i];
			}
			for (unsigned int i = offset; i < that.length; ++i) {
				ydata[i] = data[i - offset] - that.data[i];
			}
			for (unsigned int i = that.last - first + 1; i < length; ++i) {
				ydata[i + offset] = data[i];
			}
			return signal(that.first, last, ydata);
		}
		else {
			complex<double>* ydata = new complex<double>[that.length];
			for (unsigned int i = 0; i < offset; ++i) {
				ydata[i] = -that.data[i];
			}
			unsigned int endoffset = length + offset;
			for (unsigned int i = offset; i < endoffset; ++i) {
				ydata[i] = data[i - offset] - that.data[i];
			}
			for (unsigned int i = length; i < that.length; ++i) {
				ydata[i + offset] = -that.data[i];
			}
			return signal(that.first, that.last, ydata);
		}
	}
}

signal signal:: operator-(const std::complex<double>& constant) const {
	complex<double>* ydata = new complex<double>[length];
	for (unsigned int i = 0; i < length; ++i) {
		ydata[i] = data[i] - constant;
	}
	return signal(first, last, ydata);
}

signal signal::operator-() const {
	complex<double>* ydata = new complex<double>[length];
	for (unsigned int i = 0; i < length; ++i) {
		ydata[i] = -data[i];
	}
	return signal(first, last, ydata);
}

signal apl::operator-(const std::complex<double>& constant, const signal& that) {
	complex<double>* ydata = new complex<double>[that.length];
	for (unsigned int i = 0; i < that.length; ++i) {
		ydata[i] = constant - that.data[i];
	}
	return signal(that.first, that.last, ydata);
}

signal& signal::operator*=(const signal& that) {
	*this = *this * that;
	return *this;
}

signal& signal::operator*=(const complex<double>& constant) {
	for (unsigned int i = 0; i < length; ++i) {
		data[i] *= constant;
	}
	return *this;
}

signal signal::operator*(const signal& that) const {
	if (last < that.first || first > that.last) {
		complex<double>* ydata = new complex<double>[1];
		ydata[0] = 0;
		return signal(0, 0, ydata);
	}
	else if (first < that.first) {
		int ylast = last < that.last ? last : that.last;
		unsigned int ylength = ylast - that.first + 1;
		int offset = that.first - first;
		complex<double>* ydata = new complex<double>[ylength];
		for (unsigned int i = 0; i < ylength; ++i) {
			ydata[i] = data[i + offset] * that.data[i];
		}
		return signal(that.first, ylast, ydata);
	}
	else {
		int ylast = last < that.last ? last : that.last;
		unsigned int ylength = ylast - first + 1;
		int offset = first - that.first;
		complex<double>* ydata = new complex<double>[ylength];
		for (unsigned int i = 0; i < ylength; ++i) {
			ydata[i] = data[i] * that.data[i + offset];
		}
		return signal(first, ylast, ydata);
	}
}

signal signal::operator*(const complex<double>& constant) const {
	complex<double>* ydata = new complex<double>[length];
	for (unsigned int i = 0; i < length; ++i) {
		ydata[i] = constant * data[i];
	}
	return signal(first, last, ydata);
}

signal apl::operator*(const std::complex<double>& constant, const signal& that) {
	return that * constant;
}

signal& signal::operator/=(const signal& that) {
	*this = *this / that;
	return *this;
}

signal& signal::operator/=(const std::complex<double>& constant) {
	for (unsigned int i = 0; i < length; ++i) {
		data[i] /= constant;
	}
	return *this;
}

signal signal::operator/(const signal& that) const {
	if (last < that.first || first > that.last) {
		complex<double>* ydata = new complex<double>[1];
		ydata[0] = 0;
		return signal(0, 0, ydata);
	}
	else if (first <= that.first) {
		int ylast = last <= that.last ? last : that.last;
		unsigned int ylength = ylast - that.first + 1;
		int offset = that.first - first;
		complex<double>* ydata = new complex<double>[ylength];
		for (unsigned int i = 0; i < ylength; ++i) {
			ydata[i] = data[i + offset] / that.data[i];
		}
		return signal(that.first, ylast, ydata);
	}
	else {
		int ylast = last <= that.last ? last : that.last;
		unsigned int ylength = ylast - first + 1;
		int offset = first - that.first;
		complex<double>* ydata = new complex<double>[ylength];
		for (unsigned int i = 0; i < ylength; ++i) {
			ydata[i] = data[i] / that.data[i + offset];
		}
		return signal(first, ylast, ydata);
	}
}

signal signal::operator/(const std::complex<double>& constant) const {
	complex<double>* ydata = new complex<double>[length];
	for (unsigned int i = 0; i < length; ++i) {
		ydata[i] = data[i] / constant;
	}
	return signal(first, last, ydata);
}

signal apl::operator/(const std::complex<double>& constant, const signal& that) {
	complex<double>* ydata = new complex<double>[that.length];
	for (unsigned int i = 0; i < that.length; ++i) {
		ydata[i] = constant / that.data[i];
	}
	return signal(that.first, that.last, ydata);
}

signal& signal::operator^=(const std::complex<double>& constant) {
	*this = *this ^ constant;
	return *this;
}

signal signal::operator^(const std::complex<double>& constant) const {
	complex<double>* ydata = new complex<double>[length];
	for (unsigned int i = 0; i < length; ++i) {
		ydata[i] = pow(data[i], constant);
	}
	return signal(first, last, ydata);
}

signal apl::operator^(const std::complex<double>& constant, const signal& that) {
	complex<double>* ydata = new complex<double>[that.length];
	for (unsigned int i = 0; i < that.length; ++i) {
		ydata[i] = pow(constant, that.data[i]);
	}
	return signal(that.first, that.last, ydata);
}

complex<double> signal::sum() const {
	complex<double> sum = 0;
	for (unsigned int i = 0; i < length; ++i) {
		sum += data[i];
	}
	return sum;
}

double signal::lpnorm(double p) const {
	double sum = 0;
	for (unsigned int i = 0; i < length; ++i) {
		sum += pow(abs(data[i]), p);
	}
	return pow(sum, 1 / p);
}

signal signal::cconj() const {
	complex<double>* ydata = new complex<double>[length];
	for (unsigned int i = 0; i < length; ++i) {
		ydata[i] = conj(data[i]);
	}
	return signal(first, last, ydata);
}

signal signal::xcorr(const signal& that) const {
	unsigned int ylength = length + that.length - 1;
	unsigned int yfirst = that.first - last;
	complex<double>* ydata = new complex<double>[ylength];
	if (length <= that.length) {
		for (unsigned int i = 0; i < ylength; ++i) {
			complex<double> sum = 0;
			for (unsigned int j = 0; j < length; ++j) {
				int bound = i + j - (length - 1);
				if (bound >= 0 && bound < (int)that.length) {
					sum += conj(data[j]) * that.data[bound];
				}
			}
			ydata[i] = sum;
		}
	}
	else {
		for (unsigned int i = 0; i < ylength; ++i) {
			complex<double> sum = 0;
			for (unsigned int j = 0; j < that.length; ++j) {
				int bound = j - i + (length - 1);
				if (bound >= 0 && bound < (int)length) {
					sum += conj(data[bound]) * that.data[j];
				}
			}
			ydata[i] = sum;
		}
	}
	return signal(yfirst, yfirst + ylength - 1, ydata);
}

signal signal::conv(const signal& that) const {
	unsigned int ylength = length + that.length - 1;
	unsigned int yfirst = first + that.first;
	complex<double>* ydata = new complex<double>[ylength];
	if (length <= that.length) {
		#pragma omp parallel for
		for (unsigned int i = 0; i < ylength; ++i) {
			complex<double> sum = 0;
			#pragma omp simd
			for (unsigned int j = 0; j < length; ++j) {
				int bound = i - j;
				if (bound >= 0 && bound < (int)that.length) {
					sum += data[j] * that.data[bound];
				}
			}
			ydata[i] = sum;
		}
	}
	else {
		#pragma omp parallel for
		for (unsigned int i = 0; i < ylength; ++i) {
			complex<double> sum = 0;
			#pragma omp simd
			for (unsigned int j = 0; j < that.length; ++j) {
				int bound = i - j;
				if (bound >= 0 && bound < (int)length) {
					sum += that.data[j] * data[bound];
				}
			}
			ydata[i] = sum;
		}
	}

	return signal(yfirst, yfirst + ylength - 1, ydata);
}

signal signal::dft() {
	assert(first == 0);
	unsigned int padlength = 1;
	while (padlength < length) {
		padlength *= 2;
	}
	complex<double>* paddata = new complex<double>[padlength];
	for (unsigned int i = 0; i < length; ++i) {
		paddata[i] = data[i];
	}
	for (unsigned int i = length; i < padlength - 1; ++i) {
		paddata[i] = 0;
	}
	idftmode = -1;
	return dft(signal(0, padlength - 1, paddata));
}

signal signal::idft() {
	assert(first == 0);
	unsigned int padlength = 1;
	while (padlength < length) {
		padlength *= 2;
	}
	complex<double>* paddata = new complex<double>[padlength];
	for (unsigned int i = 0; i < length; ++i) {
		paddata[i] = data[i] / (double)padlength;
	}
	for (unsigned int i = length; i < padlength - 1; ++i) {
		paddata[i] = 0;
	}
	idftmode = 1;
	return dft(signal(0, padlength - 1, paddata));
}

signal signal::dft(const signal& x) {
	 if (x.length == 1) {
		 return signal(x);
	 }
	 else if (x.length == 2) {
		 complex<double>* ydata = new complex<double>[2];
		 ydata[0] = x.data[0] + x.data[1];
		 ydata[1] = x.data[0] - x.data[1];
		 return signal(0, 1, ydata);
	 }
	 else {
		 unsigned int halflength = x.length / 2;
		 complex<double>* evendata = new complex<double>[halflength];
		 complex<double>* odddata = new complex<double>[halflength];
		 for (unsigned int i = 0; i < x.length/2; ++i) {
			 evendata[i] = x.data[2*i];
			 odddata[i] = x.data[2*i + 1];
		 }
		 signal evendft = dft(signal(0, halflength - 1, evendata));
		 signal odddft = dft(signal(0, halflength - 1, odddata));
		 complex<double>* ydata = new complex<double>[x.length];
		 double w = 2*M_PI/x.length;
		 for (unsigned int i = 0; i < halflength; ++i) {
			 complex<double> c(cos(w*i), idftmode*sin(w*i));
			 ydata[i] = evendft.data[i] + c * odddft.data[i];
			 ydata[i + halflength] = evendft.data[i] - c * odddft.data[i];
		 }
		 return signal(0, x.length - 1, ydata);
	 }	 
}

signal::~signal() {
	delete[] data;
}