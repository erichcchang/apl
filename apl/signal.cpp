
#include <fstream>
#include <iostream>
#include "signal.hpp"

using namespace apl;

signal::signal(const double* array, int firstIndex, int lastIndex) {
	assert(lastIndex - firstIndex);
	first = firstIndex;
	last = lastIndex;
	length = lastIndex - firstIndex + 1;
	data = new double[length];
	for (unsigned int i = 0; i < length; ++i) {
		data[i] = array[i];
	}
}

signal::signal(double* array, int firstIndex, int lastIndex) {
	assert(lastIndex - firstIndex);
	first = firstIndex;
	last = lastIndex;
	length = lastIndex - firstIndex + 1;
	data = array;
}

signal::signal(const signal& that) {
	*this = that;
}

signal::signal(signal&& that) {
	*this = std::move(that);
}

signal& signal::operator=(const signal& that) {
	if (this != &that) {
		first = that.first;
		last = that.last;
		length = that.length;
		delete[] data;
		data = new double[length];
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

bool signal::operator==(const signal& that) const {
	if (length != that.length) {
		return false;
	}
	if (first != that.first) {
		return false;
	}
	for (unsigned int i = 0; i < length; ++i) {
		if (data[i] != that.data[i]) {
			return false;
		}
	}
	return true;
}

signal signal::xcorr(const signal& that) const {
	unsigned int ylength = length + that.length - 1;
	unsigned int yfirst = that.first - last;
	double* ydata = new double[ylength];
	if (length <= that.length) {
		for (unsigned int i = 0; i < ylength; ++i) {
			double sum = 0;
			for (unsigned int j = 0; j < length; ++j) {
				int bound = i + j - (length - 1);
				if (bound >= 0 && bound < (int)that.length) {
					sum += data[j] * that.data[bound];
				}
			}
			ydata[i] = sum;
		}
	}
	else {
		for (unsigned int i = 0; i < ylength; ++i) {
			double sum = 0;
			for (unsigned int j = 0; j < that.length; ++j) {
				int bound = j - i + (length - 1);
				if (bound >= 0 && bound < (int)length) {
					sum += data[bound] * that.data[j];
				}
			}
			ydata[i] = sum;
		}
	}
	return signal(ydata, yfirst, yfirst + ylength - 1);
}

signal signal::conv(const signal& that) const {
	unsigned int ylength = length + that.length - 1;
	unsigned int yfirst = first + that.first;
	double* ydata = new double[ylength];
	if (length <= that.length) {
		for (unsigned int i = 0; i < ylength; ++i) {
			double sum = 0;
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
		for (unsigned int i = 0; i < ylength; ++i) {
			double sum = 0;
			for (unsigned int j = 0; j < that.length; ++j) {
				int bound = i - j;
				if (bound >= 0 && bound < (int)length) {
					sum += that.data[j] * data[bound];
				}
			}
			ydata[i] = sum;
		}
	}

	return signal(ydata, yfirst, yfirst + ylength - 1);
}

void signal::print() {
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(5);
	std::cout << "[ ";
	for (unsigned int i = 0; i < length - 1; ++i) {
		std::cout << "(" << (int)i + first  << ", "  << data[i] << "), ";
	}
	std::cout << "(" << last << ", " << data[length-1] << ") ]" << std::endl;
}

signal::~signal() {
	delete[] data;
}