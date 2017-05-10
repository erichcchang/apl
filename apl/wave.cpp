#define _USE_MATH_DEFINES

#include <cassert>
#include <cmath> 
#include <fstream>
#include "wave.hpp"

using namespace apl;

wave::wave() {
	allocated = false;
}

wave::wave(const wave& that) {
	allocated = false;
	*this = that;
}

wave::wave(wave&& that) {
	allocated = false;
	*this = std::move(that);
}

wave::wave(unsigned int sampleRate, unsigned short numChannels) {
	this->sampleRate = sampleRate;
	this->numChannels = numChannels;
	this->bitDepth = 0;
	this->numSamples = 0;
	duration = 0;
	length = 0;
	allocated = false;
}

wave::wave(unsigned int sampleRate, unsigned short numChannels, unsigned int numSamples, double initValue) {
	this->sampleRate = sampleRate;
	this->numChannels = numChannels;
	this->bitDepth = 0;
	this->numSamples = numSamples;
	duration = numSamples / (double)sampleRate;
	length = numSamples*numChannels;
	data = new double[length];
	for (unsigned int i = 0; i < length; ++i) {
		data[i] = initValue;
	}
	allocated = true;
}

wave::wave(unsigned int sampleRate, unsigned short numChannels, unsigned int numSamples, const double* copy) {
	this->sampleRate = sampleRate;
	this->numChannels = numChannels;
	this->bitDepth = 0;
	this->numSamples = numSamples;
	duration = numSamples / (double)sampleRate;
	length = numSamples*numChannels;
	data = new double[length];
	for (unsigned int i = 0; i < length; ++i) {
		data[i] = copy[i];
	}
	allocated = true;
}

wave::wave(unsigned int sampleRate, unsigned short numChannels, unsigned int numSamples, double* data) {
	this->sampleRate = sampleRate;
	this->numChannels = numChannels;
	this->bitDepth = 0;
	this->numSamples = numSamples;
	duration = numSamples / (double)sampleRate;
	length = numSamples*numChannels;
	this->data = data;
	allocated = true;
}

wave::wave(const std::string& filename) {
	read(filename);
}

bool wave::read(const std::string& filename) {
	if (filename.substr(filename.length() - 4, 4) != std::string(".wav")) {
		return false;
	}
	std::fstream file(filename, std::fstream::in | std::fstream::binary);
	if (!file) {
		return false;
	}
	file.seekg(0, file.end);
	unsigned int len = (unsigned int)file.tellg();
	file.seekg(0, file.beg);
	char* temp = new char[len];
	file.read(temp, len);
	file.close();
	if (std::string(temp, temp + 4) != std::string("RIFF") && std::string(temp + 8, temp + 12) != std::string("WAVE")) {
		delete[] temp;
		return false;
	}
	if (std::string(temp + 12, temp + 16) != std::string("fmt ")) {
		delete[] temp;
		return false;
	}
	if (*((unsigned int*)(temp + 16)) != 16) {
		delete[] temp;
		return false;
	}
	if (*((unsigned short*)(temp + 20)) != 1) {
		delete[] temp;
		return false;
	}
	if (std::string(temp + 36, temp + 40) != std::string("data")) {
		delete[] temp;
		return false;
	}
	numChannels = *((unsigned short*)(temp + 22));
	sampleRate = *((unsigned int*)(temp + 24));
	bitDepth = *((unsigned short*)(temp + 34));
	unsigned int size = *((unsigned int*)(temp + 40));
	length = size * 8 / bitDepth;
	numSamples = length / numChannels;
	duration = numSamples / (double)sampleRate;
	char* start = temp + 44;
	if (bitDepth == 8) {
		data = new double[length];
		for (unsigned int i = 0; i < length; ++i) {
			data[i] = (((unsigned char*)start)[i] - 127) / 128.0;
		}
	}
	else if (bitDepth == 16) {
		data = new double[length];
		for (unsigned int i = 0; i < length; ++i) {
			data[i] = ((short*)start)[i] / 32768.0;
		}
	}
	else {
		delete[] temp;
		return false;
	}
	allocated = true;
	delete[] temp;
	return true;
}

bool wave::write(const std::string& filename, unsigned short bitDepth) {
	if (!allocated || (bitDepth != 8 && bitDepth != 16)) {
		return false;
	}
	std::fstream file(filename, std::fstream::out | std::fstream::binary);
	file.write("RIFF", 4);
	unsigned int size = 36 + length * bitDepth / 8;
	file.write((char*)&size, 4);
	file.write("WAVE", 4);
	file.write("fmt ", 4);
	size = 16;
	file.write((char*)&size, 4);
	unsigned short format = 1;
	file.write((char*)&format, 2);
	file.write((char*)&numChannels, 2);
	file.write((char*)&sampleRate, 4);
	size = sampleRate * numChannels * bitDepth / 8;
	file.write((char*)&size, 4);
	format = size/sampleRate;
	file.write((char*)&format, 2);
	file.write((char*)&bitDepth, 2);
	file.write("data", 4);
	size = format * numSamples;
	file.write((char*)&size, 4);
	char* temp = new char[size];
	if (bitDepth == 8) {
		unsigned char* start = (unsigned char*)temp;
		for (unsigned int i = 0; i < length; ++i) {
			double value = round(data[i] * 128.0 + 127.0);
			if (value > 255) {
				value = 255;
			}
			if (value < 0) {
				value = 0;
			}
			start[i] = (unsigned char)value;
		}
	}
	if (bitDepth == 16) {
		short* start = (short*)temp;
		for (unsigned int i = 0; i < length; ++i) {
			double value = round(data[i] * 32768.0);
			if (value > 32767) {
				value = 32767;
			}
			if (value < -32768) {
				value = -32768;
			}
			start[i] = (short)value;
		}
	}
	file.write(temp, size);
	file.close();
	return true;
}

wave& wave::operator=(const wave& that) {
	if (this != &that) {
		sampleRate = that.sampleRate;
		numChannels = that.numChannels;
		bitDepth = that.bitDepth;
		numSamples = that.numSamples;
		duration = that.duration;
		length = that.length;
		if (allocated) {
			delete[] data;
		}
		data = new double[length];
		for (unsigned int i = 0; i < length; ++i) {
			data[i] = that.data[i];
		}
		allocated = true;
	}
	return *this;
}

wave& wave::operator=(wave&& that) {
	if (this != &that) {
		sampleRate = that.sampleRate;
		numChannels = that.numChannels;
		bitDepth = that.bitDepth;
		numSamples = that.numSamples;
		duration = that.duration;
		length = that.length;
		if (allocated) {
			delete[] data;
		}
		data = that.data;
		allocated = that.allocated;
		that.data = nullptr;
	}
	return *this;
}

wave& wave::operator+=(const wave& that) {
	assert(sampleRate == that.sampleRate && numChannels == that.numChannels && numSamples == that.numSamples && allocated && that.allocated);
	if (bitDepth != that.bitDepth) {
		bitDepth = 0;
	}
	for (unsigned int i = 0; i < length; ++i) {
		data[i] += that.data[i];
	}
	return *this;
}

wave& wave::operator+=(double constant) {
	bitDepth = 0;
	for (unsigned int i = 0; i < length; ++i) {
		data[i] += constant;
	}
	return *this;
}

wave& wave::operator-=(const wave& that) {
	assert(sampleRate == that.sampleRate && numChannels == that.numChannels && numSamples == that.numSamples && allocated && that.allocated);
	if (bitDepth != that.bitDepth) {
		bitDepth = 0;
	}
	for (unsigned int i = 0; i < length; ++i) {
		data[i] -= that.data[i];
	}
	return *this;
}

wave& wave::operator-=(double constant) {
	bitDepth = 0;
	for (unsigned int i = 0; i < length; ++i) {
		data[i] -= constant;
	}
	return *this;
}

wave& wave::operator*=(const wave& that) {
	assert(sampleRate == that.sampleRate && numChannels == that.numChannels && numSamples == that.numSamples && allocated && that.allocated);
	if (bitDepth != that.bitDepth) {
		bitDepth = 0;
	}
	for (unsigned int i = 0; i < length; ++i) {
		data[i] *= that.data[i];
	}
	return *this;
}

wave& wave::operator*=(double constant) {
	bitDepth = 0;
	for (unsigned int i = 0; i < length; ++i) {
		data[i] *= constant;
	}
	return *this;
}

wave& wave::operator/=(const wave& that) {
	assert(sampleRate == that.sampleRate && numChannels == that.numChannels && numSamples == that.numSamples && allocated && that.allocated);
	if (bitDepth != that.bitDepth) {
		bitDepth = 0;
	}
	for (unsigned int i = 0; i < length; ++i) {
		data[i] /= that.data[i];
	}
	return *this;
}

wave& wave::operator/=(double constant) {
	bitDepth = 0;
	for (unsigned int i = 0; i < length; ++i) {
		data[i] /= constant;
	}
	return *this;
}

const wave wave::operator+(const wave& that) const {
	return wave(*this) += that;
}

const wave wave::operator+(double constant) const {
	return wave(*this) += constant;
}

const wave wave::operator-(const wave& that) const {
	return wave(*this) -= that;
}

const wave wave::operator-(double constant) const {
	return wave(*this) -= constant;
}

const wave wave::operator-() const{
	wave ans(*this);
	for (unsigned int i = 0; i < length; ++i) {
		ans.data[i] *= -1;
	}
	return ans;
}

const wave wave::operator*(const wave& that) const {
	return wave(*this) *= that;
}

const wave wave::operator*(double constant) const {
	return wave(*this) *= constant;
}

const wave wave::operator/(const wave& that) const {
	return wave(*this) /= that;
}

const wave wave::operator/(double constant) const {
	return wave(*this) /= constant;
}

void wave::append(const wave& that) {
	assert(sampleRate == that.sampleRate && numChannels == that.numChannels && that.allocated);
	unsigned int len = length + that.length;
	double* temp = new double[len];
	for (unsigned int i = 0; i < length; ++i) {
		temp[i] = data[i];
	}
	for (unsigned int i = length; i < len; ++i) {
		temp[i] = that.data[i-length];
	}
	length = len;
	numSamples += that.numSamples;
	duration += that.duration;
	data = temp;
	allocated = true;
}

wave wave::sine(double frequency, double amplitude, double phase, double duration, unsigned int sampleRate) {
	assert(2*frequency < sampleRate);
	int length = ceil(duration*sampleRate);
	double* x = new double[length];
	double w = 2*M_PI*frequency/sampleRate;
	for (unsigned int n = 0; n < length; ++n) {
		x[n] = amplitude*sin(w*n + phase);
	}
	return wave(sampleRate, 1, length, x);
}

wave::~wave() {
	if (allocated) {
		delete[] data;
	}
}