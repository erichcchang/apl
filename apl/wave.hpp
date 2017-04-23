#ifndef WAVE
#define WAVE

#include <cassert>
#include <fstream>
#include <string>

namespace apl {

	class wave {

	public:
		unsigned int sampleRate;
		unsigned short numChannels;
		unsigned short bitDepth;
		unsigned int numSamples;

		double duration;
		unsigned int length;
		double* data;
		bool allocated;

		wave();
		wave(const wave& that);
		wave(wave&& that);
		wave(unsigned int sampleRate, unsigned short numChannels, unsigned int numSamples);
		wave(unsigned int sampleRate, unsigned short numChannels, unsigned int numSamples, double initValue);
		wave(unsigned int sampleRate, unsigned short numChannels, unsigned int numSamples, const double* copy);
		wave(const std::string filename);
		bool read(const std::string filename);
		bool write(const std::string filename, unsigned short bitDepth);
		inline double wave::get(int channel, int sample) const { return data[sample*numChannels + channel]; }
		inline double wave::operator[](int i) const { return data[i]; }
		wave& operator=(const wave& that);
		wave& operator=(wave&& that);
		wave& operator+=(const wave& that);
		wave& operator+=(double constant);
		wave& operator-=(const wave& that);
		wave& operator-=(double constant);
		wave& operator*=(const wave& that);
		wave& operator*=(double constant);
		wave& operator/=(const wave& that);
		wave& operator/=(double constant);
		const wave operator+(const wave& that) const;
		const wave operator-(const wave& that) const;
		const wave operator+(double constant) const;
		const wave operator-(double constant) const;
		const wave operator-() const;
		const wave operator*(const wave& that) const;
		const wave operator/(const wave& that) const;
		const wave operator*(double constant) const;
		const wave operator/(double constant) const;
		~wave();

	};

}

#endif