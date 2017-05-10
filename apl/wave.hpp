#ifndef APL_WAVE
#define APL_WAVE

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

		wave(); //default constructor
		wave(const wave& that); // copy constructor
		wave(wave&& that); // move constructor
		wave(unsigned int sampleRate, unsigned short numChannels); // empty constructor
		wave(unsigned int sampleRate, unsigned short numChannels, unsigned int numSamples, double initValue); // constructor with scalar value
		wave(unsigned int sampleRate, unsigned short numChannels, unsigned int numSamples, const double* copy); // copy constructor of array copy
		wave(unsigned int sampleRate, unsigned short numChannels, unsigned int numSamples, double* data); // move constructor of array data
		wave(const std::string& filename); // read constructor
		bool read(const std::string& filename); // reads in wave with name string filename
		bool write(const std::string& filename, unsigned short bitDepth); // writes to filename with given bit depth
		inline double wave::get(int channel, int sample) const { return data[sample*numChannels + channel]; } // gets data from given channel and sample index
		inline double wave::operator[](int i) const { return data[i]; } // gives direct read access to data
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
		void append(const wave& that); // append that to the end of this
		static wave sine(double frequency, double amplitude, double phase, double duration, unsigned int sampleRate); // create sine wave
		~wave();

	};

}

#endif