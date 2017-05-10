#ifndef SIGNAL
#define SIGNAL

#include <cassert>

namespace apl {

	class signal {

	public:
		double* data;
		int first;
		int last;
		unsigned int length;

		signal(const double* array, int firstIndex, int lastIndex);
		signal(double* array, int firstIndex, int lastIndex);
		signal(const signal& that);
		signal(signal&& that);
		signal& operator=(const signal& that);
		signal& operator=(signal&& that);
		inline double operator[](int index) const { int temp = index + first; assert(temp >= 0 && temp < (int)length); return data[temp]; }
		bool operator==(const signal& that) const;
		signal xcorr(const signal& that) const; // real signals only
		signal conv(const signal& that) const; // real signals only
		void print();
		~signal();

	};

};

#endif SIGNAL