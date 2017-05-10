#ifndef APL_PITCH
#define APL_PITCH

#include <vector>
#include "wave.hpp"

namespace apl {

	enum note {
		A0 = -48, Bb0 = -47, B0 = -46, C1 = -45, Db1 = -44, D1 = -43, Eb1 = -42, E1 = -41, F1 = -40, Gb1 = -39, G1 = -38, Ab1 = -37,
		A1 = -36, Bb1 = -35, B1 = -34, C2 = -33, Db2 = -32, D2 = -31, Eb2 = -30, E2 = -29, F2 = -28, Gb2 = -27, G2 = -26, Ab2 = -25,
		A2 = -24, Bb2 = -23, B2 = -22, C3 = -21, Db3 = -20, D3 = -19, Eb3 = -18, E3 = -17, F3 = -16, Gb3 = -15, G3 = -14, Ab3 = -13,
		A3 = -12, Bb3 = -11, B3 = -10, C4 = -9, Db4 = -8, D4 = -7, Eb4 = -6, E4 = -5, F4 = -4, Gb4 = -3, G4 = -2, Ab4 = -1,
		A4 = 0, Bb4 = 1, B4 = 2, C5 = 3, Db5 = 4, D5 = 5, Eb5 = 6, E5 = 7, F5 = 8, Gb5 = 9, G5 = 10, Ab5 = 11,
		A5 = 12, Bb5 = 13, B5 = 14, C6 = 15, Db6 = 16, D6 = 17, Eb6 = 18, E6 = 19, F6 = 20, Gb6 = 21, G6 = 22, Ab6 = 23,
		A6 = 24, Bb6 = 25, B6 = 26, C7 = 27, Db7 = 28, D7 = 29, Eb7 = 30, E7 = 31, F7 = 32, Gb7 = 33, G7 = 34, Ab7 = 35,
		A7 = 36, Bb7 = 37, B7 = 38, C8 = 39
	};

	// add dynamics
	class chord {

	public:
		char* notes;
		unsigned char numNotes;
		unsigned char length;

		chord();
		chord(const char* notes, unsigned char numNotes, unsigned char length);
		chord(const chord& that);
		chord transpose(int semitones);
		static wave createwave(const std::vector<chord>& phrase, double bpm, int beatlength);
		~chord();

	};

}

#endif