
#include <cassert>
#include <cmath>
#include "pitch.hpp"

using namespace apl;

chord::chord() {
}

chord::chord(const char* notes, unsigned char numNotes, unsigned char length) {
	this->notes = new char[numNotes];
	for (int i = 0; i < numNotes; ++i) {
		this->notes[i] = notes[i];
	}
	this->numNotes = numNotes;
	this->length = length;
}

chord::chord(const chord& that) {
	notes = new char[that.numNotes];
	for (int i = 0; i < that.numNotes; ++i) {
		notes[i] = that.notes[i];
	}
	numNotes = that.numNotes;
	length = that.length;
}

chord chord::transpose(int semitones) {
	chord transposed;
	transposed.numNotes = numNotes;
	transposed.length = length;
	transposed.notes = new char[numNotes];
	for (int i = 0; i < numNotes; ++i) {
		int temp = notes[i] + semitones;
		assert(temp >= -48 && temp < 40);
		transposed.notes[i] = temp;
	}
	return transposed;
}

wave chord::createwave(const std::vector<chord>& phrase, double bpm, int beatlength) {
	wave music(44100, 1);
	for (chord chord : phrase) {
		wave next;
		if (chord.numNotes > 0) {
			next = wave::sine(440.0*pow(2.0, chord.notes[0] / 12.0), 0.5 / chord.numNotes, 0.0, chord.length*bpm / (60 * beatlength), 44100);
			for (int i = 1; i < chord.numNotes; ++i) {
				next += wave::sine(440.0*pow(2.0, chord.notes[i] / 12.0), 0.5 / chord.numNotes, 0.0, chord.length*bpm / (60 * beatlength), 44100);
			}
		}
		else {
			next = wave(44100, 1, (unsigned int)ceil(44100 * chord.length*bpm / (60 * beatlength)), 0.0);
		}
		music.append(next);
	}
	return music;
}

chord::~chord() {
	delete[] notes;
}