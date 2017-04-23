

#include <string>
#include <iostream>
#include "wave.hpp"

using namespace std;
using namespace apl;

int main(int argc, char* argv[]) {
	// add user interface on console
	// add plotting library
	string testDir("C:/Users/erich/Downloads/audiotestfiles/");
	string filename1 = testDir + "sine440mono.wav";
	string filename2 = testDir + "chirp440to1320mono.wav";
	string filename3 = testDir + "sine440stereoright.wav";
	wave signal1(filename1);
	for (int i = 0; i < 102; i++) {
		cout << signal1.get(0, i) << endl;
	}
	cout << endl << endl;
	wave signal2(filename2);
	for (int i = 0; i < 400; i++) {
		cout << signal2.get(0, i) << endl;
	}
	cout << endl << endl;
	wave signal3(filename3);
	for (int i = 0; i < 102; i++) {
		cout << signal3.get(0, i) << endl;
	}
	cout << endl << endl;
	for (int i = 0; i < 102; i++) {
		cout << signal3.get(1, i) << endl;
	}
	cout << endl << endl;
	wave signal4 = (signal1 + signal2) / 2.0;
	string filename4 = testDir + "sinepluschirp.wav";
	if (signal4.write(filename4, 16)) {
		wave signal5(filename4);
		for (int i = 0; i < 400; i++) {
			cout << signal5.get(0, i) << endl;
		}
	}
}