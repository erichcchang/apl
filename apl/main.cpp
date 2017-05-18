
#include <iostream>
#include <string>
#include <vector>
#include "pitch.hpp"
#include "signal.hpp"
#include "wave.hpp"

using namespace std;
using namespace apl;




int main(int argc, char* argv[]) {
	// add user interface on console
	// add plotting library
	string testDir("C:/Users/erich/Downloads/audiotestfiles/");
	/*
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

	vector<chord> voice1;
	vector<chord> voice2;
	vector<chord> voice3;
	vector<chord> voice4;
	char chordnotes1[4] = { E5,Gb5,B5,E6 };
	char chordnotes2[4] = { D5,Gb5,A5,D6 };
	char chordnotes3[4] = { Gb4,Ab4,C5,Gb5 };
	char chordnotes4[4] = { E4,Gb4,B4,E5 };
	char chordnotes5[4] = { A3,D4,Gb4,A4 };
	char chordnotes6[2] = { Bb3,D4 };
	char chordnotes7[2] = { Db4,F4 };
	char chordnotes8[4] = { A2,D3,Gb3,B3 };
	char chordnotes9[2] = { Bb2,Ab3 };
	char chordnotes10[2] = { A2,G3 };
	char chordnotes11[2] = { D1,D2 };

	voice1.push_back(chord(chordnotes1, 4, 4));
	voice1.push_back(chord(chordnotes2, 4, 2));
	voice1.push_back(chord(chordnotes3, 4, 2));
	voice1.push_back(chord(nullptr, 0, 4));
	voice2.push_back(chord(nullptr, 0, 2));
	voice2.push_back(chord(chordnotes4, 0, 2));
	voice2.push_back(chord(chordnotes5, 4, 2));
	voice2.push_back(chord(nullptr, 0, 2));
	voice2.push_back(chord(chordnotes6, 2, 2));
	voice2.push_back(chord(chordnotes7, 2, 2));
	voice3.push_back(chord(nullptr, 0, 2));
	voice3.push_back(chord(chordnotes8, 4, 4));
	voice3.push_back(chord(chordnotes6, 2, 2));
	voice3.push_back(chord(chordnotes9, 2, 2));
	voice3.push_back(chord(chordnotes10, 2, 2));
	voice4.push_back(chord(chordnotes11, 2, 6));
	voice4.push_back(chord(nullptr, 0, 6));

	char chordnotes12[3] = { Db5,Gb5,Db6 };
	char chordnotes13[4] = { B4,D5,Gb5,B5 };
	char chordnotes14[4] = { B3,D4,Gb4,B4 };
	char chordnotes15[3] = { Ab3,D4,Gb4 };
	char chordnotes16[2] = { Ab2,Gb3 };

	voice1.push_back(chord(chordnotes1, 4, 4));
	voice1.push_back(chord(chordnotes2, 4, 2));
	voice1.push_back(chord(chordnotes12, 4, 4));
	voice1.push_back(chord(chordnotes13, 4, 2));
	voice2.push_back(chord(nullptr, 0, 2));
	voice2.push_back(chord(chordnotes4, 0, 2));
	voice2.push_back(chord(chordnotes5, 4, 2));
	voice2.push_back(chord(nullptr, 0, 2));
	voice2.push_back(chord(chordnotes14, 4, 4));
	voice3.push_back(chord(nullptr, 0, 2));
	voice3.push_back(chord(chordnotes8, 4, 4));
	voice3.push_back(chord(chordnotes15, 3, 4));
	voice3.push_back(chord(chordnotes15, 3, 2));
	voice4.push_back(chord(chordnotes11, 2, 6));
	voice4.push_back(chord(nullptr, 0, 2));
	voice4.push_back(chord(chordnotes16, 2, 4));


	wave right = createwave(voice1, 50, 2)*0.4 + createwave(voice2, 50, 2)*0.2;
	wave left = createwave(voice3, 50, 2)*0.2 + createwave(voice4, 50, 2)*0.2;
	wave together = left + right;
	together.write(testDir + "testsyntheticmusic.wav", 16);
	

	complex<double> arraya[4] = {1.0, 2.0, 3.0, 4.0};
	complex<double> arrayb[4] = {2, 3, 4, 1};
	complex<double> arrayarev[4] = {4, 3, 2, 1};
	complex<double> arraybrev[4] = {1, 4, 3, 2};
	signal a(arraya, 1, 4);
	signal b(arrayb, 1, 4);
	signal arev(arrayarev, -4, -1);
	signal brev(arraybrev, -4, -1);
	signal aconvb = a.conv(b);
	signal bconva = b.conv(a);
	if (aconvb == bconva) {
		cout << "conv test passed" << endl;
	}
	else {
		cout << "conv test failed" << endl;
	}
	signal acorrb = a.xcorr(b);
	signal arevconvb = arev.conv(b);
	if (acorrb == arevconvb) {
		cout << "xcorr test passed" << endl;
	}
	else {
		cout << "xcorr test failed" << endl;
	}
	signal bcorra = b.xcorr(a);
	signal brevconva = brev.conv(a);
	if (bcorra == brevconva) {
		cout << "xcorr test passed" << endl;
	}
	else {
		cout << "xcorr test failed" << endl;
	}
	cout.setf(ios::fixed, ios::floatfield);
	cout.precision(5);
	cout << aconvb << endl;
	cout << bconva << endl;
	cout << acorrb << endl;
	cout << arevconvb << endl;
	cout << bcorra << endl;
	cout << brevconva << endl;

	signal c(arraya, 0, 3);
	signal d(arrayb, 0, 3);
	signal cdft = c.dft();
	signal ddft = d.dft();
	cout << cdft << endl;
	cout << ddft << endl;
	signal ctr = cdft.idft();
	signal dtr = ddft.idft();
	cout << ctr << endl;
	cout << dtr << endl;
	if (ctr == c && dtr == d) {
		cout << "dft and idft test passed" << endl;
	}
	else {
		cout << "dft and idft test failed" << endl;
	}
	signal cfftcconvd = (cdft * ddft).idft();
	cout << cfftcconvd << endl;
	*/
	signal e(0, 7);
	signal f(3, 5);
	signal g(4, 10);
	signal h(-5, 20);
	signal i(-7, 2);
	signal j(10, 13);
	signal k(-9, -4);
	signal edotf = 2 - e * f + 1;
	signal eplusf = e + f;	
	signal edotg = e * g;
	signal eplusg = e + g;
	signal edoth = e * h;
	signal eplush = e + h;
	signal edoti = e * i;
	signal eplusi = e + i;
	signal edotj = e * j;
	signal eplusj = e + j;
	signal edotk = e * k;
	signal eplusk = e + k;
	
	cout << edotf << endl;
	cout << eplusf << endl;
	cout << edotg << endl;
	cout << eplusg << endl;
	cout << edoth << endl;
	cout << eplush << endl;
	cout << edoti << endl;
	cout << eplusi << endl;
	cout << edotj << endl;
	cout << eplusj << endl;
	cout << edotk << endl;
	cout << eplusk << endl;
}