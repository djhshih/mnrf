/**
 * Mononucleotide repeat finder
 *
 * Author:   David JH Shih <djh.shih@gmail.com>
 * Date:     2013-11-18
 * License:  GPLv3
 *
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cstring>

using namespace std;


void mnrf(unsigned period, istream& in, ofstream& out) {

	const char delim = '\t';
	const size_t buffer_size = 1024;
	char buffer[buffer_size];
	char chromosome[8];

	// assume only one track

	/// get track name

	in.getline(buffer, buffer_size);
	if (buffer[0] != '>') {
		throw runtime_error("Error: FASTA must have a track name");
	}

	// replace first space with null
	for (size_t i = 1; i < buffer_size; ++i) {
		if (buffer[i] == ' ') {
			buffer[i] = '\0';
			break;
		}
	}
	// copy to chromosome
	strcpy(chromosome, &(buffer[1]));


	size_t start = 0, pos = 0, count = 0;
	char prev_nuc = '\0';

	
	while (!in.eof()) {

		// get next line
		in.getline(buffer, buffer_size);
		if (buffer[0] == '\0') break;

		for (size_t i = 0; i < buffer_size; ++i) {

			if (buffer[i] == '\0') break;

			if (buffer[i] == prev_nuc) {
				++count;
			} else {

				if (count >= period) {
					// record current mononucleotide stretch
					out
						<< chromosome << delim
						<< start << delim
						<< (start + count) << delim
						<< prev_nuc << delim
						<< count << endl;
				}

				// prepare for next mononucleotide stretch
				prev_nuc = buffer[i];
				start = pos;
				count = 1;
			}

			++pos;

		}

	}

	// clean up last stretch
	if (count >= period) {
		// record current mononucleotide stretch
		out << chromosome << delim << start << delim << (start + count) << delim << prev_nuc << endl;
	}

}


int main(int argc, char* argv[]) {

	if (argc < 4) {
		cout << "Usage: mnrf <period> <input.fa> <output.bed>" << endl;
		return 0;
	}

	ifstream infile;
	ofstream outfile;

	int period = atoi(argv[1]);
	infile.open(argv[2]);
	outfile.open(argv[3]);

	if (!infile.is_open()) {
		throw runtime_error("Error: could not open input file");
	}

	if (!outfile.is_open()) {
		throw runtime_error("Error: could not open output file");
	}

	if (period < 2) {
		throw domain_error("Error: period must be at least 2");
	}

	mnrf(period, infile, outfile);

	infile.close();
	outfile.close();

	return 0;
}

