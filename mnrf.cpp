/*
    Mononucleotide repeat finder.
    Copyright (C) 2013-2014  David J. H. Shih <djh.shih@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cstring>

using namespace std;


const size_t buffer_size = 1024;


bool get_fasta_header(const char* line, char* header) {

	// header line begins with '>'
	// otherwise, current line is not a headerline
	if (line[0] != '>') {
		return false;
	}

	// find first space character
	size_t count = 0;
	for (size_t i = 1; i < buffer_size; ++i) {
		if (line[i] == ' ') break;
		++count;
	}
	// copy to header
	memcpy(header, &(line[1]), count);
	header[count] = '\0';

	return true;
}


void mnrf(unsigned nrepeats, istream& in, ofstream& out) {

	const char delim = '\t';
	const unsigned track_name_size = 32;
	const char ignore_nuc = 'N';

	char track_name[track_name_size];
	char buffer[buffer_size];

	size_t count = 1;
	size_t start = 0, pos = 0;
	char prev_nuc = '\0';

	track_name[0] = '\0';

	while (!in.eof()) {

		in.getline(buffer, buffer_size);

		// ignore empty lines
		if (buffer[0] == '\0') continue;

		if (get_fasta_header(buffer, track_name)) {

			/// succeeded in getting header: new track

			// clean up last stretch
			if (count >= nrepeats) {
				// record current mononucleotide stretch
				out
					<< track_name << delim
					<< start << delim
					<< (start + count) << delim
					<< prev_nuc << delim
					<< count << endl;
			}
			
			// initialize for new track
			start = pos = 0;
			prev_nuc = '\0';
			count = 1;

		} else {

			/// continue with current track

			if (track_name[0] == '\0') {
				throw runtime_error("Error: FASTA must have at least one track name");
			}

			for (size_t i = 0; i < buffer_size; ++i) {

				char nucleotide = buffer[i];

				if (nucleotide == '\0') break;

				if (nucleotide == prev_nuc && nucleotide != ignore_nuc) {

					++count;

				} else {

					// clean up last stretch
					if (count >= nrepeats) {
						// record current mononucleotide stretch
						out
							<< track_name << delim
							<< start << delim
							<< (start + count) << delim
							<< prev_nuc << delim
							<< count << endl;
					}

					// prepare for next mononucleotide stretch
					prev_nuc = nucleotide;
					start = pos;
					count = 1;
				}

				++pos;

			}

		}

	}

	// clean up last stretch
	if (count >= nrepeats) {
		// record current mononucleotide stretch
		out
			<< track_name << delim
			<< start << delim
			<< (start + count) << delim
			<< prev_nuc << delim
			<< count << endl;
	}

}


int main(int argc, char* argv[]) {

	if (argc < 4) {
		cout << "Usage: mnrf <nrepeats> <input.fa> <output.bed>" << endl;
		return 0;
	}

	ifstream infile;
	ofstream outfile;

	int nrepeats = atoi(argv[1]);
	infile.open(argv[2]);
	outfile.open(argv[3]);

	if (!infile.is_open()) {
		throw runtime_error("Error: could not open input file");
	}

	if (!outfile.is_open()) {
		throw runtime_error("Error: could not open output file");
	}

	if (nrepeats < 2) {
		throw domain_error("Error: nrepeats must be at least 2");
	}

	mnrf(nrepeats, infile, outfile);

	infile.close();
	outfile.close();

	return 0;
}

