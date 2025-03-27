// Written by Mohamed El-Walid
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <memory>
#include <list>
#include <vector>
#include <algorithm>

using namespace std;

static char c[128];

struct MAFAlignment {
    string score; // line containing score information for alignment block
    vector<string> sName; // name of subject
    vector<int> start; // start position in the source sequence
    vector<int> size; // ungapped length of the sequence
    vector<char> strand; // strand of the sequence (+ or -)
    vector<int> srcSize; // total length of the source sequence
    vector<string> sequence; // sequence for alignment line

    // function to add a new line to the alignment block
    void addLine(const std::string& sName, int start, int size, char strand, int srcSize, const std::string& sequence) {
        this->sName.emplace_back(sName);
        this->start.emplace_back(start);
        this->size.emplace_back(size);
        this->strand.emplace_back(strand);
        this->srcSize.emplace_back(srcSize);
        this->sequence.emplace_back(sequence);
    }
};

void reverseComplement(MAFAlignment& alignment);
char complement(char);
bool compareAlignment(const MAFAlignment&, const MAFAlignment&);

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " [-t <num_threads>] <filename>" << "\n";
        return 1;
    }

    bool isThreaded = false;
    int num_threads = 1;  // default value for number of threads
    string filename;

    for (int i = 1; i < argc; i++) {
        if(string(argv[i]) == "-t") {
            isThreaded = true;
            num_threads = stoi(argv[i+1]);
            i++;
        } else {
            filename = argv[i];
        }
    }

    c['A'] = 'T'; c['a'] = 'T';
    c['C'] = 'G'; c['c'] = 'G';
    c['G'] = 'C'; c['g'] = 'C';
    c['T'] = 'A'; c['t'] = 'A';
    c['U'] = 'A'; c['u'] = 'A';
    c['M'] = 'K'; c['m'] = 'K';
    c['R'] = 'Y'; c['r'] = 'Y';
    c['W'] = 'W'; c['w'] = 'W';
    c['S'] = 'S'; c['s'] = 'S';
    c['Y'] = 'R'; c['y'] = 'R';
    c['K'] = 'M'; c['k'] = 'M';
    c['V'] = 'B'; c['v'] = 'B';
    c['H'] = 'D'; c['h'] = 'D';
    c['D'] = 'H'; c['d'] = 'H';
    c['B'] = 'V'; c['b'] = 'V';
    c['N'] = 'N'; c['n'] = 'N';
    c['-'] = '-';

    ifstream infile(filename);

    if (!infile) {
        cerr << "Error opening file " << filename << endl;
        return 1;
    }

    list<MAFAlignment> alignments;
    MAFAlignment current_alignment;

    string line;
    while (getline(infile, line)) {
        // skip comment and track lines
        if (line.empty() || line[0] == '#') {
            continue;
        }

        // start of new alignment block
        if (line.substr(0, 1) == "a") {
            if (!current_alignment.score.empty()) {
                // add current alignment to list and start a new one
                alignments.push_back(current_alignment);
                current_alignment = MAFAlignment();
            }
            current_alignment.score = line;
        }

        // s line within an alignment block
        else if (line.substr(0, 1) == "s") {
            string sName;
            int start, size, srcSize;
            char strand, tmp;
            string sequence;

            //create stringstream of line for input into MAFAlignment object
            istringstream lineStream(line);
            lineStream >> tmp >> sName >> start >> size >> strand >> srcSize >> sequence;
            current_alignment.addLine(sName, start, size, strand, srcSize, sequence);
        }
    }

    // add last alignment block to list
    if (!current_alignment.score.empty()) {
        alignments.push_back(current_alignment);
    }

    // Count SNPS and aligned BP
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (int i = 0; i < alignments.size(); i++) {
        auto it = alignments.begin();
        std::advance(it, i);
        if (it->strand[1] == '-') {
            reverseComplement(*it);
            it->strand[1] = '+';
            it->strand[0] = '-';
        }
    }

    alignments.sort(compareAlignment);

    for(auto it = alignments.begin(); it != alignments.end(); ++it){
        cout << it->score << '\n';
        cout << "s\t" << it->sName[1] << '\t' << it->start[1] << '\t' << it->size[1] << '\t' << it->strand[1] << '\t' << it->srcSize[1] << '\t' << it->sequence[1] << '\n';
        cout << "s\t" << it->sName[0] << '\t' << it->start[0] << '\t' << it->size[0] << '\t' << it->strand[0] << '\t' << it->srcSize[0] << '\t' << it->sequence[0] << "\n\n";
    }
    
    return 0;
}

void reverseComplement(MAFAlignment& alignment) {
    for (auto& str : alignment.sequence) {
        int n = str.length();
        // Reverse the string
        reverse(str.begin(), str.end());

        // Apply complement to each character
        for (char& ch : str) {
            ch = complement(ch);
        }
    }
}

char complement(char base) {
    return c[base];
}

bool compareAlignment(const MAFAlignment& a, const MAFAlignment& b) {
    // Compare sName[1] first
    if (a.sName[1] != b.sName[1])
        return a.sName[1] < b.sName[1];
    
    // If sName[1] is equal, compare start[1]
    return a.start[1] < b.start[1];
}