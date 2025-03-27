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

struct MAFAlignment {
    string score; // line containing score information for alignment block
    vector<string> sName; // name of subject
    vector<int> start; // start position in the source sequence
    vector<int> size; // ungapped length of the sequence
    vector<char> strand; // strand of the sequence (+ or -)
    vector<int> srcSize; // total length of the source sequence
    vector<string> sequence; // sequence for alignment line
    
    vector<int> stop;

    bool chainIsGenerated = false;
    vector<string> chainStrings;

    int alignmentIndex;

    // function to add a new line to the alignment block
    void addLine(const std::string& sName, int start, int size, char strand, int srcSize, const std::string& sequence) {
        this->sName.emplace_back(sName);
        this->start.emplace_back(start);
        this->size.emplace_back(size);
        this->strand.emplace_back(strand);
        this->srcSize.emplace_back(srcSize);
        this->sequence.emplace_back(sequence);
        this->stop.emplace_back(start + size);
    }

    bool isEmpty(){
        return this->sequence[0].empty();
    }

    void fillFromSubset(MAFAlignment& alignment, int sStart, int sStop, int rStart, 
                        int qStart, int rSize, int qSize){
        this->sName.clear();
        this->start.clear();
        this->size.clear();

        this->score = alignment.score;
        this->sName = alignment.sName;

        this->start.push_back(rStart);
        this->start.push_back(qStart);

        this->size.push_back(rSize);
        this->size.push_back(qSize);

        this->strand = alignment.strand;
        this->srcSize = alignment.srcSize;
        
        this->sequence.push_back(alignment.sequence[0].substr(sStart, sStop - sStart + 1));
        this->sequence.push_back(alignment.sequence[1].substr(sStart, sStop - sStart + 1));
    }

    void generateChain(){
        if(this->chainIsGenerated)
            return;
        
        for(int i = 1; i < this->sequence.size(); i++){
            string chainString;
            string headerLine = "chain\t";
            headerLine.append(this->score.substr(8,this->score.length() - 8) + "\t");
            headerLine.append(this->sName[0] + "\t" + to_string(this->srcSize[0]) + "\t" + this->strand[0] + "\t");
            headerLine.append(to_string(this->start[0]) + "\t" + to_string(this->stop[0]) + "\t");
            
            if(this->strand[i] == '+'){
                headerLine.append(this->sName[i] + "\t" + to_string(this->srcSize[i]) + "\t" + this->strand[i] + "\t");
                headerLine.append(to_string(this->start[i]) + "\t" + to_string(this->stop[i]) + "\t");
            } else {
                headerLine.append(this->sName[i] + "\t" + to_string(this->srcSize[i]) + "\t" + this->strand[i] + "\t");
                headerLine.append(to_string(this->srcSize[i] - this->stop[i]) + "\t" + to_string(this->srcSize[i] - this->start[i]) + "\t");
            }
            chainString.append(headerLine + to_string(this->alignmentIndex) + "\n");

            string &rSeq = this->sequence[0];
            string &qSeq = this->sequence[i];

            int aligned = 0;
            int qGap = 0;
            int rGap = 0;

            bool inQGap = false;
            bool inRGap = false;

            for(int j = 0; j < rSeq.length(); j++){
                if(rSeq[j] != '-' && qSeq[j] != '-'){
                    if(inQGap || inRGap){
                        chainString.append(to_string(aligned) + "\t" + to_string(qGap) 
                        + "\t" + to_string(rGap) + "\n");
                        aligned = qGap = rGap = 0;
                        inQGap = inRGap = false;
                    }
                    aligned++;
                } else if(rSeq[j] == '-'){
                    rGap++;
                    inRGap = true;
                } else {
                    qGap++;
                    inQGap = true;
                }
            }
            chainString.append(to_string(aligned) + "\n");
            this->chainStrings.emplace_back(chainString);
        }
        this->chainIsGenerated = true;
    }

    void printChain(){
        this->generateChain();
        for(auto it = this->chainStrings.begin(); it != this->chainStrings.end(); ++it){
            cout << *(it) << "\n";
        }
        cout << "\n";
    }

    MAFAlignment(const std::string& sName, int start, int size, int stop, int index) {
        this->sName.emplace_back(sName);
        this->start.emplace_back(start);
        this->size.emplace_back(size);
        this->stop.emplace_back(stop);
        this->alignmentIndex = index;
    }
    
    MAFAlignment(int index){
        this->alignmentIndex = index;
    }
};

bool compareAlignment(const MAFAlignment&, const MAFAlignment&);

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " [-t <num_thchainStringsreads>] <MAF filename>" << "\n";
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

    ifstream infile(filename);

    if (!infile) {
        cerr << "Error opening alignment " << filename << '\n';
        return 1;
    }
    
    list<MAFAlignment> alignments;

    string line;
    int index = 0;
    MAFAlignment current_alignment = MAFAlignment(index);

    while (getline(infile, line)) {
        // skip comment and track lines
        if (line.empty() || line[0] == 'a') {
            if (!current_alignment.score.empty()) {
                // add current alignment to list and start a new one
                alignments.emplace_back(current_alignment);
                index++;
                current_alignment = MAFAlignment(index);
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
        alignments.emplace_back(current_alignment);
    }

    alignments.sort(compareAlignment);

    infile.close();
    
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (int i = 0; i < alignments.size(); i++) {
        auto it = alignments.begin();
        std::advance(it, i);
        
        // Find which alignment from the MAF contains the bed coords.
        it->alignmentIndex = i;
        it->generateChain();
    }

    for(auto it = alignments.begin(); it != alignments.end(); ++it){
        it->printChain();
    }
    
    return 0;
}

bool compareAlignment(const MAFAlignment& a, const MAFAlignment& b) {
    // Compare sName[1] first
    if (a.sName[0] != b.sName[0])
        return a.sName[0] < b.sName[0];
    
    // If sName[1] is equal, compare start[1]
    return a.start[0] < b.start[0];
}