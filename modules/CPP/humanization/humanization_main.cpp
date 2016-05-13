* Written in 2014 by Tong Li of the Costas Maranas Lab in the Chemical
* Engineering Department of the Pennsylvania State University

*This file contains classes for deimmunization during antibody design */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Deimmunization.h"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

using namespace std;

int main(int argc, char * argv[])
{

    string seqFile = argv[1];
    string spotsFile = argv[2];
    string kindsFile = argv[3];
    Deimmunization de = Deimmunization();
    string databaseName = "Human_9mer_Sequences.txt";
    ifstream in;
    in.open(seqFile.c_str());
    string sequence;
    in >> sequence;
    boost::trim(sequence);
    in.close();
    //cout << sequence << endl;
    in.open(spotsFile.c_str());
    vector<int> positions;
    string line;
    in >> line;
    while(!in.eof())
    {
        boost::trim(line);
        positions.push_back(atoi(line.c_str()));
        in >> line;
    }
    in.close();
    //cout << positions.size() << endl;
    //for (vector<int>::iterator it = positions.begin(); it != positions.end(); it++)
    //{
    //    cout << *it << endl;
    //}

    in.open(kindsFile.c_str());
    vector< vector<string> > kinds;
    string kindsLine;
    in >> kindsLine;
    while(!in.eof())
    {
        boost::trim(kindsLine);
        string st = kindsLine.substr(0, kindsLine.length()-1);
        vector<string> strs;
        boost::split(strs, st, boost::is_any_of(":"));
        kinds.push_back(strs);
        in >> kindsLine;
    }
    in.close();
    //cout << "sequence:" << sequence << endl;
    //cout << "positions: " << endl;
    for(vector<int>::iterator it = positions.begin(); it!=positions.end(); it++)
     {
         cout << *it << endl;
     }
    for(vector< vector<string> >::iterator it = kinds.begin(); it!=kinds.end(); it++)
     {
         cout << (*it)[0] << endl;
     }
    de.make_humanization_cuts(sequence, positions, kinds);

    return 0;
}

