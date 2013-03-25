#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <cstring>
#include <map>
#include <cassert>
#include "super.h"

using namespace std;

void read_file(string name, bool ** &conn, int &size) {
	ifstream f(name.c_str());
	char buffer[80];
	assert(f.is_open());
	set<int> v;
	multimap<int, int> e;

	while(!f.eof()) {
		f.getline (buffer, 70);
		if(buffer[0] == 'e') {
			int vi, vj;
			sscanf(buffer, "%*c %d %d", &vi, &vj);
			v.insert(vi);
			v.insert(vj);
			e.insert(make_pair(vi, vj));
		}
	}
	size = v.size();

	conn = new bool *[size];
	for(size_t i = 0; i < size; i++) {
		conn[i] = new bool[size];
		memset(conn[i], 0, size * sizeof(bool));
	}
	for(multimap<int, int>::iterator it = e.begin(); it != e.end(); it++) {
		conn[it->first - 1][it->second -1] = true;
		conn[it->second - 1][it->first - 1] = true;
	}

	cout << "|E| = "<< e.size() << " |V| = " << v.size() << " p = " << (double) e.size() / (v.size() * (v.size() - 1) /2) <<endl;
	f.close();
}

int main(int argc, char *argv[]) {
	assert(argc == 2);
	cout<<"args = "<<argv[1]<<endl;
	bool **conn;
	int size;
	read_file(argv[1], conn, size);
	
	clock_t start1, start2;
	int *qmax;
	int qsize;
	
	start1 = time(NULL);
	start2 = clock();
	Maxclique m(conn, size);
	m.mcqdyn(qmax, qsize);

	sort(qmax, qmax + qsize);
	cout << "Maximum clique: ";
	for (int i = 0; i < qsize; i++) 
		cout << qmax[i] + 1 << " ";
	cout<<endl;
	for(int i = 0; i < qsize; ++i) {
		for(int j = i + 1; j < qsize; ++j)
			if(!conn[qmax[i]][qmax[j]])
				cout<<"Failure!!!!!"<<endl;
	}

	cout << "Size = " << qsize << endl;
	cout << "Number of steps = " << m.steps() << endl;
	cout << "Time = " << difftime(time(NULL), start1) << endl;
	cout << "Time (precise) = " << ((double) (clock() - start2)) / CLOCKS_PER_SEC << endl <<endl;
	delete [] qmax;

	for(int i = 0; i < size; i++)
		delete [] conn[i];
	delete [] conn;

	return 0;
}
