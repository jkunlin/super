#ifndef	SUPER_H 
#define SUPER_H

#include <iostream>
#include <algorithm>
#include <assert.h>
#ifdef DBG
using namespace std;
#endif

class Maxclique {
	const bool* const* e;
	int pk, level;
	const float Tlimit;
	class Vertices {
		class Vertex {
			int i, d;
			public:
			void set_i(const int ii)  { i = ii; }
			int get_i() const { return i; }
			void set_degree(int dd) { d = dd; }
			int get_degree() const { return d; }
		};
		Vertex *v;
		int sz;
		static bool desc_degree(const Vertex vi, const Vertex vj) { return (vi.get_degree() > vj.get_degree()); }
		public:
#ifdef DBG
		void dbg_v(const string msg="") const {
			std::cout << msg << " Vertices: [";
			for (int i=0; i < sz; i++) 
				std::cout << "(" << v[i].get_i() << "," << v[i].get_degree() << ") ";
			std::cout << "]" << std::endl;
		}
#endif
		Vertices(int size) : sz(0) { v = new Vertex[size]; }
		~Vertices () {}
		void dispose() { if (v) delete [] v; }
		void sort() { std::sort(v, v+sz, desc_degree); }
		void init_colors();
		void set_degrees(Maxclique&);
		int size() const { return sz; }
		void resize(const int size) { sz = size; }
		void push(const int ii) { v[sz++].set_i(ii); };
		void pop() { sz--; };
		Vertex& at(const int ii) const { return v[ii]; };
		Vertex& end() const { return v[sz - 1]; };
	};
	class ColorClass {
		int *i;
		int sz;
		public:
#ifdef DBG
		void dbg_i(const string msg="") const {
			std::cout << msg << " Class: [";
			for (int ii=0; ii < sz; ii++) 
				std::cout << i[ii] << " ";
			std::cout << "]" << std::endl;
		}
#endif
		ColorClass() : sz(0), i(0) {}
		ColorClass(const int sz) : sz(sz), i(0) { init(sz); }
		~ColorClass() { if (i) delete [] i;}
		void init(const int sz) { i = new int[sz]; rewind(); }
		void push(const int ii) { i[sz++] = ii; };
		void pop() { sz--; };
		void pop(const int p_ci) { i[p_ci] = -1; }
		void rewind() { sz = 0; };
		int size() const { return sz; }
		int& at(const int ii) const { return i[ii]; }
		int& end() const { return i[sz - 1]; }
		ColorClass& operator=(const ColorClass& dh) {
			for (int j = 0; j < dh.sz; j++) i[j] = dh.i[j];
			sz = dh.sz;
			return *this;
		}
	};
	Vertices V;
	ColorClass *C, QMAX, Q;
#ifdef SAT
	std::set<int> *sat;
#endif
	class StepCount {
		int i1, i2;
		public:
		StepCount() : i1(0), i2(0) {}
		void set_i1(const int ii)  { i1 = ii; }
		int get_i1() const { return i1; }
		void set_i2(const int ii)  { i2 = ii; }
		int get_i2() const { return i2; }
		void inc_i1()  { i1++; }
	};
	StepCount *S;
	bool connection(const int i, const int j) const { return e[i][j]; }
	bool conflict(const int, const ColorClass&);
	void cut(const Vertices&, Vertices&);
	void cut_new(Vertices&, Vertices&, int);
	void color_sort(Vertices&, bool = false);
	void re_color(int, int);
	int count_conflict(int, const ColorClass&, int&);
	bool mcs_conflict(int, const ColorClass&);

	void re_color_sort(Vertices&, Vertices&, bool = true);
#ifdef SAT
	void sat_color_sort(Vertices&);
#endif
	void expand(Vertices);
	void expand_dyn(Vertices);
	void expand_dyn(Vertices, Vertices);
	void _mcq(int*&, int&, bool);
	void degree_sort(Vertices &R) { R.set_degrees(*this); R.sort(); }
	public:
#ifdef DBG
	void dbg_C() const {
		for (int i=0; i < V.size(); i++) {
			std::cout << "C["<< i << "] : ";
			C[i].dbg_i();
		}
	}
	void dbg_conn() const {
		for (int i=0; i < V.size(); i++) {
			for (int j=0; j < V.size(); j++) {
				std::cout <<e[i][j];
			}
			std::cout<< std::endl;
		}
	}
#endif
	Maxclique(const bool* const*, const int, const float=0.025);
	int steps() const { return pk; }
	void mcq(int* &maxclique, int &sz) { _mcq(maxclique, sz, false); }
	void mcqdyn(int* &maxclique, int &sz) { _mcq(maxclique, sz, true); }
	~Maxclique() {
		if (C) delete [] C;
		if (S) delete [] S;
		V.dispose();
	};
};

Maxclique::Maxclique (const bool* const* conn, const int sz, const float tt) : pk(0), level(1), Tlimit(tt), V(sz), Q(sz), QMAX(sz) {
	assert(conn!=0 && sz>0);
	for (int i=0; i < sz; i++) V.push(i);
	e = conn;
	C = new ColorClass[sz + 1];
	for (int i=0; i < sz + 1; i++) C[i].init(sz);
	S = new StepCount[sz];
#ifdef SAT	
	sat = new std::set<int> [sz];
#endif
}

void Maxclique::_mcq(int* &maxclique, int &sz, bool dyn) { 
	V.set_degrees(*this);
	V.sort();
	V.init_colors();
	if (dyn) {
		/**********construct funtion have done the init
		  for (int i=0; i < V.size() + 1; i++) {
		  S[i].set_i1(0);
		  S[i].set_i2(0);
		  }
		 */
		expand_dyn(V, V);
	}
	else
		expand(V);
	maxclique = new int[QMAX.size()]; 
	for (int i=0; i<QMAX.size(); i++) { 
		maxclique[i] = QMAX.at(i);
	}
	sz = QMAX.size();
}

void Maxclique::Vertices::init_colors() { 
	const int max_degree = v[0].get_degree();
	for (int i = 0; i < max_degree; i++)
		v[i].set_degree(i);
	for (int i = max_degree; i < sz; i++)
		v[i].set_degree(max_degree);
}

void Maxclique::Vertices::set_degrees(Maxclique &m) { 
	for (int i=0; i < sz; i++) {
		int d = 0;
		for (int j=0; j < sz; j++)
			if (m.connection(v[i].get_i(), v[j].get_i())) d++;
		v[i].set_degree(d);
	}
}

bool Maxclique::conflict(const int pi, const ColorClass &A) {
	for (int i = 0; i < A.size(); i++)
		if (connection(pi, A.at(i)))
			return true;
	return false;
}

void Maxclique::cut(const Vertices &A, Vertices &B) {
	int pi = A.end().get_i();
	for (int i = 0; i < A.size() - 1; i++) {
		if (connection(pi, A.at(i).get_i())){
			B.push(A.at(i).get_i());
			//			B.end().set_degree(A.at(i).get_degree());
		}
	}
}
/*
void Maxclique::cut_new(Vertices &Va, Vertices &Vp, int pi) {
	bool find = false;
	for(int i = 0; i < Va.size(); ++i) {
		if(Va.at(i).get_i() == pi) {
			find = true;
			continue;
		}
		if(connection(pi, Va.at(i).get_i()))
			Vp.push(Va.at(i).get_i());
		if(find) {
			Va.at(i - 1).set_i(Va.at(i).get_i());
		}
	}
	Va.pop();
}
*/
void Maxclique::cut_new(Vertices &Va, Vertices &Vp, int pi) {
	for(int i = 0; i < Va.size(); ++i) {
		if(Va.at(i).get_i() == pi)
			Va.at(i).set_i(-1);
		if(Va.at(i).get_i() != -1 && connection(pi, Va.at(i).get_i()))
			Vp.push(Va.at(i).get_i());
	}
}

#ifdef SAT
void Maxclique::sat_color_sort(Vertices &R) {
	int maxno = 1;
	int min_k = QMAX.size() - Q.size() + 1;
	C[1].rewind();
	C[2].rewind();
	int k = 1;
	int max_sat;
	int max_sat_i;

	//	R.set_degrees(*this);
	/*	
		for(int i = R.size() - 1; i > 0; --i) {
		if(R.at(i).get_degree() > R.at(i - 1).get_degree()) {
		std::swap(R.at(i), R.at(i - 1));
		}

		sat[i].clear();
		}
	 */	
	/*
	   max_sat = R.at(0).get_degree();
	   max_sat_i = 0;
	   for(int i = 0; i < R.size(); ++i) {
	   if(R.at(i).get_degree() > max_sat) {
	   max_sat = R.at(i).get_degree();
	   max_sat_i = i;
	   }
	   sat[i].clear();
	   }
	   std::swap(R.at(0), R.at(max_sat_i));
	 */
	C[1].push(R.at(0).get_i());
	R.at(0).set_degree(k < min_k ? 0 : 1);


	int j = 1;
	int pi;
	for (int i=1; i < R.size(); i++) {

		pi = R.at(i - 1).get_i();
		//update sat, set R.at(i) the largest sat vertex
		max_sat = 0;
		max_sat_i = i;
		for(int si = i; si < R.size(); ++si) {
			if(connection(pi, R.at(si).get_i())) {
				sat[si].insert(1);
			}
			if(sat[si].size() > max_sat) {
				max_sat = sat[si].size();
				max_sat_i = si;
			}
		}
		// connection component checked
		//if(max_sat = 0) {
		//do something
		//	}

		if(max_sat_i != i) {
			std::swap(R.at(i), R.at(max_sat_i));
			sat[i].swap(sat[max_sat_i]);
		}

		int pi = R.at(i).get_i();
		k = 1;
		while (conflict(pi, C[k]))
			k++;
		if (k > maxno) {
			maxno = k;
			C[maxno + 1].rewind();
		}
		C[k].push(pi);
		if (k < min_k) {
			R.at(j++).set_i(pi);
		}
	}

	if (j > 1) R.at(j-1).set_degree(0);
	if (min_k <= 0) min_k = 1;
	for (k = min_k; k <= maxno; k++)
		for (int i = 0; i < C[k].size(); i++) {
			R.at(j).set_i(C[k].at(i));
			R.at(j++).set_degree(k);
		}

}
#endif

int Maxclique::count_conflict(int pi, const ColorClass& A, int &q_ci) {
	int num = 0;
	bool first = true;
	for(int i = 0; i < A.size(); ++i) {
		if(A.at(i) != -1 && connection(pi, A.at(i))) {
			if(first) {
				first = false;
				q_ci = i;
			}
		}
			++num;
	}
	return num;
}
bool Maxclique::mcs_conflict(int pi, const ColorClass& A) {
	for(int i = 0; i < A.size(); ++i) {
		if(A.at(i) != -1 && connection(pi, A.at(i))) {
			return true;
		}
	}
	return false;
}
void Maxclique::re_color(int k, int min_k) {
	int pi = C[k].end();
	int q_ci = -1;
	int qi = -1;
	for(int k1 = 1; k1 < min_k - 1; ++k1) {
		if(count_conflict(pi, C[k1], q_ci) == 1) {
			qi = C[k1].at(q_ci);
			for(int k2 = k1 + 1; k2 <= min_k - 1; ++k2) {
				if(!mcs_conflict(qi, C[k2])) {
					C[k].pop();
					C[k1].push(pi);
					C[k2].push(qi);
					C[k1].pop(q_ci);
				}
			}
		}
	}
}

void Maxclique::re_color_sort(Vertices &Va, Vertices &R, bool recolor) {
	int maxno = 1;
	int min_k = QMAX.size() - Q.size() + 1;
	C[1].rewind();
	C[2].rewind();
	int k = 1;
	
	for (int i=0; i < Va.size(); i++) {
		int pi = Va.at(i).get_i();
		k = 1;
		while (mcs_conflict(pi, C[k]))
			k++;
		if (k > maxno) {
			maxno = k;
			C[maxno + 1].rewind();
		}
		C[k].push(pi);
		// assert maxno >= min_k
		
		if (recolor && k >= min_k && k == maxno) {
			re_color(k, min_k);

			if(C[maxno].size() == 0)
				maxno -= 1;
		}
		
	}
	int j = Va.size() - 1;
	if(min_k <= 0) min_k = 1;
	//last color maybe better than this one, try it when no bug
	for(k = maxno; k >= min_k; --k) {
		for(int i = 0; i < C[k].size(); ++i) {
			if(C[k].at(i) == -1)
				continue;
			R.at(j).set_i(C[k].at(i));
			R.at(j--).set_degree(k);
		} 
	}
	if(j >= 0) {
		R.at(j).set_i(C[k - 1].at(C[k - 1].size() - 1)); //it maybe -1, but no problem i think
		R.at(j).set_degree(min_k - 1);
	}

}
//origin color,modify branch strage, it is right
void Maxclique::color_sort(Vertices &R, bool sorted) {
	int j = 0;
	int maxno = 1;
	int min_k = QMAX.size() - Q.size() + 1;
	C[1].rewind();
	C[2].rewind();
	int k = 1;
	for (int i=0; i < R.size(); i++) {
		int pi = R.at(i).get_i();
		k = 1;
		while (conflict(pi, C[k]))
			k++;
		if (k > maxno) {
			maxno = k;
			C[maxno + 1].rewind();
		}
		C[k].push(pi);
		if (k < min_k) {
			R.at(j++).set_i(pi);
		}
	}
	if (j > 0) R.at(j-1).set_degree(0);
	if (min_k <= 0) min_k = 1;
	for (k = min_k; k <= maxno; k++)
		for(int i = C[k].size() - 1; i >= 0; --i) {
			R.at(j).set_i(C[k].at(i));
			R.at(j++).set_degree(k);
		}

	//	   for (int i = 0; i < C[k].size(); i++) {
	//	   R.at(j).set_i(C[k].at(i));
	//	   R.at(j++).set_degree(k);
	//	   }
}

void Maxclique::expand(Vertices R) {
	while (R.size()) {
		if (Q.size() + R.end().get_degree() > QMAX.size()) {
			Q.push(R.end().get_i());
			Vertices Rp(R.size());
			cut(R, Rp);
			if (Rp.size()) {
				color_sort(Rp);
				pk++;
				expand(Rp);
			}
			else if (Q.size() > QMAX.size()) { 
				std::cout << "step = " << pk << " current max. clique size = " << Q.size() << std::endl; 
				QMAX = Q;
			}    
			Rp.dispose();
			Q.pop();
		}
		else {
			return;
		}
		R.pop();
	}
}

void Maxclique::expand_dyn(Vertices R) {
	S[level].set_i1(S[level].get_i1() + S[level - 1].get_i1() - S[level].get_i2());
	S[level].set_i2(S[level - 1].get_i1());
	while (R.size()) {
		if (Q.size() + R.end().get_degree() > QMAX.size()) {
			Q.push(R.end().get_i());
			Vertices Rp(R.size());
			cut(R, Rp);
			if (Rp.size()) {
				if ((float)S[level].get_i1()/++pk < Tlimit) {
					degree_sort(Rp);
					//		sat_color_sort(Rp);
					//	if(level == 1) std::cout<<"color number ("<<pk<<") = "<<Rp.end().get_degree()<<std::endl;
				}
				/*		else if ((float)S[level].get_i1()/++pk < Tlimit) {
						degree_sort(Rp);
						color_sort(Rp);
						}
				 */
				color_sort(Rp);
				S[level].inc_i1();
				level++;
				expand_dyn(Rp);
				level--;
			}
			else if (Q.size() > QMAX.size()) { 
				std::cout << "step = " << pk << " current max. clique size = " << Q.size() << std::endl; 
				QMAX = Q;
			}    
			Rp.dispose();
			Q.pop();
		}
		else {
			return;
		}
		R.pop();
	}
}

void Maxclique::expand_dyn(Vertices Va, Vertices R) {
	S[level].set_i1(S[level].get_i1() + S[level - 1].get_i1() - S[level].get_i2());
	S[level].set_i2(S[level - 1].get_i1());
	while (R.size()) {
		if (Q.size() + R.end().get_degree() > QMAX.size()) {
			Q.push(R.end().get_i());
			Vertices Vp(R.size());
			cut_new(Va, Vp, R.end().get_i());
			if (Vp.size()) {
				Vertices Rp(Vp.size());
				Rp.resize(Vp.size());
				if ((float)S[level].get_i1()/++pk < Tlimit) {
					degree_sort(Vp);
					//		sat_color_sort(Rp);
			//		if(level == 1) std::cout<<"color number ("<<pk<<") = "<<Rp.end().get_degree()<<std::endl;
				}
				if ((float)S[level].get_i1()/pk < Tlimit / 15) 
					re_color_sort(Vp, Rp);
				else 
					re_color_sort(Vp, Rp, false);

				S[level].inc_i1();
				level++;
				expand_dyn(Vp, Rp);
				level--;
				Rp.dispose();
			}
			else if (Q.size() > QMAX.size()) { 
				std::cout << "step = " << pk << " current max. clique size = " << Q.size() << std::endl; 
				QMAX = Q;
			}    
			Vp.dispose();
			Q.pop();
		}
		else {
			return;
		}
		R.pop();
	//	Va.pop();
	}
}
#endif
