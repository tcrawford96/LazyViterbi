#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <list>
#include <tuple>
#include <stack>

using namespace std;


class LazyViterbi {
  int N, E, M; // N is number of states in Markov Model; E is number of Emmision states; M is number of quantile for branch metric
  int minEMScore, maxEMScore, minTMScore, maxTMScore; //These variable are used for scaling scores
  int deltaT; // If node is from before this parameter, ignore it
  int stdev; // Standard deviation of emmision distribution
  int tp; // Transition probability
  int current_index, max_t; //Stores index for circular list
  
  vector<char> data; // This vector holds data from 0-255 as character data for efficient storage
  vector<vector<int>> TM; // This matrix holds quantized values of the transition probability (-1*log(P))
  vector<vector<int>> EM; // This matrix holds quantized values of the emmision probability (-1*log(P))
  map<pair<int,int>, int> trellis; // Holds trellis info: time, state, previous state
  vector<stack<tuple <int,int,int>>> PQ; // Holds priority queueu info: time, state, previous state

  void loadData(string);
  void constructEM();
  void constructTM();
  tuple<int,int,int> retNextNode();
  void expandNode(tuple <int,int,int>);
	void traceback(string);
public:
  LazyViterbi(string, int, int, int, int);
};

LazyViterbi::LazyViterbi(string filename, int t_N, int t_E, int t_M, int t_deltaT) {
	N = t_N;
  E = t_E;
  M = t_M;
  deltaT = t_deltaT;
  stdev = 10;
  tp = 4;
  current_index = 0;
  max_t = 0;
  
  loadData(filename);
  constructEM();
  constructTM();
  // Initialize PQ
  for(int i=0; i<M; i++) {
    stack<tuple <int,int,int>> subList;
    PQ.push_back(subList);
  }

  // Add first nodes to PQ
  int datum = int(data[0]);
  for(int i=0; i<N; i++){
    PQ[round(float(EM[i][datum] - minEMScore)/(maxEMScore - minEMScore)*M)].push(make_tuple(1,i,0));
  }
	int t=0;
  while(max_t < data.size()-1) {
    auto node = retNextNode();
    expandNode(node);
  }

	traceback(filename);
  cout << "Done" << endl;
}

void LazyViterbi::expandNode(tuple <int,int,int> node) {
  int t = get<0>(node);
  if (t+1 > max_t)
    max_t = t+1;
  int s = get<1>(node);
  int s_0 = get<2>(node);
  trellis[make_pair(t,s)] = s_0;
  int datum = int(data[t+1]);
//cout << "Expand t=" << t << " s=" << s << " data=" << datum << " c_index=" << current_index << endl;

  for(int i=0;i<N;i++) {
    int score = round(float(EM[i][datum] - minEMScore + TM[i][s])/(maxEMScore - minEMScore + 10)*M) + current_index;
    if (score > M)
      score -= M;
    PQ[score].push(make_tuple(t+1, i, s));
  }
}

tuple<int,int,int> LazyViterbi::retNextNode() {

  bool bad_flag = true;
  tuple<int,int,int> shadow;
  
  while(bad_flag) {
    while(PQ[current_index].empty()){
      current_index+=1;
      if(current_index == M)
				current_index = 0;
    }
    shadow = PQ[current_index].top();
    PQ[current_index].pop();
    if(get<0>(shadow) > max_t - deltaT) {
      pair<int,int> key = {get<0>(shadow), get<1>(shadow)};
      if (trellis.find(key) == trellis.end())
			bad_flag = false;
    }
  }
  return shadow;
}
  
void LazyViterbi::loadData(string datafile) {
  ifstream infile(datafile);
  int temp;

  while (infile >> temp) {
    data.push_back(char(temp));
  }
}

void LazyViterbi::constructEM() {
  minEMScore = 100, maxEMScore = 0;
  for(int i = 0; i < N; i++) {
    vector<int> subVec;
    EM.push_back(subVec);
    for(int j = 0; j < E; j++) {
      int w = round(log10(sqrt(2.0*M_PI*stdev*stdev)) + pow((j-2*i),2)/2/stdev/stdev*log10(M_E));
      EM[i].push_back(w);
      if(EM[i][j] > maxEMScore)
	maxEMScore = EM[i][j];
      if(EM[i][j] < minEMScore)
	minEMScore = EM[i][j];
    }
  }
  cout << minEMScore << " " << maxEMScore << endl;
}

void LazyViterbi::constructTM() {
  for(int i = 0; i < N; i++) {
    vector<int> subVec;
    TM.push_back(subVec);
    for(int j = 0; j < N;j++) {
      if(i == j) 
	TM[i].push_back(0);
      else if ((i-j) == 1 || (j-i) == 1)
	TM[i].push_back(10);
      else
	TM[i].push_back(5);
    } 
  }
  minTMScore = 0;
  maxTMScore = 10;
}

void LazyViterbi::traceback(string datafile){
	string outName = datafile.substr(0, datafile.find_last_of(".")) + "Out.txt";	
	ofstream outfile(outName);	

	int t = max_t;
	int s = (--trellis.end())->second;
	double err = 0.0;

//	outfile.open(outName);
	while(t >= 0){
			outfile << trellis[make_pair(t,s)] << '\n';
cout <<	trellis[make_pair(t,s)] << ' ';
			err += abs(trellis[make_pair(t,s)]-data[t]);
			s = trellis[make_pair(--t,s)];
	}
//	outfile.close();
	cout << err/data.size() << endl;
}

int main() {
  LazyViterbi("ExampleHMM.txt", 127, 255, 9, 30);
}
