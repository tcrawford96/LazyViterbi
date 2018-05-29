#include <math.h>
#include <iostream>
#include <vector>
#include <map>
#include <list>
#include <utility>
#include <forward_list>
#include <tuple>
//#include "LazyViterbi.h"

using namespace std;

//for a given node (state at a given time t), returns the metric to next_state at time t+1
int compute_metric(tuple <int,int,int> shadow, int nextState, vector<vector<int>> TM, vector<vector<int>> EM){
//	int metric = TM[get<0>(shadow)][nextState] + EM[nextState][get<1>(shadow)];
	return 1; //metric;
}

//inserts a shadow node into the PQ based on its metric
void insert_PQ(int metric,tuple <int,int,int> shadow, vector<list<tuple <int,int,int>>> PQ){
/*	for(int i = 0;i<PQ.size();i++){
		if(i == metric){
			PQ[i].push_front(shadow);
			break;
		}
	}
*/}

//inserts a node into the trellis using a key of pair<int state, int time> and value (int prev_st)
void insert_Trellis(tuple <int,int,int> shadow, map<pair<int,int>, int> trellis){
/*	pair<int,int> key = {get<0>(shadow), get<1>(shadow)};
	int value = get<2>(shadow);
	trellis.insert({key,value});
*/}

//given a shadow node from the PQ, expands the node by adding it as a real node in the trellis and creating shadow nodes for each state at time t+1
void expand(tuple <int,int,int> shadow, vector<list<tuple <int,int,int>>> PQ, vector<vector<int>> TM, vector<vector<int>> EM, map<pair<int,int>, int> trellis){
	int t = get<1>(shadow)+1;
	insert_Trellis(shadow,trellis);

	//find minimum metric to scale other metrics by
	int minMetric = compute_metric(shadow,0,TM,EM);
	for(int i=0; i < 127; i++){
		int metric = compute_metric(shadow,i,TM,EM);
		if(metric < minMetric)
			minMetric = metric;
	}

	//calculate relative metrics and create shadow nodes in the PQ
	for(int i=0; i < 127; i++){
		int metric = compute_metric(shadow,i,TM,EM)-minMetric;
		int prev_st = get<0>(shadow);
		tuple <int,int,int> nextNode {i, t, prev_st};
		insert_PQ(metric,nextNode,PQ);		
	}		
}


//returns true if shadow is in trellis, false otherwise
bool trellis_Lookup(tuple <int,int,int> shadow, map<pair<int,int>, int> trellis){
/*	if(trellis.find(pair<int,int> {get<0>(shadow),get<1>(shadow)}) == trellis.end())
		return false;
	else	
		return true;
*/
	return true;
}

//performs a traceback starting at finalNode, stepping backwards for the length of the time Buffer, and only recovering the oldest few states for the solution vector
void traceback(map<pair<int,int>, int> trellis, pair <int,int> finalNode, int timeBuffer, int statesRecovered, int (*soln)){
/*	int finalTime = finalNode.second;
	pair<int,int> currNode = finalNode;
	for(int t = finalTime; t >= (finalTime-timeBuffer); t--){
		currNode = {trellis[currNode],t-1};
		if(t <= statesRecovered){
			soln[t] = currNode.first;
		}
	}
*/}





void LazyViterbi(vector<int> data, vector<int> soln,int time) {

vector<vector<int>> TM;
vector<vector<int>> EM;
map<pair<int,int>, int> trellis;
vector<list<tuple <int,int,int>>> PQ;

int N = 127; //number of states in markov model
int M = 255; //number of values allowed in data (0-255)

//contstruct metric matrices

for(int i = 0; i < N; i++) {
	vector<int> subVec;
	TM.push_back(subVec);
	for(int j = 0; j < N;j++) {
		if(i != j) 
			TM[i].push_back(3);
		else
			TM[i].push_back(0);
	}
}


// j corresponds to the variation from the state i
for(int i = 0; i < N; i++) {
	vector<int> subVec;
	EM.push_back(subVec);
	for(int j = 0; j < M; j++) {
		EM[i].push_back(round(log10((sqrt(2.0*M_PI*100.0))/exp(-pow((j-i),2)/200.0))));
	cout << j << '-' << i << ' ';
	}
}

cout << 'a';

//run decoder for first time buffer

int timeBuffer = 500;
int t = 0;
tuple <int,int,int> firstNode {data[0],0,0};
expand(firstNode, PQ, EM, TM ,trellis);

while(t < timeBuffer){
	tuple<int,int,int> currNode = (PQ[0]).front();
	(PQ[0]).pop_front();
	
	if( trellis_Lookup(currNode, trellis) == false){
		expand(currNode,PQ,EM,TM,trellis);
		t = get<1>(currNode);
	}
} 


}

int main(){

vector<int> data = {0,0,0};
vector<int> soln = {0,0,0};
int time = 0;

LazyViterbi(data,soln,time);
return 1;
}
