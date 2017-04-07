#include <iostream>
#include <cmath>
#include <stack>
using namespace std;


struct node
{
	string name;
	//a=0, c=1, t=2, g=3
	double blen;
	node *left;
	node *right;
	double *L_vals;


	node(string n,double blen_in, double l1, double l2, double l3, double l4) : name(n), blen(blen_in),
		left(NULL), right(NULL)
	{
    	L_vals = new double[4];
    	L_vals[0]=l1;
	    L_vals[1]=l2;
    	L_vals[2]=l3;
  		L_vals[3]=l4;
	}
	node(double blen_in, node *left_in, node *right_in) : blen(blen_in), left(left_in),
 		right(right_in), L_vals(NULL)
	{/*empty*/}
  
	~node()
	{
		delete[] L_vals;
	}
	
	int hits = 0;
	
};

string print_newick(node* root);

double pij_JC(double branchLen, int ltr1, int ltr2);

double calc_likelihood_i(node* leaf, double* base_freq);

double* calc_likelihood_r(node *leaf);

double calc_likelihood_r(node *leaf, double *base_freq);

