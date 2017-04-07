#include "phylotree.h"
#include <sstream>
using namespace std;

// does a post-order traversal of a given phylogenetic tree and prints the newick string
string print_newick(node* root)
{
	string ret;
    node* curr = root;
    stack<node*> nodes;
    node* lastNodeVisited = NULL;

    while(!nodes.empty() || curr!=NULL)
    {
		if(curr!=NULL)
		{
			nodes.push(curr);
			curr = curr->left;
			if(curr!=NULL)
				ret+="(";
		}
		else
		{
			node* peek = nodes.top();
			if(peek->left==NULL && peek->right==NULL)
			{
				ret+=peek->name+":";
				std::ostringstream strs;
				strs << peek->blen;
				ret += strs.str();
			}
			else if(lastNodeVisited==peek->left)
			ret+=",";

			if(peek->right!=NULL && lastNodeVisited!=peek->right)
			curr = peek->right;
			else
			{
			if(lastNodeVisited==peek->right&&lastNodeVisited!=NULL)
				ret+=")";
			lastNodeVisited = nodes.top();
			nodes.pop();
			}	
		}
    }
    ret+=";";
    
    return ret;
}
/*
node* read_newick(string newick)
{
	int ctr=0;
	stack<node*> nodes;
	while(newick[ctr]!=';')
	{
		if(newick[ctr]=='(')
		{
			
		}
		ctr++;
	}
}
*/

// Given a branch length and two nucleobases, calculates the jukes cantor probability
double pij_JC(double branchLen, int ltr1, int ltr2)
{
    if(ltr1 == ltr2)
		return 1/4.0 + (3/4.0)*exp(-((4/3.0)*branchLen));
    else
		return 1/4.0 + (1/4.0)*exp(-((4/3.0)*branchLen));
}	

// an iterative approach to calculating the likelihood at a given leaf with a given transition matrix
double calc_likelihood_i(node* leaf, double* base_freq)
{
    node* curr = leaf;
    stack<node*> nodes;
    node* lastNodeVisited = NULL;

    while(!nodes.empty() || curr!=NULL)
    {
		if(curr!=NULL)
		{
	    	nodes.push(curr);
	    	curr = curr->left;
		}
		else
		{
		    node* peek = nodes.top();
		    if(peek->right!=NULL && lastNodeVisited!=peek->right)
				curr = peek->right;
		    else
		    {
			//ins visit
				if(lastNodeVisited==peek->right && peek->left!=NULL)
				{
				    if(peek->L_vals == NULL)
						peek->L_vals = new double[4];
        	   	    for(int k=0;k<4;k++)
				    {
			     		double sum_lft = 0, sum_rgt = 0;
			    		for(int lft = 0; lft<4; lft++){
			           		sum_lft += pij_JC(peek->left->blen,k,lft) * peek->left->L_vals[lft];
			           		}
			      		for(int rgt = 0; rgt<4; rgt++){
      			            sum_rgt += pij_JC(peek->right->blen,k,rgt) * peek->right->L_vals[rgt];
							}
				  	  	//a=1, c=2, t=3, g=4
				    	peek->L_vals[k] = sum_lft*sum_rgt;
	 			    }
	    		}
	 			peek->hits += 1;
				lastNodeVisited = nodes.top();
				nodes.pop();
	    	}	
		}
    }

    double total_l = 0.0;

    for(int k=0;k<4;k++)
    {
        double freq = base_freq[k];
        double l = leaf->L_vals[k]; 
		total_l += freq*l;

    }
    
    return log(total_l);
    
}

// a recursive method to calculate the likelihood of a given phylogenetic tree
double calc_likelihood_r(node *leaf, double *base_freq)
{
    double total_l = 0.0;
    double *L_node_vals = calc_likelihood_r(leaf);
    for(int k=0;k<4;k++)
    {
        double freq = base_freq[k];
        double l = L_node_vals[k]; 
		total_l += freq*l;
    }
    return log(total_l);
}

//recursive helper function
double* calc_likelihood_r(node *leaf)
{
    //check for end of tree
    if(leaf->right==NULL && leaf->left==NULL)
    {
        leaf->hits += 1;
		return leaf->L_vals;
	}

    else
    {
        node *left_l = leaf->left;
        node *right_l = leaf->right;
		double *left_likelihood = calc_likelihood_r(left_l);
		double *right_likelihood = calc_likelihood_r(right_l);	

        if(leaf->L_vals==NULL)
            leaf->L_vals = new double[4];

        for(int k=0;k<4;k++)
        {

            double sum_lft = 0, sum_rgt = 0;

            for(int lft = 0; lft<4; lft++)
                sum_lft += pij_JC(left_l->blen,k,lft) * left_likelihood[lft];
            for(int rgt = 0; rgt<4; rgt++)
                sum_rgt += pij_JC(right_l->blen,k,rgt) * right_likelihood[rgt];

            //a=1, c=2, t=3, g=4
            leaf->L_vals[k] = sum_lft*sum_rgt;
        }
        leaf->hits += 1;
        return leaf->L_vals;
    }
}

