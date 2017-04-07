#include "phylotree.h"
#include "gtest/gtest.h"

TEST(tree1, NodesHitI) {
    node chicken("chicken",0.0537,0,0,0,1);
    node saimo("saimo",0.0561,1,0,0,0);
    node xenopus("xenopus",0.0691,0,1,0,0);
    node iguana("iguana",0.0483,0,0,0,1);
    node turtle("turtle",0.0151,0,1,0,0);
    node int_3(0.0144,&saimo,&xenopus);
    node int_2(0.0097,&chicken,&int_3);
    node int_1(0.0151,&int_2,&iguana);
    node root(0,&int_1,&turtle);
    
    double *base_freq = new double[4];
    base_freq[0] = 0.25;
    base_freq[1] = 0.25;
    base_freq[2] = 0.25;
    base_freq[3] = 0.25;
    
    calc_likelihood_i(&root, base_freq);
    EXPECT_EQ(1, chicken.hits);
    EXPECT_EQ(1, saimo.hits);
    EXPECT_EQ(1, xenopus.hits);
    EXPECT_EQ(1, iguana.hits);
    EXPECT_EQ(1, turtle.hits);
    EXPECT_EQ(1, int_3.hits);
    EXPECT_EQ(1, int_2.hits);
    EXPECT_EQ(1, int_1.hits);
    EXPECT_EQ(1, root.hits);
    
    delete[] base_freq;
}

TEST(tree1, NodesHitR) {
    node chicken("chicken",0.0537,0,0,0,1);
    node saimo("saimo",0.0561,1,0,0,0);
    node xenopus("xenopus",0.0691,0,1,0,0);
    node iguana("iguana",0.0483,0,0,0,1);
    node turtle("turtle",0.0151,0,1,0,0);
    node int_3(0.0144,&saimo,&xenopus);
    node int_2(0.0097,&chicken,&int_3);
    node int_1(0.0151,&int_2,&iguana);
    node root(0,&int_1,&turtle);
    
    double *base_freq = new double[4];
    base_freq[0] = 0.25;
    base_freq[1] = 0.25;
    base_freq[2] = 0.25;
    base_freq[3] = 0.25;
    
    calc_likelihood_r(&root, base_freq);
    EXPECT_EQ(1, chicken.hits);
    EXPECT_EQ(1, saimo.hits);
    EXPECT_EQ(1, xenopus.hits);
    EXPECT_EQ(1, iguana.hits);
    EXPECT_EQ(1, turtle.hits);
    EXPECT_EQ(1, int_3.hits);
    EXPECT_EQ(1, int_2.hits);
    EXPECT_EQ(1, int_1.hits);
    EXPECT_EQ(1, root.hits);
    
    delete[] base_freq;
}

TEST(tree1, sameBetweenRecursiveIterative) {
    node chicken("chicken",0.0537,0,0,0,1);
    node saimo("saimo",0.0561,1,0,0,0);
    node xenopus("xenopus",0.0691,0,1,0,0);
    node iguana("iguana",0.0483,0,0,0,1);
    node turtle("turtle",0.0151,0,1,0,0);
    node int_3(0.0144,&saimo,&xenopus);
    node int_2(0.0097,&chicken,&int_3);
    node int_1(0.0151,&int_2,&iguana);
    node root(0,&int_1,&turtle);
    
    double *base_freq = new double[4];
    base_freq[0] = 0.25;
    base_freq[1] = 0.25;
    base_freq[2] = 0.25;
    base_freq[3] = 0.25;
    
    EXPECT_EQ(calc_likelihood_i(&root, base_freq),calc_likelihood_r(&root, base_freq));
    
    delete[] base_freq;
}

TEST(tree1, correctPrint) {
	node chicken("chicken",0.0537,0,0,0,1);
    node saimo("saimo",0.0561,1,0,0,0);
    node xenopus("xenopus",0.0691,0,1,0,0);
    node iguana("iguana",0.0483,0,0,0,1);
    node turtle("turtle",0.0151,0,1,0,0);
    node int_3(0.0144,&saimo,&xenopus);
    node int_2(0.0097,&chicken,&int_3);
    node int_1(0.0151,&int_2,&iguana);
    node root(0,&int_1,&turtle);
	string string1 = "(((chicken:0.0537,(saimo:0.0561,xenopus:0.0691)),iguana:0.0483),turtle:0.0151);";
	string string2 = print_newick(&root);
	
    EXPECT_TRUE(string1==string2);
    
}

TEST(tree2, NodesHitI) {
    node a("A",0.1,1,0,0,0);
    node c("C",0.1,0,1,0,0);
    node t("T",0.1,0,0,1,0);
    node g("G",0.1,0,0,0,1);
    node ca(0.1,&c,&a);
    node tca(0.1,&t,&ca);
    node root(0,&g,&tca);
    
    double *base_freq = new double[4];
    base_freq[0] = 0.25;
    base_freq[1] = 0.25;
    base_freq[2] = 0.25;
    base_freq[3] = 0.25;
    
    calc_likelihood_i(&root, base_freq);
    EXPECT_EQ(1, a.hits);
    EXPECT_EQ(1, c.hits);
    EXPECT_EQ(1, t.hits);
    EXPECT_EQ(1, g.hits);
    EXPECT_EQ(1, ca.hits);
    EXPECT_EQ(1, tca.hits);
    EXPECT_EQ(1, root.hits);
    
    delete[] base_freq;
}
TEST(tree2, NodesHitR) {
    node a("A",0.1,1,0,0,0);
    node c("C",0.1,0,1,0,0);
    node t("T",0.1,0,0,1,0);
    node g("G",0.1,0,0,0,1);
    node ca(0.1,&c,&a);
    node tca(0.1,&t,&ca);
    node root(0,&g,&tca);
    
    double *base_freq = new double[4];
    base_freq[0] = 0.25;
    base_freq[1] = 0.25;
    base_freq[2] = 0.25;
    base_freq[3] = 0.25;
    
    calc_likelihood_i(&root, base_freq);
    EXPECT_EQ(1, a.hits);
    EXPECT_EQ(1, c.hits);
    EXPECT_EQ(1, t.hits);
    EXPECT_EQ(1, g.hits);
    EXPECT_EQ(1, ca.hits);
    EXPECT_EQ(1, tca.hits);
    EXPECT_EQ(1, root.hits);
    
    delete[] base_freq;
}

TEST(tree2, sameBetweenRecursiveIterative) {
    node a("A",0.1,1,0,0,0);
    node c("C",0.1,0,1,0,0);
    node t("T",0.1,0,0,1,0);
    node g("G",0.1,0,0,0,1);
    node ca(0.1,&c,&a);
    node tca(0.1,&t,&ca);
    node root(0,&g,&tca);
    
    double *base_freq = new double[4];
    base_freq[0] = 0.25;
    base_freq[1] = 0.25;
    base_freq[2] = 0.25;
    base_freq[3] = 0.25;
    
    EXPECT_EQ(calc_likelihood_i(&root, base_freq), calc_likelihood_r(&root, base_freq));
    
    delete[] base_freq;
}

TEST(tree2, correctPrint) {
    node a("A",0.1,1,0,0,0);
    node c("C",0.1,0,1,0,0);
    node t("T",0.1,0,0,1,0);
    node g("G",0.1,0,0,0,1);
    node ca(0.1,&c,&a);
    node tca(0.1,&t,&ca);
    node root(0,&g,&tca);
    
	string string1 = "(G:0.1,(T:0.1,(C:0.1,A:0.1)));";
	string string2 = print_newick(&root);
	
    EXPECT_TRUE(string1==string2);
    
}
