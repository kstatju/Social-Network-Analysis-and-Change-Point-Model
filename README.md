# Network-Influence-Measures
Python Functions for Network Influence Measures
# Import Netwrok_Centrality functions 
import Network_Centrality as cent

# Call Dynamic Centrality function and it will return a matrix
# Input data (x) should be three dimentional and dimension one represent the time variable

B = cent.grindrod(x = Atensor, alpha = 0.1)

# Call Katz Centrality function and it will return a matrix
# Input data (x) should be square matrix

kz = cent.katz(x = Atensor, .015)

# Call function "dy_rank_degree" to find the rank of nodes for the Dynamic Centrality
cent.dy_rank_degree(B)

# Call function "kz_rank_degree" to find the rank of nodes for the Katz Centrality
cent.kz_rank_degree(kz)
