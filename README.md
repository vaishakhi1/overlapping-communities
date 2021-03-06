# overlapping-communities
spectral clustering when you have multiple community memberships in a network

Community detection for real-world networks (an ongoing project)

Project Description.
The purpose of these codes is to study the behavior of spectral detection for community detection for stochastic block models. The ability to identify the community memberships of the nodes depends on the parameters of the SBM. If the separation between the in-degree and the out-degree is too small, compared to the average degree and size of the graph, no polynomial time algorithm can recover the community memberships of the nodes. In the case of two communities, a detectability threshold marks a second order transition beyond which the performance of the spectral clustering algorithm improves with the separation. Below the threshold, the algorithm does not do better than random chance.

In this work, we explore the relationship between the SBM parameters and the performance of the spectral clustering algorithm in the presence of overlap in communities. In the case of SBMs with two communities, the overlapping nodes are in both communities. The specific spectral clustering algorithm uses the nonbacktracking matrix. 

We explore two possible approaches to classifying overlapping communities. In one approach, we treat the overlapping nodes as a separate community, and cluster the nodes into three communities. In an alternate, hierarchical approach, the overlapping nodes are separated first and the remaining nodes are classified into two communities. 

The inputs to the code specify the size of the graph, the average degree of the nonoverlapping graph. The range of the separation is also specified (in-degree – out-degree).

The code outputs a matrix with the performance of the spectral clustering algorithm as a function of the separation and overlap.  We also compare the observed threshold with the theoretical detectability threshold based on the Kesten Stigum bound.

Dependency Requirements -  Matlab – R2015b and above


Run SpectralClusteringhierarchical.m for the hierarchical clustering or SpectralClusteringKMM.m for the approach of clustering into three groups
