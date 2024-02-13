# EgoDist: a measure for comparing networks

EgoDist is a Matlab function that computes a series of measures that quantify the dissimilarity between two undirected, unweighted networks.
It is defined by processing the distributions, on the network, of three egonet features, namely the degree, the clustering coefficient and the egonet persistence.
The method exploits the statistics of the three features to define one- or multi-dimensional distribution functions, which are then compared to define a distance between networks.
It does not require the alignment of the two networks being compared, which can also have different size.

See the paper for more details:
C. Piccardi, Metrics for network comparison using egonet feature distribution, Scientific Reports, 13, 14657, 2023, https://doi.org/10.1038/s41598-023-40938-4

![image](https://github.com/CarloPiccardi/EgoDist-a-measure-for-comparing-networks/assets/159918290/b260515d-cd26-4424-a5b3-8a4e2d48f91c)

# Usage


function distance = EgoDist(A1,A2,delta,cap,DistanceType)

INPUTS:

A1,A2: Binary undirected adjacency matrices (possibly with different size)

delta: discretization interval for discrete distributions (0<delta<1)

cap: upper bound for discrete distributions (0<cap<=1, cap>>delta)

DistanceType: {'D','C','P','SUM','CP','DC','DP','DCP'}

[Note: To speed up computations, no check is performed on the correctness and consistency of the inputs.]

OUTPUT:   

distance: distance between networks A1, A2

Example of usage (the two nets are distributed in the folder "networks"):

load('net_SFBA_n1000_d000_1.mat'); A1=A; %loading A1

load('net_GEO_n2000_d000_10.mat'); A2=A; %loading A2

distance=EgoDist(A1,A2,0.01,1,'DCP')

distance =

   672.5809
   
 
