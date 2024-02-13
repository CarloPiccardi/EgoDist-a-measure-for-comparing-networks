# EgoDist: a measure for comparing networks

EgoDist is a Matlab function that computes a series of measures that quantify the dissimilarity between two undirected, unweighted networks.
It is defined by processing the distributions, on the network, of three egonet features, namely the degree, the clustering coefficient and the egonet persistence.
The method exploits the statistics of the three features to define one- or multi-dimensional distribution functions, which are then compared to define a distance between networks.
It does not require the alignment of the two networks being compared, which can also have different size.

