TreeChopper

TreeChopper clusters tree leaf nodes according to phylogenetic distance.

Algorithm: A graph is constructed from the tree like so:  all leaves are visited, and from each leaf, all neighboring leaves within a specified distance threshold are added to a graph with an edge placed between them.  After building this graph, each edge connecting pairs of nodes is examined and a Jaccard similarity coefficient is computed (see http://www.biomedcentral.com/1741-7007/3/7 for details).  Those edges that loosely connect nodes as defined by this similarity coefficient are removed.  The nodes connected by the remaining edges are clustered by transitive closure (single linkage clustering) and reported as OTUs.  

The minimum phylogenetic distance between clustered nodes, and the minimum similarity coefficient between nodes in the graph are tuneable parameters.

Example usage and sample data are provided in the sample_data/ directory.



INSTALLATION REQUIREMENTS:

1.  Bioperl:  Download and install Bioperl from http:://bioperl.org

2.  slclust:  Install the slclust utility and copy it into the local util/ directory here.  slclust can be obtained here:  https://sourceforge.net/projects/slclust


Questions, comments, etc?  Contact Brian Haas bhaas@broadinstitute.org

