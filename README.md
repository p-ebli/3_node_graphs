This R script generates all **unique 3-node directed Boolean network graphs** with self-loops allowed, ensuring each graph has at least two edges and all nodes are connected to at least one other (excluding isolated nodes). It uses **BoolNet** and **igraph** to create and mutate networks, filters out duplicates, and uses **triad_census()** to group the graphs by structural motifs. The resulting graphs and their motif classifications are saved to a text file for further analysis in systems biology, network science, or machine learning.
