### testing

library(igraph)
library(BoolNet)
library(Matrix)


mutate_graph_edges_and_nodes <- function(adjmat) {
  adjmat <- as.matrix(adjmat)  
  
  num_nodes <- nrow(adjmat)
  edges <- which(adjmat == 1, arr.ind = TRUE)
  
  if (nrow(edges) == 0) return(adjmat)
  
  num_mutations <- sample(1:min(3, nrow(edges)), 1)
  edges_to_switch <- sample(1:nrow(edges), num_mutations)
  
  for (i in edges_to_switch) {
    if (runif(1) > 0.5) {
      adjmat[sample(1:num_nodes, 1), edges[i, 2]] <- 1
    } else {
      adjmat[edges[i, 1], sample(setdiff(1:num_nodes, edges[i, 1]), 1)] <- 1
    }
    adjmat[edges[i, 1], edges[i, 2]] <- 0
  }
  
  return(adjmat)
}


cluster_mapping <- list(
  "1000000000000000" = "1  A, B, C, the empty graph",
  "0100000000000000" = "2  A -> B, C, the graph with a single directed edge",
  "0010000000000000" = "3  A <-> B, C, the graph with a mutual connection between two vertices",
  "0001000000000000" = "4  A <- B -> C, the out-star",
  "0000100000000000" = "5  A -> B <- C, the in-star",
  "0000010000000000" = "6  A -> B -> C, directed line",
  "0000001000000000" = "7  A <-> B <- C",
  "0000000100000000" = "8  A <-> B -> C",
  "0000000010000000" = "9  A -> B <- C, A -> C",
  "0000000001000000" = "10 A <- B <- C, A -> C",
  "0000000000100000" = "11 A <-> B <-> C",
  "0000000000010000" = "12 A <- B -> C, A <-> C",
  "0000000000001000" = "13 A -> B <- C, A <-> C",
  "0000000000000100" = "14 A -> B -> C, A <-> C",
  "0000000000000010" = "15 A -> B <-> C, A <-> C",
  "0000000000000001" = "16 A <-> B <-> C, A <-> C, the complete graph"
)


combined_graphs <- list()
clustered_results <- list()


for (i in 1:1000) {
  adjmat <- as.matrix(
    as_adjacency_matrix(plotNetworkWiring(generateRandomNKNetwork(3, sample(c(1, 2, 3), 1))))
  )
  
  adjmat[adjmat > 1] <- 1
  
  if (!any(sapply(combined_graphs, function(x) identical(x, adjmat)))) {
    combined_graphs[[length(combined_graphs) + 1]] <- adjmat
  }
  
  max_attempts <- 50
  mutation_count <- 0
  
  while (mutation_count < max_attempts) {
    mutation_count <- mutation_count + 1
    mutated_adjmat <- mutate_graph_edges_and_nodes(combined_graphs[[sample(1:length(combined_graphs), 1)]])
    
    if (!any(sapply(combined_graphs, function(x) identical(x, mutated_adjmat)))) {
      combined_graphs[[length(combined_graphs) + 1]] <- mutated_adjmat
      mutation_count <- 0  
      
      
      graph_representation <- graph_from_adjacency_matrix(mutated_adjmat, mode = "directed")
      
      
      if (length(E(graph_representation)) == 0) {
        #print(paste("Graph", length(combined_graphs), "has no edges. Skipping triad census."))
        next  
      }
      
      
      cluster_info <- tryCatch(
        triad_census(graph_representation),
        error = function(e) {
          #print(paste("Skipping triad census for graph:", length(combined_graphs), "due to error:", e$message))
          return(NULL)  
        }
      )
      
      
      if (!is.null(cluster_info)) {
        # Convert the triad census result to binary string
        cluster_key <- paste(as.integer(cluster_info > 0), collapse = "")  # Create binary string (1 if count > 0, else 0)
  
        
        if (cluster_key %in% names(cluster_mapping)) {
          #print("Cluster Key Found in Mapping")
          
          header <- paste("Cluster", cluster_mapping[[cluster_key]])  
          #print(paste("Header Retrieved:", header))  
        #} else {
          #print("Cluster Key NOT Found in Mapping")
          #print(paste("Generated Key: ", cluster_key))  
          #header <- paste("Cluster:", cluster_key)  
        }
        
        # Assign graph to its cluster in clustered_results
        if (!cluster_key %in% names(clustered_results)) {
          clustered_results[[cluster_key]] <- list()
        }
        clustered_results[[cluster_key]] <- append(clustered_results[[cluster_key]], list(mutated_adjmat))
      }
    }
  }
}

#save in .txt file

file_conn <- file("clustered_graphs.txt", open = "w")


for (cluster_key in names(clustered_results)) {
  
  if (cluster_key %in% names(cluster_mapping)) {
    header <- paste("Cluster", cluster_mapping[[cluster_key]])
  #} else {
    #header <- paste("Cluster:", cluster_key) 
  }
  
  
  writeLines(header, file_conn)
  
  
  for (adjmat in clustered_results[[cluster_key]]) {
    write.table(adjmat, file = file_conn, row.names = FALSE, col.names = FALSE, append = TRUE)
    writeLines("\n", file_conn)  
  }
}

close(file_conn)


print(paste("Total number of unique graphs in combined_graphs:", length(combined_graphs)))
print("Clustered graphs successfully saved to 'clustered_graphs.txt'")
