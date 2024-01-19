**Title: Novel Graph Clustering Methodology Based on Maximal Clique**

**Abstract:**
This paper introduces a new graph clustering methodology utilizing Maximal Clique, analyzed against Ground Truth using the F-measure. It focuses on identifying Maximal Cliques in gene datasets, refining these into initial clusters by removing overlaps, and specifically addressing the clustering of single-node groups.

**1. Introduction:**
The paper begins by defining graphs and clustering, emphasizing graph clustering techniques. It discusses the significance of edge structures in graph clustering and introduces two primary forms: point clustering and edge clustering. The paper references key works in the field, including Density-Based Method for graph clustering.

**2. Clustering Algorithms:**
- **2.2 M1-Graph Algorithm:** 
  - Describes the M1-Graph algorithm, which uses Maximal Cliques for initial clustering and then merges unclustered nodes using a tree-like structure based on edge Jaccard values. 
  - Highlights the algorithm's tendency to form large, dense clusters.

- **2.3 Limitations of M1-Graph:** 
  - Discusses the M1-Graph's limitation in its exclusion of internal cluster nodes in Jaccard calculations, leading to biased, oversized clusters.

- **2.4 M2-Graph Algorithm:** 
  - Introduces the M2-Graph algorithm as an improvement over M1-Graph, including internal cluster nodes in Jaccard value calculations to prevent large cluster formations.
  - Shows more evenly distributed clusters without data loss.

- **2.5 Limitations of M2-Graph:** 
  - Notes a slight decrease in F1-score compared to M1-Graph but an overall improvement in cluster distribution.

- **2.6 M3-Graph Algorithm:** 
  - Presents M3-Graph, a further attempt to refine clustering by stopping the process when encountering clusters larger than two nodes.
  - Aims for more meaningful clustering but results in a lower F1-score.

- **2.7 Limitations of M3-Graph:** 
  - Acknowledges the lower effectiveness of M3-Graph in capturing significant clustering patterns compared to M1 and M2 algorithms.

**3. Results and Interpretation of Graph Clustering Algorithms:**
- Compares the F1-scores, cluster numbers, and maximum cluster sizes of the three algorithms.
- Concludes that M2-Graph is the most effective, balancing cluster distribution and meaningfulness.

**4. Conclusion:**
- Summarizes the findings, highlighting the M2-Graph algorithm's balance in avoiding large clusters while maintaining meaningful clustering.
- Discusses the potential of these methodologies in gene dataset clustering, emphasizing the importance of cluster density and connectivity.
