![image](https://github.com/Flydragonet/Density_Decomposition_Computation/assets/104469258/87b99a68-c17e-49e3-9561-3573c602f452)# Efficient Algorithms for Density Decomposition on Large Static and Dynamic Graphs

# Datasets

Because some datasets used in the paper are too large to be uploaded to GitHub, we have summarized the download links for the dataset in the table below.

Datasets used for performance studies.

| Dataset | Link |
| --- | --- |
| DBLP | http://www.konect.cc/networks/com-dblp/ |
| Citeseer | http://www.konect.cc/networks/citeseer/ |
| Yahoo | http://www.konect.cc/networks/lasagne-yahoo/ |
| Skitter | http://www.konect.cc/networks/as-skitter/ |
| Weibo | https://networkrepository.com/soc-sinaweibo.php |
| UKlink | http://www.konect.cc/networks/dimacs10-uk-2002/ |
| Twitter | http://www.konect.cc/networks/twitter/ |
| Wiki | http://www.konect.cc/networks/wikipedia_link_en/ |

Datasets used for case studies are uploaded to the “Case Study Datasets” folder. 

# Preprocess

The datasets need to be preprocessed, since our paper only focuses on undirected simple graph. Some datasets may be directed graphs, or contain loops, or contain multiple edges. A preprocessing is needed to convert such datasets into undirected simple graph by ignoring the direction of edges, deleting loops, and considering multiple edges as one edge. The preprocessing can be done with the simple program “preprocess.cpp” in the “Graph” folder. We preprocessed the SK dataset in the paper as an example, stored in "Graphs/sk.txt".

# Compile

```
g++ main.cpp -o main -std=c++11 -O3
```

# Usage of Static Algorithms

Datasets need to be stored in the “Graphs” folder. In “dataset.txt”, enter the graph name as our example shown. To run the program, run the following command. 

```
./main
```

# Usage of Dynamic Algorithms

Our dynamic algorithms take input in the form of unreversible orientations rather than undirected graphs. An unreversible orientation is stored in a file format where each line "from to" represents a directed edge <from, to> (see example in Graphs/sk_unreversible.txt). To compute an unreversible orientation from an undirected graph, one can use the get_unreversible.cpp program located in the Graphs folder. This program implements the proposed INDEGREE+Retest algorithm. The compilation process is the same as that for static algorithms, and it is executed using the command ``./get_unreversible <undirected_graph_address> <output_address>``. Once the unreversible orientation is obtained, it should be placed in the Graphs folder, and then the orientation name should be entered in "dataset.txt" as demonstrated in our example. Finally, the dynamic algorithms can be run by executing the command ``./main``.

Note that for dynamic algorithms, the number of edges to be updated can be changed by the "insert_num" variable, and this number must be more than the number of edges of the whole graph. 
