# Efficient Algorithms for Density Decomposition on Large Static and Dynamic Graphs

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

The datasets need to be preprocessed, since our paper only focuses on undirected simple graph. Some datasets may be directed graphs, or contain loops, or contain multiple edges. A preprocessing is needed to convert such datasets into undirected simple graph by ignoring the direction of edges, deleting loops, and considering multiple edges as one edge. The preprocessing can be done with the simple program “preprocess.cpp” in the “Graph” folder. We uploaded a preprocessed dataset as an example, stored in **"example_input.txt"**.

# Compile

```
g++ main.cpp -o main -std=c++11 -O3
```

# Usage

Datasets need to be stored in the “Graphs” folder. In “dataset.txt”, enter the graph name as our example shown. To run the program, run the following command. 

```
./main
```
