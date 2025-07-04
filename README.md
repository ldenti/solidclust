# isONclust3.5
isONclust3.5 is a method for the de novo clustering of third-generation transcriptomic datasets. isONclust3.5 is heavily inspired by [isONclust3](https://github.com/aljpetri/isONclust3). Although they both follow the same greedy methodology, isONclust3.5 introduces the notion of _solid_ high confidence minimizer. A high confidence minimizer (as introduced in isONclust3) is considered during the clustering step iff its count in the cluster under investigation is a tangible fraction of all the counts associated to that minimizer among all already created clusters. In other words, before computing the Jaccard similarity between a read and a cluster, we prefilter the high-confidence minimizers based on the past history of the clustering process. This allows a clustering more robust to sequencing errors.


``` sh
git clone https://github.com/ldenti/isonclust3.5.git
cd isonclust3.5
make
./build/isonclust -i example/small.fq -o example/small.clustering.out -t 0.5 --weighted 0.1
```
