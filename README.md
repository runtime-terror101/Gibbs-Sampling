# Gibbs-Sampling
This python script <br>
* Generates 100 random 1 kilobase sequences (assuming they are upstream regions of some genes). <br>

* Randomly generates a motif of length 10 bases.<br>

* Plants the motif after introducing 0,1, or 2 random mutations (uniformly distributed). <br>

* Then implements the *Gibbs Sampler* algorithm which is used to find common motifs in DNA sequences; to identify the motif locations and the consensus motif. This probablistic search algorithm runs in polynomial time, which is an improvement over the brute force algorithm's exponential running time.<br>
