# `matroid09_bases`

The file `matroid09_bases` contains a list, with exactly one representative of each isomorphism class of ordinary matroids on at most 9 elements. Two matroids are isomorphic if they have the same number of elements, and can be transformed into each other by relabeling the elements.

Each row of the file encodes a matroid. Each line consists of the following, separated by whitespace:
- an identification number of the given matroid
- the number of elements of the matroid
- rank of the matroid
- number of bases
- a sequence of bases, each base listed between single quotes.

The rows are lexicographically ordered by the number of elements and the rank.

**Warning:** due to internet connectivity problems while downloading, the file is incomplete: it is truncated around line 332831/383172(?), at a rank 5 matroid on 9 elements and 116 bases. The missing cases will not be relevant for our computations anyway.

**Source:** https://research-repository.uwa.edu.au/en/datasets/matroids-on-9-elements