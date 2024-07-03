import math
import pandas as pd
from datetime import datetime, timedelta
from itertools import combinations, takewhile, chain
from functools import cache
from typing import Callable
from scipy.special import comb
from pathlib import Path
import os.path

INPUT = "resources\\matroid_sets\\matroids09_bases"
UPDATE_FREQUENCY = timedelta(seconds=3)
CUTOFF = (2,2)

def timeprint(string: str):
    print(f"[{datetime.now()}] {string}")

@cache
def base_indices(R: int, N: int) -> dict[tuple[int], int]:
    """Returns a dictionary, which tells for each basis its position in the lexicographic ordering."""
    return { base: i for i, base in enumerate(combinations(range(N), R)) }

def read_MayhewRoyle_representatives(verbose=True) -> pd.DataFrame:
    global last_time
    read_data = []
    with open(INPUT, "r") as f:
        if verbose: timeprint(f"Opening {INPUT} ...")
        last_time = datetime.now()
        for i, line in enumerate(f.readlines()):
            split = line.split()
            read_data.append({
                'id': split[0],
                'N': int(split[1]),
                'R': int(split[2]),
                '#bases': int(split[3]),
                'bases': tuple(tuple(int(c) for c in basis[1:-1]) for basis in split[4:])
            })
            if read_data[-1]['N'] == CUTOFF[1] and read_data[-1]['R'] == CUTOFF[0]:
                read_data.pop()
                break # Early cutoff, because we don't really care about (4,9) and above.
            if i % 1000 != 999 or not verbose: continue
            now = datetime.now()
            if now - last_time > UPDATE_FREQUENCY:
                timeprint(f"Processed {i+1} lines.")
                last_time = now
    if verbose: timeprint("Reading finished, conversion into pandas DataFrame starting...")
    return pd.DataFrame(read_data)

def calculate_char_vectors(database: pd.DataFrame, verbose=True):
    """Take a database produced by `read_MayhewRoyle_representatives`,
    and append the additional column `char_vector`, where the
    characteristic vector of the set of bases is added for each 
    matroid (bases listed in the lexicographic order)."""
    if verbose: timeprint("Calculating characteristic vectors...")
    def char_vector(line):
        bases = [False] * comb(line['N'], line['R'], exact=True)
        for base in line['bases']:
            bases[base_indices(line['R'], line['N'])[base]] = True
        return ''.join('1' if t else '0' for t in bases)
    database['char_vector'] = database.apply(char_vector, axis=1)

def iterate_over_matroid_batches(database: pd.DataFrame, verbose=True):
    """Groups the `database` by `R`, `N`, and `#bases`, and returns these
    groups one by one, together with the above three parameters. If `verbose`,
    then also prints status updates of these parameters if enough time has 
    passed."""
    gd = database.groupby(by=['R', 'N', '#bases'])
    for R in range(1, 10):
        for N in range(R, 10):
            if R == CUTOFF[0] and N == CUTOFF[1]: return
            if verbose: timeprint(f"R = {R}, N = {N}")
            last_time = datetime.now()
            for b in range(1, comb(N,R,exact=True)+1):
                if (R, N, b) not in gd.groups: yield(R, N, b, database.iloc[0:0])
                else: yield (R, N, b, gd.get_group((R, N, b)))
                if not verbose: continue
                now = datetime.now()
                if now - last_time > UPDATE_FREQUENCY:
                    last_time = now
                    timeprint(f"R = {R}, N = {N} ...bases: {b}")

def write_char_vectors_to_files(database: pd.DataFrame, verbose=True):
    """Writes the characterstic vectors of all rank `R` matroids
    on `N` elements and `B` bases, one vector per line, as a
    `0-1` sequence, into the file `resources\\matroid_sets\\r{R}n{N}\\{b}_bases_representatives.txt`."""
    if verbose: timeprint("Writing characteristic vectors to files...")
    for R, N, b, matroids in iterate_over_matroid_batches(database):
        path = Path(f"resources\\matroid_sets\\r{R}n{N}\\{b}_bases_representatives.txt")
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            f.writelines(f"{idx}:{cv}\n" for idx, cv in zip(matroids['id'],matroids["char_vector"]))
    
def heaps_algorithm_with_action(number_of_elements: int, action, element, k: int|None = None):
    """Take the list `[1,...,number_of_elements]` and performs a sequence
    of `number_of_elements!-1` transpositions on it, so that each permutation
    appears exactly once. Given moreover is an `element` of a set equipped 
    with a (left) `action` of the symmetric group on `number_of_elements` many 
    elements: the first argument of `action` is the transposition `(i,j)` (`i<j`)
    which should act on the second argument. Every time Heap's algorithm performs 
    a transposition, that transposition acts on the current value of `element`,
    and updates it. The result is yielded every time."""
    if k == None: k = number_of_elements
    if k <= 1: 
        yield element
        return
    for e in heaps_algorithm_with_action(number_of_elements, action, element, k-1):
        element = e # this is to ensure that e is "passed by reference"
        yield e
    for i in range(k-1):
        if k % 2 == 0:
            element = action((i, k-1), element)
        else:
            element = action((0, k-1), element)
        for e in heaps_algorithm_with_action(number_of_elements, action, element, k-1):
            element = e # this is to ensure that e is "passed by reference"
            yield e

@cache
def base_permutations(R: int, N: int) -> list[list[tuple[int]]]:
    """`base_permutations(R,N)[i1][i2][t]` tells the index (in the lexicographic
    ordering) of the basis obtained from the `t`th basis by exchanging the labels
    of the `i1`th and `i2`th elements."""
    return [ 
        [
            tuple(
                base_indices(R, N)[tuple(sorted(t))] 
                for t in combinations(chain(range(i1),(i2,),range(i1+1,i2),(i1,),range(i2+1,N)), R)
            ) if i1 < i2 else None
            for i2 in range(N)
        ]
        for i1 in range(N)    
    ]

def generate_permuted_matroids(R: int, N: int, matroid_iterator):
    """Iterate over all matroids of rank `R` on `N` elements obtained by
    permuting the elements of the representatives of equivalence classes of
    matroids (given as characteristic vectors) in `matroid_iterator`."""
    for matroid in matroid_iterator:
        yield from ("".join(cv) for cv in set(
            tuple(l) for l in heaps_algorithm_with_action(
                N,
                lambda transposition, cv: [cv[base_permutations(R, N)[transposition[0]][transposition[1]][t]] for t in range(len(cv))],
                matroid
            )
        ))

def generate_all_matroids(database: pd.DataFrame, verbose=True):
    """Generates and writes to files all matroids whose isomorphism class
    is represented by some element of `database`. The set of resulting
    rank `R` matroids on `N` elements and with `b` basis will be written to
    the file `resources\\matroid_sets\\r{R}n{N}\\{b}_bases_all.txt`."""
    if verbose: timeprint("Generating all matroids...")
    for R, N, b, matroids in iterate_over_matroid_batches(database, verbose):
        path = Path(f"resources\\matroid_sets\\r{R}n{N}\\{b}_bases_all.txt")
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, 'w') as f:
            for cv in generate_permuted_matroids(R, N, matroids['char_vector']):
                f.write(cv+"\n")


if __name__ == "__main__":
    database = read_MayhewRoyle_representatives()
    calculate_char_vectors(database)
    write_char_vectors_to_files(database)
    generate_all_matroids(database)
    timeprint("Program complete.")