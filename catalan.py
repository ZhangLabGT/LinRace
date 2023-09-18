import math
import pandas as pd
import functools


@functools.cache
def num_bintrees(n):
    if n == 0:
        return 1
    total = 0
    for i in range(0, n):
        total += num_bintrees(i) * num_bintrees(n - 1 - i)
    return total


# total number of binary trees and cell arrangements
N = [i for i in range(10)]
counts = [math.factorial(n) * num_bintrees(n) for n in N]

df = pd.DataFrame({'N': N, 'count': counts})
df.to_csv('num_possibilities.csv')
