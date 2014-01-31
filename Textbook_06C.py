#!/usr/bin/env python
'''
A solution to a code challenges that accompanies Bioinformatics Algorithms: An Active-Learning Approach by Phillip Compeau & Pavel Pevzner.
The textbook is hosted on Stepic and the problem is listed on ROSALIND under the Textbook Track.

Problem Title: 2-Break Distance Problem
Rosalind ID: 6C
URL: http://rosalind.info/problems/6c/
'''

from collections import defaultdict


def two_break_dist(P, Q):
    '''Returns the 2-Break Distance of Circular Chromosomes P and Q.'''

    # Construct the break point graph of P and Q.
    graph = defaultdict(list)
    for perm_cycle in P+Q:
        L = len(perm_cycle)
        for i in xrange(len(perm_cycle)):
            # Add the edge between consecutive items (both orders since the breakpoint graph is undirected).
            # Note: Modulo L in the higher index for the edge between the last and first elements.
            graph[perm_cycle[i]].append(-1*perm_cycle[(i+1) % L])
            graph[-1*perm_cycle[(i+1) % L]].append(perm_cycle[i])

    # BFS to find the number of connected components in the breakpoint graph.
    component_count = 0
    remaining = set(graph.keys())
    while len(remaining) > 0:
        component_count += 1
        queue = [remaining.pop()]  # Components are cyclic, so starting point is unimportant.
        while queue:
            current = queue.pop(0)
            queue += filter(lambda node: node in remaining, graph.get(current, []))
            remaining -= set(queue)  # Overkill, but it's nice and concise!

    # Theorem: d(P,Q) = blocks(P,Q) - cycles(P,Q)
    return sum(map(len,P)) - component_count


def main():
    '''Main call. Reads, runs, and saves problem specific data.'''

    # Read the input data.
    with open('data/textbook/rosalind_6c.txt') as input_data:
        P, Q = [line.strip().lstrip('(').rstrip(')').split(')(') for line in input_data.readlines()]
        P = [map(int, perm_cycle.split()) for perm_cycle in P]
        Q = [map(int, perm_cycle.split()) for perm_cycle in Q]

    # Get the 2-Break Distance.
    dist = two_break_dist(P, Q)

    # Print and save the answer.
    print str(dist)
    with open('output/textbook/Textbook_06C.txt', 'w') as output_data:
        output_data.write(str(dist))

if __name__ == '__main__':
    main()
