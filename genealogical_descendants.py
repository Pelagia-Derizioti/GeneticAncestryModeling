from typing import List
from auxiliary_functions import new_matrix, mean_non_zero
from individual import Individual

def get_genealogical_descendants(pop: List[List[Individual]]) -> List[float]:
    num_gen = len(pop)
    pop_size = len(pop[0])

    num_descendants = new_matrix(num_gen, pop_size)

    for j in range(pop_size):
        q = set()
        q.add(j)

        for i in range(num_gen):
            num_descendants[i][j] = len(q)
            new_q = set()
            for k in q:
                for child in pop[i][k].children:
                    new_q.add(child)
            q = new_q

    mean_num_descendants = [float(0) for i in range(num_gen)]
    for i in range(num_gen):
        mean_num_descendants[i] = mean_non_zero(num_descendants[i])

    return mean_num_descendants


