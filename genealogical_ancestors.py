from typing import List, Tuple
from auxiliary_functions import new_matrix, mean_non_zero, mean_list_of_float, max_list_of_int
from individual import Individual
import numpy as np

def get_genealogical_ancestors(pop: List[List[Individual]]) -> Tuple[List[float], int, int]:
    num_gen = len(pop)
    pop_size = len(pop[0])

    num_ancestors = new_matrix(num_gen, pop_size)
    num_desc_in_gen0 = new_matrix(num_gen, pop_size)

    for j in range(pop_size):
        q = set()
        q.add(j)

        for i in range(num_gen - 1):
            num_ancestors[i][j] = len(q)
            new_q = set()
            for k in q:
                new_q.add(pop[i][k].mom_id)
                new_q.add(pop[i][k].dad_id)
                num_desc_in_gen0[i][k] += 1
            q = new_q
        for k in q:
            num_desc_in_gen0[num_gen - 1][k] += 1
        num_ancestors[num_gen - 1][j] = len(q)

    mean_num_ancestors = [float(0.0) for i in range(num_gen)]
    TMRCA = -1
    IAP = -1
    for i in range(num_gen):
        mean_num_ancestors[i] = mean_list_of_float(num_ancestors[i])
        if TMRCA == -1 and (max_list_of_int(num_desc_in_gen0[i]) == pop_size):
            TMRCA = i
        if IAP == -1 and (int(mean_non_zero(num_desc_in_gen0[i])) == pop_size):
            IAP = i

    return mean_num_ancestors, TMRCA, IAP


