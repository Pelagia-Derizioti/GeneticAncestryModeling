import numpy as np
import matplotlib.pyplot as plt
from typing import List
from chrom_pair import ChromPair
from individual import Individual

def calculate_ibd_proportion(pop: List[List[Individual]], chrom_lengths: List[int]) -> List[float]:
    num_gen = len(pop)
    pop_size = len(pop[0])
    seq_len = sum(chrom_lengths)  # Total genome length
    ibd_proportions = []

    # Initialize chromosomes for the first generation
    prev_gen = [ChromPair() for _ in range(pop_size)]
    for j in range(pop_size):
        for chrom_length in chrom_lengths:
            prev_gen[j].chrom_pair[0].add_segment(j, chrom_length)  # Paternal chromosome
            prev_gen[j].chrom_pair[1].add_segment(j, chrom_length)  # Maternal chromosome

    ibd_proportions.append(0.0)

    # Loop through generations starting from Generation 1
    for i in range(1, num_gen):
        total_ibd_length = 0
        curr_gen = [ChromPair() for _ in range(pop_size)]

        # Update chromosomes based on parents
        for j in range(pop_size):
            mom_id, dad_id = pop[i][j].mom_id, pop[i][j].dad_id
            mom_break_pos, dad_break_pos = pop[i][j].mom_break_pos, pop[i][j].dad_break_pos

            curr_gen[j].chrom_pair[0] = prev_gen[mom_id].get_chrom(mom_break_pos)
            curr_gen[j].chrom_pair[1] = prev_gen[dad_id].get_chrom(dad_break_pos)

        # Compute IBD across all individuals in the generation
        for j in range(pop_size):
            for k in range(j + 1, pop_size):  # Avoid duplicate comparisons
                ibd_length = curr_gen[j].get_ibd_length(curr_gen[k])
                total_ibd_length += ibd_length

        print(f"\n--- Generation {i} ---")
        print(f"Total IBD Length: {total_ibd_length}")
        print(f"Population Size: {pop_size}")
        print(f"Sequence Length: {seq_len}")
        print(f"Denominator (pop_size * seq_len): {(pop_size * (pop_size - 1) // 2) * seq_len}")


        # Store the IBD proportion for this generation
        #ibd_proportions.append(total_ibd_length / (pop_size * seq_len))
        ibd_proportions.append(total_ibd_length / ((pop_size * (pop_size - 1) // 2) * seq_len))

        prev_gen = curr_gen  # Update for next generation

    return ibd_proportions
