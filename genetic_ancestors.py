from typing import List
from chrom_pair import ChromPair
from individual import Individual
import numpy as np
from auxiliary_functions import mean_list_of_int


def calculate_chromosomes_for_generation(
    pop: List[List[Individual]], prev_gen: List[ChromPair], chrom_lengths: List[int], gen_index: int, pop_size: int
) -> List[ChromPair]:
    
 #Updates chromosome data for the current generation based on parental information.
    
    curr_gen = [ChromPair() for _ in range(pop_size)]
    for j in range(pop_size):
        mom_id, dad_id = pop[gen_index][j].mom_id, pop[gen_index][j].dad_id
        mom_break_pos, dad_break_pos = pop[gen_index][j].mom_break_pos, pop[gen_index][j].dad_break_pos
        curr_gen[j].chrom_pair[0] = prev_gen[mom_id].get_chrom(mom_break_pos)
        curr_gen[j].chrom_pair[1] = prev_gen[dad_id].get_chrom(dad_break_pos)
    return curr_gen


def get_unique_ancestors(prev_gen: List[ChromPair], pop_size: int) -> np.ndarray:
    
    #Extracts unique ancestors for each individual in the population.
   
    num_ancestors = np.zeros(pop_size, dtype=int)
    for j in range(pop_size):
        ancestors = set(prev_gen[j].chrom_pair[0].ids + prev_gen[j].chrom_pair[1].ids)
        num_ancestors[j] = len(ancestors)
    return num_ancestors


def get_genetic_ancestors(gen_minus: int, pop: List[List[Individual]], chrom_lengths: List[int]) -> float:
   
    #Calculates the mean number of unique genetic ancestors for a population after a given number of generations.
   
    pop_size = len(pop[0])
    # Initialize chromosome data for the initial generation
    prev_gen = [ChromPair() for _ in range(pop_size)]
    for j in range(pop_size):
        for chrom_length in chrom_lengths:
            prev_gen[j].chrom_pair[0].add_segment(j, chrom_length)
            prev_gen[j].chrom_pair[1].add_segment(j, chrom_length)

    # Backtrack through generations to update chromosome data
    for i in range(gen_minus - 1, -1, -1):
        prev_gen = calculate_chromosomes_for_generation(pop, prev_gen, chrom_lengths, i, pop_size)

    # Calculate the number of unique ancestors for each individual
    num_ancestors = get_unique_ancestors(prev_gen, pop_size)

    # Return the mean number of ancestors
    return mean_list_of_int(num_ancestors)

