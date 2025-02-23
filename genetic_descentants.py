from typing import List, Tuple
import numpy as np
from individual import Individual
from chrom_pair import ChromPair
from auxiliary_functions import sum_list_of_int, mean_list_of_int, mean_non_zero

def get_genetic_descendants(pop: List[List[Individual]], chrom_lengths: List[int]) -> Tuple[List[float], List[float], List[int], List[float], List[int]]:
    num_gen = len(pop)
    pop_size = len(pop[0])
    seq_len = sum_list_of_int(chrom_lengths)

    mean_num_descendants = [float(0.0) for i in range(num_gen)]
    mean_num_descendants[0] = 1

    mean_segment_count = [float(0.0) for i in range(num_gen)]
    mean_segment_count[0] = float(len(chrom_lengths))

    segment_len = [int(0) for i in range(num_gen)]
    segment_len[0] = int(mean_list_of_int(chrom_lengths))

    roh_freq = [float(0.0) for i in range(num_gen)]
    roh_freq[0] = 0

    mean_roh_len = [int(0) for i in range(num_gen)]
    mean_roh_len[0] = 0

    # Run recombination simulation
    prev_gen = [ChromPair() for i in range(pop_size)]
    for j in range(pop_size):
        for k in range(len(chrom_lengths)):
            prev_gen[j].chrom_pair[0].add_segment(j, chrom_lengths[k])
            prev_gen[j].chrom_pair[1].add_segment(j, chrom_lengths[k])

    for i in range(1, num_gen):
        curr_gen = [ChromPair() for j in range(pop_size)]
        num_descendants = [int(0) for j in range(pop_size)]
        num_segments = 0
        roh_lengths = [int(0) for j in range(pop_size)]
        roh_counts = [int(0) for j in range(pop_size)]
        for j in range(pop_size):
            curr_gen[j].chrom_pair[0] = prev_gen[pop[i][j].mom_id].get_chrom(pop[i][j].mom_break_pos)
            curr_gen[j].chrom_pair[1] = prev_gen[pop[i][j].dad_id].get_chrom(pop[i][j].dad_break_pos)

            # Get number of descendants
            ancestors = set()
            for k in range(len(curr_gen[j].chrom_pair[0].ids)):
                ancestors.add(curr_gen[j].chrom_pair[0].ids[k])
            for k in range(len(curr_gen[j].chrom_pair[1].ids)):
                ancestors.add(curr_gen[j].chrom_pair[1].ids[k])

            for k in ancestors:
                num_descendants[k] += 1

            # Get number of segments
            num_segments += (len(curr_gen[j].chrom_pair[0].ids) + len(curr_gen[j].chrom_pair[1].ids))

            # Get lengths of ROH
            roh_lengths[j], roh_counts[j] = curr_gen[j].get_roh(seq_len)

        # Get average number of descendants
        mean_num_descendants[i] = mean_non_zero(num_descendants)

        # Get average segment count
        mean_segment_count[i] = float(num_segments) / float(pop_size)

        # Get average segment length
        segment_len[i] = 2 * seq_len * pop_size // num_segments

        # Get ROH frequency
        roh_freq[i] = float(mean_list_of_int(roh_lengths)) / float(seq_len)
        s = sum_list_of_int(roh_counts)
        if s > 0:
            mean_roh_len[i] = sum_list_of_int(roh_lengths) // s
        else:
            mean_roh_len[i] = 0

        prev_gen = curr_gen

    return mean_num_descendants, mean_segment_count, segment_len, roh_freq, mean_roh_len

