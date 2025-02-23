import argparse
import random
import numpy as np
import time
import sys
import json
import os
from datetime import datetime
from tqdm import tqdm
from auxiliary_functions import read_file, sum_list_of_int, print_json
from population import get_forward_population_gender_based_monoamorous_couples
from genealogical_ancestors import get_genealogical_ancestors
from genetic_ancestors import get_genetic_ancestors
from genealogical_descendants import get_genealogical_descendants
from genetic_descentants import get_genetic_descendants
from ibd_analysis import calculate_ibd_proportion


def main():
    # Command-line argument parser
    parser = argparse.ArgumentParser(description='Simulation Parameters')
    parser.add_argument('-p', '--pop_size', type=int, default=100, help='Population size')
    parser.add_argument('-g', '--num_gen', type=int, default=50, help='Number of generations')
    parser.add_argument('-r', '--recomb_rate', type=float, default=1e-8, help='Recombination rate')
    parser.add_argument('-f', '--chrom_len_file', type=str, default='', help='Filepath of chromosome lengths file')
    parser.add_argument('-s', '--seed', type=int, default=0, help='Seed for random number generation')
    parser.add_argument('-k', '--kappa_parameter', type=float, default=1e-08, help='Kappa parameter')
    parser.add_argument('-a', '--avoid_relatives', type=str, default='', help='Avoid siblings or first cousins or both or leave it empty for not inbreeding avoidance. Possible values:"siblings", "cousins", "siblingscousins"')


    args = parser.parse_args()

    # Set seed for random number generator
    if args.seed != 0:
        random.seed(args.seed)
        np.random.seed(args.seed)

    # Read chromosome lengths from file (or use default)
    chrom_lengths = [3000000000]
    if args.chrom_len_file:
        if len(args.chrom_len_file):
            print("Read file", file=sys.stderr)
            chrom_lengths = read_file(args.chrom_len_file)

    # Initialize data dictionary
    data = {
        'pop_size': args.pop_size,
        'num_gen': args.num_gen,
        'recomb_rate': args.recomb_rate,
        'seed': args.seed,
        'chrom_lengths': chrom_lengths,
        'seq_len': sum_list_of_int(chrom_lengths),
        'kappa_parameter': args.kappa_parameter,
        'avoid_relatives': args.avoid_relatives,
    }
    print(f'script started at: {datetime.now().strftime("%Y_%m_%d_%H_%M_%S")}', file=sys.stderr)
    print(data, file=sys.stderr)
    # Simulate backward-in-time population
    #print("Compute get_backward_population", file=sys.stderr )
    #start_time = time.perf_counter()
    #backward_pop = get_backward_population(args.pop_size, args.num_gen, chrom_lengths, args.recomb_rate)
    #passed_time = time.perf_counter() - start_time
    #print(f"Finished get_backward_population in {str(round(passed_time, 2))} seconds", file=sys.stderr)
    # Simulate forward-in-time population
    
    
    print("Compute get_forward_population", file=sys.stderr)
    start_time = time.perf_counter()
    forward_pop,mom_time_list_per_gen, mom_time_list_per_gen_count, average_mom_time_all = get_forward_population_gender_based_monoamorous_couples(args.pop_size, args.num_gen, chrom_lengths, args.recomb_rate, args.kappa_parameter,args.avoid_relatives)
    passed_time = time.perf_counter() - start_time
    print(f"Finished get_forward_population in {str(round(passed_time, 2))} seconds", file=sys.stderr)
    data['runtime_create_population_with_forward_method'] = passed_time
    reversed_pop = list(reversed(forward_pop))
    
    data["Average_time_mom_finds_a_partner_per_gen"] = mom_time_list_per_gen
    data["Number_of_moms_searched_for_a_partner_per_gen"] = mom_time_list_per_gen_count
    data["average_mom_time_all"] = average_mom_time_all
    print("Compute get_genealogical_ancestors", file=sys.stderr)
    start_time = time.perf_counter()
    genealogical_ancestors, TMRCA, IAP = get_genealogical_ancestors(reversed_pop)

    data['genealogical_ancestors'] = genealogical_ancestors
    data['TMRCA'] = TMRCA
    data['IAP'] = IAP
    passed_time = time.perf_counter() - start_time
    data['runtime_genealogical_ancestors'] = passed_time
    print(f"Finished get_genealogical_ancestors in {str(round(passed_time, 2))} seconds", file=sys.stderr)

    print(" get_genetic_ancestors", file=sys.stderr)
    start_time = time.perf_counter()
    genetic_ancestors = [float(0) for i in range(args.num_gen)]
    for i in tqdm(range(args.num_gen), bar_format='{n_fmt}/{total_fmt} | {l_bar}{bar} | elapsed: {elapsed} - remaining: {remaining} | {rate_fmt}{postfix}'):
        genetic_ancestors[i] = get_genetic_ancestors(i,reversed_pop, chrom_lengths)
    data['genetic_ancestors'] = genetic_ancestors
    passed_time = time.perf_counter() - start_time
    data['runtime_genetic_ancestors'] = passed_time
    print(f"Finished get_genetic_ancestors in {str(round(passed_time, 2))} seconds", file=sys.stderr)
    

  

    # Get number of genealogical descendants
    print("Compute get_genealogical_descendants", file=sys.stderr)
    start_time = time.perf_counter()
    genealogical_descendants = get_genealogical_descendants(forward_pop)
    data['genealogical_descendants'] = genealogical_descendants
    passed_time = time.perf_counter() - start_time
    data['runtime_genealogical_descendants'] = passed_time
    print(f"Finished get_genealogical_descendants in {str(round(passed_time, 2))} seconds", file=sys.stderr)

        # Simulate recombinations and get number of genetic descendants, segment count, ROH frequencies and lengths
    print("Compute get_genetic_descendants", file=sys.stderr)
    start_time = time.perf_counter()
    genetic_descendants, segment_count, segment_len, roh_freq, roh_len = get_genetic_descendants(forward_pop, chrom_lengths)
    data['genetic_descendants'] = genetic_descendants
    data['segment_count'] = segment_count
    data['segment_len'] = segment_len
    data['roh_freq'] = roh_freq
    data['roh_len'] = roh_len
    passed_time = time.perf_counter() - start_time
    data['runtime_genetic_descendants'] = passed_time
    print(f"Finished get_genetic_descendants in {str(round(passed_time, 2))} seconds", file=sys.stderr)
    
    
    # Assuming `pop` (population data) and `chrom_lengths` are already created
    ibd_proportions = calculate_ibd_proportion(forward_pop, chrom_lengths)

    print("IBD Proportions:", ibd_proportions)
    data['ibd_proportions'] = ibd_proportions
    

    # Print the results as JSON
    print_json(data)
    file_path = "/home/people/s222822/thesis/res/"
    file_current_time = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    slurm_job_id = os.getenv("SLURM_JOB_ID", "no_job_id")
    file_name = f"results_{slurm_job_id}_job_{file_current_time}.json"
    with open(file_path+file_name, 'w') as json_file:
        json.dump(data, json_file, indent=4, sort_keys=True)


if __name__ == '__main__':
    main()

