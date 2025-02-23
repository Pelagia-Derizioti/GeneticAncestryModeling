import numpy as np
from typing import List
from individual import Individual
from chrom_break_pos import get_chrom_break_pos
from scipy.stats import vonmises

def get_backward_population(pop_size: int, num_gen: int, chrom_lengths: List[int], recomb_rate: float) -> List[List[Individual]]:
    chrom_break_pos = get_chrom_break_pos(chrom_lengths)
    pop = [[] for i in range(num_gen)]

    for i in range(num_gen - 1):
        for j in range(pop_size):
            pop[i].append(Individual())
            pop[i][j].set_parents(pop_size)
            pop[i][j].set_break_pos(chrom_break_pos, recomb_rate)
    
    for j in range(pop_size):
        pop[num_gen-1].append(Individual())

    return pop

def get_forward_population(pop_size: int, num_gen: int, chrom_lengths: List[int], recomb_rate: float) -> List[List[Individual]]:
    chrom_break_pos = get_chrom_break_pos(chrom_lengths)
    pop = [[] for i in range(num_gen)]

    for j in range(pop_size):
        pop[0].append(Individual())

    for i in range(1, num_gen):
        for j in range(pop_size):
            pop[i].append(Individual())
            mom_id, dad_id = pop[i][j].set_parents(pop_size)
            pop[i][j].set_break_pos(chrom_break_pos, recomb_rate)
            pop[i - 1][mom_id].children.append(j)
            pop[i - 1][dad_id].children.append(j)

    return pop

def get_forward_population_gender_based_monoamorous_couples(pop_size: int, num_gen: int, chrom_lengths: List[int], recomb_rate: float, kappa_parameter: int, avoid_relatives: str=None) -> List[List[Individual]]:
    chrom_break_pos = get_chrom_break_pos(chrom_lengths)
    mu = 0  # mean direction (in radians)
    kappa = 2  # concentration parameter
    # High kappa (Large Concentration):
    # When k is high, the distribution becomes sharply concentrated around the mean direction mu.
    # It indicates that the values are tightly clustered around mu meaning there is less dispersion. In other words, most of the data points are close to mu
    # Visually, it looks like a narrow, tall peak centered around mu on the circle.
    # Low kappa  (Small Concentration):
    # When kappa is low, the distribution is more spread out around the circle, meaning there is more dispersion.
    # A very low kappa (close to 0) makes the distribution nearly uniform, indicating that the values are almost equally likely to occur at any point on the circle, with little preference for the mean direction mu.
    # Visually, this looks like a broader, flatter distribution centered around mu.
    # High mu (Changing the Mean Direction)
    # A high or low value of mu simply shifts the "center" or peak of the distribution around the circle.
    # In practical terms, mu  determines the angle (or direction) around which the data is centered.
    # Low mu (Different Mean Direction):
    # If mu is set to a different value (e.g., a smaller angle), the distribution will shift its peak to this new angle.
    # The relative spread and concentration of the data are still governed by kappa mu simply determines where the "center" is on the circle.
    # 
    # SO 0 kappa eliminates totally the mu!!!! (perfect uniform on circle) 

    dist_param=2
    # WE MUST SET THE dist_param
    # The dist_param is the kappa parameter of the von mises distribution which will select a random male partner. 
    # Has this to be the same with the kappa of the von mises distibution which allow children to select moms? I do not know.
    # Has this to be different from the von mises distibution which allow children to select moms? I do not know.
    # Either way how much should they be?
    
    pop = [[] for i in range(num_gen)]
    mom_time_list_per_gen=[]
    mom_time_list_per_gen_count=[]
    mom_time_all=[]
 
    print('[DEBUG] : GEN 0')

    for j in range(pop_size):
        print(f'[DEBUG] : ind == {j}')
        ind = Individual()
        ind.location = np.random.uniform(0, 2 * np.pi)  # Assign location uniformly
        pop[0].append(ind)

    mom_time_list_curr_gen=[]
    print(f'[DEBUG] : GEN 1')
    for j in range(pop_size):
        print(f'[DEBUG] : ind == {j}')
        ind = Individual()
        pop[1].append(ind)
        mom_id, dad_id , mom_time= pop[1][j].set_parents_straight_monoamorous_couple(pop_size, pop[0],None, kappa_parameter,avoid_relatives)
        if mom_time > -1 :
            mom_time_list_curr_gen.append(mom_time)
            mom_time_all.append(mom_time)
        pop[1][j].location = Individual.uniform_location_on_circle(pop[0][mom_id],pop[0][dad_id],pop[1])
        pop[1][j].set_break_pos(chrom_break_pos, recomb_rate)
        pop[0][mom_id].children.append(j)
        pop[0][dad_id].children.append(j)

    mom_time_list_per_gen.append(np.average(mom_time_list_curr_gen))
    mom_time_list_per_gen_count.append(len(mom_time_list_curr_gen))
    for i in range(2, num_gen):
        print(f'[DEBUG] : GEN {i}')

        mom_time_list_curr_gen=[]
        for j in range(pop_size):
            print(f'[DEBUG] : ind == {j}')


            ind = Individual()
            pop[i].append(ind)
            mom_id, dad_id , mom_time = pop[i][j].set_parents_straight_monoamorous_couple(pop_size, pop[i-1],pop[i-2], kappa_parameter,avoid_relatives)
            if mom_time > -1 :
                mom_time_list_curr_gen.append(mom_time)
                mom_time_all.append(mom_time)
            pop[i][j].location = Individual.uniform_location_on_circle(pop[i - 1][mom_id],pop[i - 1][dad_id],pop[i])
            pop[i][j].set_break_pos(chrom_break_pos, recomb_rate)
            pop[i - 1][mom_id].children.append(j)
            pop[i - 1][dad_id].children.append(j)
        mom_time_list_per_gen.append(np.average(mom_time_list_curr_gen))
        mom_time_list_per_gen_count.append(len(mom_time_list_curr_gen))

    return pop, mom_time_list_per_gen, mom_time_list_per_gen_count, np.average(mom_time_all)




