
from typing import List, Tuple
import math
import random
import numpy as np
import sympy as sp
import sys
import matplotlib.pyplot as plt
import time 
from datetime import datetime
from scipy.stats import vonmises

class Individual:
    def __init__(self,location_degrees=None) -> None:
        self.mom_id: int = -1
        self.dad_id: int = -1
        self.pair_id: int = -1
        self.mom_break_pos: List[int] = []
        self.dad_break_pos: List[int] = []
        self.children: List[int] = []
        self.location: float = np.radians(location_degrees) if location_degrees is not None else 0.0
    
    @staticmethod
    def _location_cl_str(rad) -> str:
        return Individual._radians_to_fraction(rad)
    @staticmethod
    def uniform_location_between_parents(mom_ind, dad_ind, current_gen, tries=1000):
        loc1, loc2 = mom_ind.location, dad_ind.location
        # Calculate clockwise and counterclockwise distances
        clockwise_dist = (loc2 - loc1) % (2 * np.pi)
        counterclockwise_dist = (loc1 - loc2) % (2 * np.pi)

        
        # Choose the shorter path
        if clockwise_dist <= counterclockwise_dist:
            start, end = loc1, loc1 + clockwise_dist
        else:
            start, end = loc2, loc2 + counterclockwise_dist
        # Generate a location uniformly between start and end
        other_locations=[ind.location for ind in current_gen  ]
        child_location=None
        count=0
        while count<tries:
            count+=1
            child_location = np.random.uniform(start, end) % (2 * np.pi)
            if child_location not in other_locations:
                return child_location
        
        #print(f"[DEBUG] No valid location found for (loc1, loc2) = {(loc1, loc2)} and (clockwise_dist, counterclockwise_dist) = {(clockwise_dist, counterclockwise_dist)}.")
        #print(f"[DEBUG] Other loations:")
        #print(f"{other_locations}")
        raise ValueError(f"No valid location found for (loc1, loc2) = {(loc1, loc2)} and (clockwise_dist, counterclockwise_dist) = {(clockwise_dist, counterclockwise_dist)}.")
    @staticmethod
    def uniform_location_on_circle(mom_ind, dad_ind, current_gen, tries=1000):

        other_locations=[ind.location for ind in current_gen  ]
        child_location=None
        count=0
        while count<tries:
            count+=1
            child_location = np.random.uniform(0, 2 * np.pi) % (2 * np.pi)
            if child_location not in other_locations:
                return child_location
        
        #print(f"[DEBUG] Other loations:")
        #print(f"{other_locations}")
        raise ValueError(f"No valid location found.")
    
    @staticmethod
    def in_jupyter() -> bool:
        # Check if 'ipykernel' or 'jupyter_client' is in sys.modules
        return 'ipykernel' in sys.modules or 'jupyter_client' in sys.modules
    @staticmethod
    def _radians_to_fraction(rads: float) -> str:

        radians = rads % (2 * np.pi)

        multiple_of_pi = sp.Rational(radians / np.pi).limit_denominator()
        
        if multiple_of_pi == 0:
            return "0"
        elif multiple_of_pi == 1:
            return "π"
        elif multiple_of_pi == -1:
            return "-π"
        else:
            numerator = multiple_of_pi.p
            denominator = multiple_of_pi.q
            if Individual.in_jupyter():
                return f"\\frac{{{'' if numerator==1 else numerator}π}}{{{denominator}}}"  # LaTeX format
            else:
                return f"({numerator}*π)/{denominator}"  # Console format
    @staticmethod
    def _location_cl_str(rad) -> str:
        return Individual._radians_to_fraction(rad)
    @staticmethod
    def distance_between(ind1: "Individual", ind2: "Individual") -> float:
        distance = abs(ind1.location - ind2.location)
        if distance > np.pi:
            # distance is the shortest path so if the distance is clockwise we need the distance counter-clockwise
            distance = 2 * np.pi - distance
        return distance
    
    
    @staticmethod
    def midpoint_between(ind1: "Individual", ind2: "Individual") -> float:
        loc1, loc2 = ind1.location, ind2.location

        shortest_distance = Individual.distance_between(ind1, ind2)
        
        if (loc2 - loc1) % (2 * np.pi) <= np.pi:
            # Clockwise is the shorter path
            midpoint = (loc1 + shortest_distance / 2) % (2 * np.pi)
        else:
            # Counterclockwise is the shorter path
            midpoint = (loc1 - shortest_distance / 2) % (2 * np.pi)
        
        return midpoint
    


    def location_str(self) -> str:
        return self._radians_to_fraction(self.location)
    def __repr__(self) -> str:
        return (f"Individual(mom_id={self.mom_id}, dad_id={self.dad_id}, location={self.location})")
    def __str__(self) -> str:
        return f"Individual(mom_id={self.mom_id}, dad_id={self.dad_id}, location={self.location_str()})"

    def set_parents(self, pop_size: int) -> Tuple[int, int]:
        self.mom_id = random.randint(0, pop_size - 1)
        while True:
            self.dad_id = random.randint(0, pop_size - 1)
            if self.dad_id != self.mom_id:
                break
        return self.mom_id, self.dad_id
    
    #def set_location_middle_point_of_parents(self, mom_ind, dad_ind):
    #    return Individual.midpoint_between(mom_ind, dad_ind)
    def is_sibling(self, other: "Individual") -> bool:
        # Two individuals are siblings if they share at least one parent

        return (self.mom_id == other.mom_id and self.mom_id != -1) or \
           (self.dad_id == other.dad_id and self.dad_id != -1)
    
    def old_is_cousin(self, other: "Individual", previous_gen: List["Individual"]) -> bool:
    # Check if their mothers or fathers are siblings (meaning they share a grandparent)
        mom_self = self.mom_id
        dad_self = self.dad_id
        mom_other = other.mom_id
        dad_other = other.dad_id

        if previous_gen is None:
            return False

    # Check if their parents are siblings
        return ((previous_gen[mom_self].is_sibling(previous_gen[mom_other]) or previous_gen[dad_self].is_sibling(previous_gen[dad_other])) or (previous_gen[mom_self].is_sibling(previous_gen[dad_other]) or previous_gen[dad_self].is_sibling(previous_gen[mom_other]))) 

    def is_cousin(self, other: "Individual", previous_gen: List["Individual"]) -> bool:
        # Get grandparents by accessing parents' mom_id and dad_id
        if self.mom_id == -1 or self.dad_id == -1 or other.mom_id == -1 or other.dad_id == -1:
            return False

        my_grandparents = {
            previous_gen[self.mom_id].mom_id, previous_gen[self.mom_id].dad_id,
            previous_gen[self.dad_id].mom_id, previous_gen[self.dad_id].dad_id
        }

        other_grandparents = {
            previous_gen[other.mom_id].mom_id, previous_gen[other.mom_id].dad_id,
            previous_gen[other.dad_id].mom_id, previous_gen[other.dad_id].dad_id
        }

        # Count shared grandparents, exclude -1 which represents "no ancestor"
        shared_grandparents = my_grandparents.intersection(other_grandparents) - {-1}

        return len(shared_grandparents) > 0

    def find_single_male_partner(self, mu: float, pop_size: int, previous_gen: List["Individual"], preprevious_gen: List["Individual"], kappa_parameter: int,avoid_relatives: str=None) -> int:
        # Step 1: Filter only unpaired fathers (and exclude siblings and cousins)
        if avoid_relatives is None:
            single_males = [
                (i, previous_gen[i].location)
                for i in range(math.floor(pop_size / 2), pop_size)
                if previous_gen[i].pair_id == -1
            ]
        elif avoid_relatives == '':
            single_males = [
                (i, previous_gen[i].location)
                for i in range(math.floor(pop_size / 2), pop_size)
                if previous_gen[i].pair_id == -1
            ]    
        elif avoid_relatives == 'siblings':
            single_males = [
                (i, previous_gen[i].location)
                for i in range(math.floor(pop_size / 2), pop_size)
                if previous_gen[i].pair_id == -1
                and not self.is_sibling(previous_gen[i])  # Exclude siblings
                
            ]

        elif avoid_relatives=='cousins':
            single_males = [
                (i, previous_gen[i].location)
                for i in range(math.floor(pop_size / 2), pop_size)
                if previous_gen[i].pair_id == -1
                and not self.is_cousin(previous_gen[i], preprevious_gen)  # Exclude cousins
            ]
        elif avoid_relatives=='siblingscousins':
            single_males = [
                (i, previous_gen[i].location)
                for i in range(math.floor(pop_size / 2), pop_size)
                if previous_gen[i].pair_id == -1
                and not self.is_sibling(previous_gen[i])  # Exclude siblings
                and not self.is_cousin(previous_gen[i], preprevious_gen)  # Exclude cousins
            ]
        else:
            raise ValueError(f'Invalid avoid_relatives value : {avoid_relatives}')

        # Step 2: Handle the case of no available mates gracefully
        if not single_males:
            #print("[DEBUG] No valid male candidates found.")
            return -9  # Return a signal for no available partner

        # Step 3: Extract indices and locations of unpaired fathers
        male_indices, male_locations = zip(*single_males)

        # Step 4: Calculate distances from the mother to each unpaired father
        distances = [Individual.distance_between(self, previous_gen[i]) for i in male_indices]

        # Step 5: Calculate weights using von Mises distribution
        center_mu = 0 if mu is None else mu
        weights = [vonmises.pdf(dist, kappa=kappa_parameter, loc=center_mu) for dist in distances]

        # Step 6: Normalize weights to make them probabilities
        normalized_weights = [w / sum(weights) for w in weights]

        # Step 7: Select a father based on the normalized weights
        selected_index = np.random.choice(male_indices, p=normalized_weights)
        selected_father = previous_gen[selected_index]

        return selected_index
    def old2_find_single_male_partner(self, mu: float, pop_size: int, previous_gen: List[int], kappa_parameter:int) -> int:
    # Step 1: Filter only unpaired fathers (and exclude siblings)
        single_males = [
            (i, previous_gen[i].location)
            for i in range(math.floor(pop_size / 2), pop_size)
            if previous_gen[i].pair_id == -1 and not self.is_sibling(previous_gen[i])
    ]

    # Step 2: Handle the case of no available mates gracefully
        if not single_males:
            #print("[DEBUG] No valid male candidates found.")
            return -9  # Return a signal for no available partner

    # Step 3: Extract indices and locations of unpaired fathers
        male_indices, male_locations = zip(*single_males)
    
    # Step 4: Calculate distances from the mother to each unpaired father
        distances = [Individual.distance_between(self, previous_gen[i]) for i in male_indices]
    
    # Step 5: Calculate weights using von Mises distribution
        center_mu = 0 if mu is None else mu
        weights = [vonmises.pdf(dist, kappa=kappa_parameter, loc=center_mu) for dist in distances]

    # Step 6: Normalize weights to make them probabilities
        normalized_weights = [w / sum(weights) for w in weights]

    # Step 7: Select a father based on the normalized weights
        selected_index = np.random.choice(male_indices, p=normalized_weights)
        selected_father = previous_gen[selected_index]

        #print(f"[DEBUG] Mother ID: {self.mom_id}, Selected Father: {selected_index}")
        return selected_index

    def set_parents_straight_monoamorous_couple(self, pop_size: int, previous_gen: List["Individual"],preprevious_gen:List["Individual"], kappa_parameter: int, avoid_relatives: str=None) -> Tuple[int, int]:
        

        # Randomly select a mother
        remaining_moms=list(range(0,math.floor(pop_size / 2)))
        #print(f'Remaing moms:remaining_moms={remaining_moms} ')
    
        while len(remaining_moms)>0 :
            mom_id = random.choice(remaining_moms)
            start_time = time.perf_counter()
            remaining_moms.remove(mom_id)

            

             # Proceed only if the mother isn't already paired
            if previous_gen[mom_id].pair_id == -1:
                dad_id = previous_gen[mom_id].find_single_male_partner(None, pop_size, previous_gen, preprevious_gen, kappa_parameter,avoid_relatives)
                
                # Retry if no valid male was found (-9)
                if dad_id == -9:
                    #print(f"[DEBUG] No mate found for Mother ID {mom_id}.")
                    continue
                
                passed_time = time.perf_counter() - start_time
                passed_seconds=round(passed_time, 2)


                # Pair the two if successful
                self.dad_id = dad_id
                self.mom_id = mom_id
                previous_gen[self.mom_id].pair_id = dad_id
                previous_gen[self.dad_id].pair_id = mom_id
                #print(f'[DEBUG] Mom_id : {mom_id} selected father : {dad_id}')
                return self.mom_id, self.dad_id, passed_seconds
            else:
                self.mom_id = mom_id
                self.dad_id = previous_gen[mom_id].pair_id
                #print(f'[DEBUG] Mom_id : {mom_id} has already father : {previous_gen[mom_id].pair_id}')
                return self.mom_id, self.dad_id, -1
        if len(remaining_moms) == 0 :
            raise ValueError(f'At this point,the function hasnt return so this means that this individual couldnt select a legit mom either because all moms couldnt find a free partner or either because all moms are sisters or cousins with all men.(a very extreme case)')

            
    def old2_set_parents_straight_monoamorous_couple(self, pop_size: int, previous_gen: List[int],kappa_parameter:int) -> Tuple[int, int]:
        while True:
        # Randomly select a mother
            self.mom_id = random.randint(0, math.floor(pop_size / 2) - 1)

        # Proceed only if the mother isn't already paired
            if previous_gen[self.mom_id].pair_id == -1:
                dad_id = previous_gen[self.mom_id].find_single_male_partner(None, pop_size, previous_gen,kappa_parameter)
            
            # Retry if no valid male was found (-9)
                if dad_id == -9:
                    #print("[DEBUG] No male found for this mother. Retrying...")
                    continue  # Loop to reselect a new mother
            
            # Pair the two if successful
                self.dad_id = dad_id
                previous_gen[self.mom_id].pair_id = dad_id
                previous_gen[self.dad_id].pair_id = self.mom_id
                break  # Exit the loop on successful pairing
            else:
                self.dad_id = previous_gen[self.mom_id].pair_id
                break  # Pair already exists

        return self.mom_id, self.dad_id

    
    def old_find_single_male_partner(self, mu: float, pop_size: int, previous_gen: List[int],kappa_parameter: int) -> int:
    # Step 1: Filter only unpaired fathers
        single_males = [
        (i, previous_gen[i].location) 
        for i in range(math.floor(pop_size / 2), pop_size)
        if previous_gen[i].pair_id == -1  # Only unpaired fathers
        ]

    # Step 2: Handle edge case where no unpaired fathers are available
        if not single_males:
            raise ValueError("No single male partners available!")

    # Step 3: Extract indices and locations of unpaired fathers
        male_indices, male_locations = zip(*single_males)

    # Step 4: Calculate distances from the mother to each unpaired father
        distances = [Individual.distance_between(self, previous_gen[i]) for i in male_indices]

    # Step 5: Calculate weights using von Mises distribution centered on the mother's location (`mu`)
        center_mu = 0 if mu is None else mu
        weights = [vonmises.pdf(dist, kappa=kappa_parameter, loc=center_mu) for dist in distances]

    # Step 6: Normalize weights to make them probabilities
        normalized_weights = [w / sum(weights) for w in weights]

    # Step 7: Select a father based on the normalized weights
        selected_index = np.random.choice(male_indices, p=normalized_weights)
        return selected_index
    
    def old_set_parents_straight_monoamorous_couple(self, pop_size: int, previous_gen: List[int],kappa_parameter: int) -> Tuple[int,int]:
        self.mom_id = random.randint(0, math.floor(pop_size/2)-1)
        
        # uniform random
        # self.dad_id = random.randint(round(pop_size/2), pop_size-1)
        
        # Normal distribution based
        if previous_gen[self.mom_id].pair_id == -1:
            self.dad_id = previous_gen[self.mom_id].find_single_male_partner(None, pop_size, previous_gen,kappa_parameter)
            previous_gen[self.mom_id].pair_id = self.dad_id
            previous_gen[self.dad_id].pair_id = self.mom_id
        else:
            self.dad_id = previous_gen[self.mom_id].pair_id

        return self.mom_id, self.dad_id
    
    def set_break_pos(self, chrom_break_pos: List[int], recomb_rate: float) -> None:
        seq_len = chrom_break_pos[-1]

        for i in range(2):
            break_pos = []
            for j in range(len(chrom_break_pos) - 1):
                if random.randint(0, 1) == 1:
                    break_pos.append(chrom_break_pos[j])
            
            num_recomb = np.random.default_rng().binomial(seq_len, recomb_rate)
            for j in range(num_recomb):
                break_pos.append(random.randint(1, seq_len - 2))

            break_pos.sort()
            break_pos.append(seq_len)

            if i == 0:
                self.mom_break_pos = break_pos
            else:
                self.dad_break_pos = break_pos
    
