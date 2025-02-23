from typing import List

def get_chrom_break_pos(chrom_lengths: List[int]) -> List[int]:
    chrom_break_pos = [0 for i in range(len(chrom_lengths) + 1)]
    for i in range(len(chrom_lengths)):
        chrom_break_pos[i + 1] = chrom_break_pos[i] + chrom_lengths[i]
    return chrom_break_pos
