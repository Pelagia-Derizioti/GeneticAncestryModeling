from typing import List, Tuple
from sequence import Sequence

class ChromPair:
    def __init__(self):
        # chrom_pair is expected to be a list with two Sequence objects
        self.chrom_pair = [Sequence(), Sequence()]

    def get_chrom(self, break_pos: List[int]) -> Sequence:
        chrom_index = 0
        # Access the first sequence's chromosome data
        seq = self.chrom_pair[chrom_index].get_sequence(0, break_pos[0])
        chrom_index = 1 - chrom_index

        for i in range(len(break_pos) - 1):
            if break_pos[i] != break_pos[i + 1]:
                seq.extend(self.chrom_pair[chrom_index].get_sequence(break_pos[i], break_pos[i + 1]))
                chrom_index = 1 - chrom_index

        return seq
    
    def get_ibd_length_whole_segments(self, other: "ChromPair") -> int:
        """
        Computes the total length of Identical By Descent (IBD) segments shared 
        between two ChromPair objects (two individuals).
        """
        ibd_length = 0

        # Compare both chromosome copies
        for chrom_index in range(2):  # 0 = paternal, 1 = maternal
            my_seq = self.chrom_pair[chrom_index]
            other_seq = other.chrom_pair[chrom_index]

            # Compare segment-by-segment
            min_length = min(len(my_seq.ids), len(other_seq.ids))
            for i in range(min_length):
                if my_seq.ids[i] == other_seq.ids[i]:  # Check if segments are identical
                    ibd_length += my_seq.lengths[i]  # Sum the length of shared segments

        return ibd_length
    
    def get_ibd_length(self, other: "ChromPair") -> int:
        ibd_length = 0

        # Compare the two sequences (chromosome pairs)
        for my_seq, other_seq in zip(self.chrom_pair, other.chrom_pair):
            i, j = 0, 0  # Pointers for my_seq and other_seq segments
            my_start, other_start = 0, 0  # Track start positions in each sequence

            while i < len(my_seq.ids) and j < len(other_seq.ids):
                # Check if the IDs match (same ancestor)
                if my_seq.ids[i] == other_seq.ids[j]:
                    # Compute the overlap between the two segments
                    my_end = my_start + my_seq.lengths[i]
                    other_end = other_start + other_seq.lengths[j]

                    overlap_start = max(my_start, other_start)
                    overlap_end = min(my_end, other_end)

                    # Add the overlapping length to the IBD total (if there is overlap)
                    if overlap_start < overlap_end:
                        ibd_length += overlap_end - overlap_start

                # Move to the next segment in the sequence with the smaller end position
                if my_start + my_seq.lengths[i] <= other_start + other_seq.lengths[j]:
                    my_start += my_seq.lengths[i]
                    i += 1
                else:
                    other_start += other_seq.lengths[j]
                    j += 1

        return ibd_length


    def get_roh(self, seq_len: int) -> Tuple[int, int]:
        # Merge break positions from both chromosome sequences
        break_pos = [0]
        break_pos.append(self.chrom_pair[0].lengths[0])

        for k in range(1, len(self.chrom_pair[0].lengths)):
            break_pos.append(break_pos[k - 1] + self.chrom_pair[0].lengths[k])

        break_pos.append(self.chrom_pair[1].lengths[0])
        for k in range(1, len(self.chrom_pair[1].lengths)):
            break_pos.append(break_pos[k - 1] + self.chrom_pair[1].lengths[k])

        break_pos.sort()
        break_pos.append(seq_len)

        # Find ROH for each block
        roh_length = 0
        roh_count = 0
        for i in range(len(break_pos) - 1):
            if break_pos[i] != break_pos[i + 1]:
                # Compare sequence IDs to determine if ROH exists
                if self.chrom_pair[0].get_sequence(break_pos[i], break_pos[i + 1]).ids[0] == self.chrom_pair[1].get_sequence(break_pos[i], break_pos[i + 1]).ids[0]:
                    roh_length += (break_pos[i + 1] - break_pos[i])
                    roh_count += 1
        return roh_length, roh_count