from typing import List

class Sequence:
    def __init__(self) -> None:
        self.ids: List[int] = []
        self.lengths: List[int] = []

    def add_segment(self, id: int, length: int) -> None:
        if length > 0:
            self.ids.append(id)
            self.lengths.append(length)

    def extend(self, other: 'Sequence') -> None:
        self.ids.extend(other.ids)
        self.lengths.extend(other.lengths)

    def get_sequence(self, start_pos: int, end_pos: int) -> 'Sequence':
        i = 0
        pointer = 0
        out_seq = Sequence()
        while True:
            pointer += self.lengths[i]
            if start_pos >= pointer:
                i += 1
            elif end_pos > pointer:
                out_seq.add_segment(self.ids[i], pointer - start_pos)
                i += 1
                start_pos = pointer
            else:
                out_seq.add_segment(self.ids[i], end_pos - start_pos)
                return out_seq
