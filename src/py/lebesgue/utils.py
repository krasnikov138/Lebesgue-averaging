import sys
import math
from time import perf_counter

class ByteSize:
    unit_mapping = {
        'b': 0, 'Kb': 1, 'Mb': 2, 'Gb': 3
    }

    def __init__(self, size):
        self._bytes = size
        self._unit = self.get_unit(size)
        self._normalized = size / self.unit_multiplier(self._unit)

    @staticmethod
    def get_unit(size):
        degree = min(3, int(math.log2(size) / 10))
        return ['b', 'Kb', 'Mb', 'Gb'][degree]
        
    @classmethod
    def unit_multiplier(cls, unit):
        return 2 ** (10 * cls.unit_mapping[unit])

    @property
    def normalized(self):
        return self._normalized

    @property
    def unit(self):
        return self._unit
    
    def __str__(self):
        return f'{self._normalized:.2f} {self._unit}'


class ProgressBar:
    def __init__(self, total_length, msg=None):
        self.current_length = 0
        self.total_length = total_length
        self.msg = msg
        self._max_text_length = 0
        self._start = perf_counter()

    def update(self, chunk_size):
        self.current_length += chunk_size
        progress = self.current_length / self.total_length
        ready = int(50 * progress)
        text = '[' + '=' * ready + '>' + ' ' * (50 - ready) + ']: '
        if self.msg:
            text += self.msg + ' - ' 
        text += f'{ByteSize(self.current_length)} / {ByteSize(self.total_length)}'

        self.print(text)        

    def print(self, text):
        self._max_text_length = max(self._max_text_length, len(text))
        text += ' ' * (self._max_text_length - len(text)) + '\r'
        sys.stdout.write(text)
        sys.stdout.flush()

    def complete(self, name):
        time = perf_counter() - self._start
        text = f'{name}: {ByteSize(self.total_length)} was download in {time:.1f} s.'
        self.print(text)
        print()
