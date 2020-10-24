import pandas as pd
from lebesgue.bindings import PyEndfFile

class EndfFile(PyEndfFile):

    @property
    def table_of_content(self):
        return pd.DataFrame(super().table_of_content, columns=['mf', 'mt'])

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()
