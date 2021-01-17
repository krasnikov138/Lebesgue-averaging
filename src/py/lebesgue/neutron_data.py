import re
import os
import posixpath
import requests
from zipfile import ZipFile
from collections import defaultdict

import pandas as pd

from .utils import ProgressBar


def parse_index(index):
    index = index.replace('\r', '')
    start = index.find('----\n') + 5
    start = index.find('----\n', start) + 5
    end = index.find('----', start)
    table = index[start:end]
    df = defaultdict(list)
    for row in table.split('\n'):
        if not row:
            continue
        df['mat'].append(int(row[5:11]))
        df['material'].append(row[11:24].strip())
        df['laboratory'].append(row[24:37].strip())
        df['date'].append(row[37:50].strip())
        df['authors'].append(row[50:85].strip())
        href = re.search(r'<a href="([^"]+)">', row).group(1)
        df['href'].append(href)
    df = pd.DataFrame(df)
    material = df.material.str.extract(r'(\d+)-([a-zA-Z]+)-(\d+)(M?)', expand=True)
    material.columns = ['Z', 'atom', 'A', 'modified']
    material = material.astype({'Z': 'int', 'A': 'int'})
    material.modified = (material.modified == 'M')
    df = pd.concat([df, material], axis=1).drop(columns='material')
    return df[['mat', 'atom', 'Z', 'A', 'modified', 'laboratory', 'date', 'authors', 'href']]


def download_file(url, data_path):
    fname = os.path.join(data_path, posixpath.basename(url))
    with requests.get(url, stream=True) as resp:
        total = int(resp.headers['Content-Length'])
        bar = ProgressBar(total, os.path.basename(fname))
        with open(fname, 'wb') as f:
            for chunk in resp.iter_content(chunk_size=8192):
                bar.update(len(chunk))
                f.write(chunk)
        bar.complete(os.path.basename(fname))
    return fname


def extract_file(fname, data_path):
    name = os.path.basename(fname)
    with ZipFile(fname) as zfile:
        names = zfile.namelist()
        if len(names) > 1:
            print(f'Warning: archive {name} contains more than 1 file.')
        else:
            name = names[0]
        zfile.extractall(data_path)
    print(f'{name} was extracted to {data_path}')


class EndfIndex:
    url = 'https://www-nds.iaea.org/public/download-endf/'

    def __init__(self, index, catalog):
        self.index = index
        self.catalog = catalog

    def query(self, mat=None, atom=None, Z=None, A=None, modified=None):
        mask = pd.Series(True, index=self.index.index)
        for name, value in zip(['mat', 'atom', 'Z', 'A', 'modified'], [mat, atom, Z, A, modified]):
            if value is None:
                continue
            if self.is_iterable(value):
                mask &= self.index[name].isin(value)
            else:
                mask &= (self.index[name] == value)
        return EndfIndex(self.index.loc[mask], self.catalog)

    @staticmethod
    def is_iterable(value):
        try:
            it = iter(value)
            return True
        except TypeError:
            return False

    def __repr__(self):
        if self.index.empty:
            return 'Empty endf-index'
        return repr(self.index.loc[:, self.index.columns.drop('href')])

    def show_all(self):
        if self.index.empty:
            print('Empty endf-index')
        else:
            print(self.index.loc[:, self.index.columns.drop('href')].to_string())

    def download(self, data_path):
        os.makedirs(data_path, exist_ok=True)
        for href in self.index.href:
            url = posixpath.join(self.url, self.catalog, href)
            fname = download_file(url, data_path)
            extract_file(fname, data_path)
            os.remove(fname)


class EndfDatabase(EndfIndex):
    versions = {
        7: 'ENDF-B-VII.1',
        8: 'ENDF-B-VIII.0'
    }

    def __init__(self, version=8):
        assert version in self.versions.keys(), f"Supported versions are: {list(self.versions.keys())}"
        self.catalog = self.versions[version]
        index_url = posixpath.join(self.url, self.catalog, 'n-index.htm')
        self.index = parse_index(requests.get(index_url).text)
