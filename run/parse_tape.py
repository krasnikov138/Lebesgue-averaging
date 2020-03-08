import os
import sys
import pandas as pd

def load(fname):
    content = pd.Series(open(fname).read().split('\n'))
    content = content[content.str.slice(71, 75) == '3  1'].iloc[3:]
    content = content.str.slice(0, 66)
    content = pd.Series(''.join(content).split())
    extracted = content.str.extract(r'(\d\.\d+)([\+-])(\d)')
    na = extracted[0].isna()
    res = pd.Series(index=content.index)
    res.loc[na] = content.loc[na]
    res = res.astype('float')
    extracted.loc[:, [0, 2]] = extracted.loc[:, [0, 2]].astype('float')
    extracted.loc[extracted[1] == '-', 2] *= -1
    extracted = extracted[0] * 10**extracted[2]
    res.loc[~na] = extracted.loc[~na]
    df = pd.DataFrame({'energy': res.iloc[::2].reset_index(drop=True), 
                       'cross_section': res.iloc[1::2].reset_index(drop=True)})
    return df


if __name__ == '__main__':
	if len(sys.argv) not in [2, 3]:
		print("Usage: python parse_tape.py TAPE_PATH [OUT_NAME]")
		exit(0)
	fname = sys.argv[1]
	if len(sys.argv) == 3:
		out_name = sys.argv[2]
	else:
		out_name = os.path.dirname(fname) + '/out'

	df = load(fname)
	df.to_csv(out_name, index=False, header=False, sep=' ')
