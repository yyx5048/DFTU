import re
import pandas as pd
import numpy as np
from collections import defaultdict

df = pd.read_csv("Final_DMREF_Materials_List_Duplicate.csv")


def list_duplicates(seq):
    #--corresponding to the csv list idex with +2
    tally = defaultdict(list)
    for i, item in enumerate(seq):
        tally[item].append(i)
    return ((key,(np.array(locs)+2).tolist()) for key,locs in tally.items()
                            if len(locs)>1)


formula_list = []
for idx,val in enumerate(df["Formula"].values):

    #--split
    tmp_form = re.findall('[A-Z][^A-Z]*', re.sub(" ", "" ,val))
    tmp_form.sort()

    formula_list.append("".join(tmp_form))

dup_idx = 0
for dup in list_duplicates(formula_list):
    dup_idx += 1
    print(dup)
print(dup_idx)
