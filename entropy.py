import sys
from collections import Counter
import numpy as np
from Bio import AlignIO

entropy = []
print("site,protein,entropy,mutations")
for prot in sys.argv[1:]:
    alignment = AlignIO.read(open(prot), "fasta")
    prot_name = prot.split(".")[0]
    for pos in range(alignment.get_alignment_length()):
        state_counts = Counter(str(alignment[:, pos]))
        ref = alignment[0, pos]
        state_counts.pop('X', None)
        state_counts.pop('-', None)
        nsb_array = np.array(list(state_counts.values()))
        fg_entropy = -np.sum((nsb_array[nsb_array != 0] / nsb_array[nsb_array != 0].sum()) * np.log2(nsb_array[nsb_array != 0] / nsb_array[nsb_array != 0].sum()))
        mutation = ""
        if fg_entropy == 0.0:
            fg_entropy = 0
            mutation = "NA;"
        else:
            for key in state_counts:
                if not key == ref:
                    mutation += f"{ref}->{key}({state_counts[key]});"
        print(f"{pos+1},{prot_name},{fg_entropy:.8f},{mutation[:-1]}")
