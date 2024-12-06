import pandas as pd
import pysam
import logomaker
import matplotlib.pyplot as plt
import h5py
import hdf5plugin
import csv
import numpy as np
from matplotlib.lines import Line2D



# data = pd.read_csv('hits_unique.tsv', sep='\t')
# pd.set_option('display.max_columns', None)

# def fetch_sequence(row):
#     return fasta.fetch(row['chr'], row['start_untrimmed'], row['end_untrimmed'])

# fasta = pysam.FastaFile('hg38.fa')
# data['sequence'] = data.apply(fetch_sequence, axis=1)
# data = data.sort_values(by='hit_importance', ascending=False)
# print(data.head())

# sequence = data['sequence'][78336]
# motif_data = pd.read_csv('motif_data.tsv', sep='\t')
# motif_scale = motif_data[motif_data['motif_name'] == 'neg_patterns.pattern_0']['motif_scale'].values[0]
# motif_start = motif_data[motif_data['motif_name'] == 'neg_patterns.pattern_0']['motif_start'].values[0]
# motif_end = motif_data[motif_data['motif_name'] == 'neg_patterns.pattern_0']['motif_end'].values[0]
# motif_df = pd.DataFrame(0, columns=['A', 'C', 'G', 'T'], index=range(len(sequence)))

# for i, nucleotide in enumerate(sequence):
#     motif_df.at[i, nucleotide] = 1 if (i < motif_start or i > motif_end) else motif_scale

# Open the H5 file (read mode)
file_path = 'variant_shap.counts.h5'
variant_list = pd.read_csv('variant_list.tsv', sep="\t", header=None)
with h5py.File(file_path, 'r') as h5_file:
    # Print all top-level groups and datasets
    def print_structure(name, obj):
        print(f"{name}: {obj}")

    h5_file.visititems(print_structure)

    shap_seq = h5_file['projected_shap/seq']
    print("Shape of projected_shap/seq:", shap_seq.shape)
    print("Data type:", shap_seq.dtype)

    print(shap_seq[0])   
    print(shap_seq[0].shape)  
    variant = input("Variant id (chr[n]_[num]_[a1]_[a2]_b38): ")
    variant_num = variant_list[variant_list[4] == variant].index[0]
    print(f'Variant row number: {variant_num}')
    
    data1 = shap_seq[variant_num][:,957:1157]
    data2 = shap_seq[variant_num + 13279][:,957:1157]



# Data for the first plot
contribution_df1 = pd.DataFrame(
    data1,
    index=['A', 'C', 'G', 'T'],
    columns=[i for i in range(data1.shape[1])]
)

logo1 = logomaker.Logo(contribution_df1.T)
logo1.style_spines(visible=False)
logo1.style_spines(spines=['left', 'bottom'], visible=True)
logo1.ax.set_ylabel("Contribution")
logo1.ax.set_xlabel("Position")
logo1.ax.set_title(f"Predicted Contribution with Reference Allele for Variant {variant}")


print(contribution_df1)

#Show the first plot
plt.figure(1)
plt.xticks(np.linspace(0, 200, 5), labels=np.linspace(-100, 100, 5))

logo1.ax.axvline(x=100, color='red', linestyle='--', linewidth=2)
line = Line2D([100, 100], [0, 1], color='red', linestyle='-', linewidth=40, alpha=0.2)  
logo1.ax.add_line(line)
plt.ylim(min(data1.min(), data2.min()) - 0.01, max(data1.max(), data2.max()) + 0.01)

plt.show(block=False)

# Data for the second plot
contribution_df2 = pd.DataFrame(
    data2,
    index=['A', 'C', 'G', 'T'],
    columns=[i for i in range(data2.shape[1])]
)
logo2 = logomaker.Logo(contribution_df2.T)
logo2.style_spines(visible=False)
logo2.style_spines(spines=['left', 'bottom'], visible=True)
logo2.ax.set_ylabel("Contribution")
logo2.ax.set_xlabel("Position")
logo2.ax.set_title(f"Predicted Contribution with Alternate Allele for Variant {variant}")

print(contribution_df2)
# Show the second plot
plt.figure(2)
logo2.ax.axvline(x=100, color='red', linestyle='--', linewidth=2)
line2 = Line2D([100, 100], [0, 1], color='red', linestyle='-', linewidth=40, alpha=0.2)  
logo2.ax.add_line(line2)
plt.xticks(np.linspace(0, 200, 5), labels=np.linspace(-100, 100, 5))
plt.ylim(min(data1.min(), data2.min()) - 0.01, max(data1.max(), data2.max()) + 0.01)
plt.show()

