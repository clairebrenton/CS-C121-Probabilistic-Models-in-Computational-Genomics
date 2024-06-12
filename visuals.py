import pandas as pd
import matplotlib.pyplot as plt


# Creating the DataFrame manually since the file loading failed
data = pd.DataFrame({
    'chromosome': ['chr1:1000000-2000000'] * 7,
    'position': [172741, 325026, 375797, 423797, 518726, 568632, 868896],
    'putative_genotype': ['AA', 'GG', 'AA', 'AA', 'GG', 'TT', 'TT'],
    'posterior_probability': [0.9947035619462352, 0.23212160875773435, 0.9956845526575606, 0.44434615725878507, 0.013658578683693861, 0.13421011701263563, 0.7885189537040443],
    'n_reads': [5, 3, 8, 7, 7, 6, 5]
})

# Posterior Probability vs. Genomic Position
plt.figure(figsize=(10, 6))
plt.scatter(data['position'], data['posterior_probability'], c='blue', marker='o')
plt.title('Posterior Probability vs. Genomic Position')
plt.xlabel('Genomic Position')
plt.ylabel('Posterior Probability')
plt.grid(True)
plt.savefig('posterior_probability_vs_position.png')
plt.show()

# Number of Reads vs. Genomic Position
plt.figure(figsize=(10, 6))
plt.bar(data['position'].astype(str), data['n_reads'], color='green')  # Convert positions to string for better labeling
plt.title('Number of Reads vs. Genomic Position')
plt.xlabel('Genomic Position')
plt.ylabel('Number of Reads')
plt.grid(True)
plt.savefig('n_reads_vs_position.png')
plt.show()