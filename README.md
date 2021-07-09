# GeneEssentiality
Extract and calculate features of nucleotide sequences and their longest translated CDS. Then use these features to train a machine learning model

To install the package use the command:
```
devtools::install_github("g1o/GeneEssentiality")
```
To extract the features from a nucleotide sequence, write the path to the fasta file to the Extract function. Gziped fastas can also be read:
```
Extract(FASTA_PATH="/home/USER/seq.fasta.gz",CPU=5) 
```
*Warning: Do not use too many cores, the RAM increases linearly with the number of used threads. It is also recomended to extract features in a new R session, with only the package loaded to avoid RAM issues if using over 10 threads.  
Using 10 threads consumes about 12GB in a clean session. Tested with the procariote data.*

Using reproduce_bacteria_results(CPU=10) will take about 4 hours to end. This does the complete process, from the extraction of features, training and testing to plotting the results from the procariote data. 

Using reproduce_insects_results(CPU=10) will reproduce the insects results using the features already extracted and present in this package. 

reproduce_insects_results(set1=First,set2=Second,CPU=10) will do the same comparisons as it does with the insects data, training models with the fist set and testing with the second, vice-versa. 





