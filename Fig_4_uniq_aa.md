Make a output directory

```
mkdir -p ./uniq_aa_results/merged_fasta
```

Merge FASTA files except for a query group (NRC2, NRC3, NRC4a, NRC4b, or NRC4other)

```
cat $(ls *.afa | grep -v 'NRC2.afa') > ./uniq_aa_results/merged_fasta/excl_NRC2.afa
cat $(ls *.afa | grep -v 'NRC3.afa') > ./uniq_aa_results/merged_fasta/excl_NRC3.afa
cat $(ls *.afa | grep -v 'NRC4a.afa') > ./uniq_aa_results/merged_fasta/excl_NRC4a.afa
cat $(ls *.afa | grep -v 'NRC4b.afa') > ./uniq_aa_results/merged_fasta/excl_NRC4b.afa
cat $(ls *.afa | grep -v 'NRC4other.afa') > ./uniq_aa_results/merged_fasta/excl_NRC4other.afa
```

We prepared a Python script ```uniq_aa.py``` to determine unique amino acids in a query group based on frequency and ratio thresholds. The script requires the Biopython library to run.


```
Usage: uniq_aa.py <QUERY_FASTA> <SUBJECT_FASTA> <MIN_FREQUENCY> <MIN_RATIO>
```

**Arguments:**

**QUERY_FASTA:** Path to the FASTA file containing the query group sequences.

**SUBJECT_FASTA:** Path to the FASTA file containing the subject sequences (excluding the query group).

**MIN_FREQUENCY:** Minimum frequency threshold for an amino acid to be considered unique in the query group.

**MIN_RATIO:** Minimum ratio threshold for an amino acid to be considered unique in the query group. The ratio is calculated as the frequency of the amino acid in the query group divided by its frequency in the subject group.


```
./uniq_aa.py NRC2.afa ./uniq_aa_results/merged_fasta/excl_NRC2.afa 0.8 50 > ./uniq_aa_results/NRC2_to_others_f0.8_r50.txt
./uniq_aa.py NRC3.afa ./uniq_aa_results/merged_fasta/excl_NRC3.afa 0.8 50 > ./uniq_aa_results/NRC3_to_others_f0.8_r50.txt
./uniq_aa.py NRC4a.afa ./uniq_aa_results/merged_fasta/excl_NRC4a.afa 0.8 50 > ./uniq_aa_results/NRC4a_to_others_f0.8_r50.txt
./uniq_aa.py NRC4b.afa ./uniq_aa_results/merged_fasta/excl_NRC4b.afa 0.8 50 > ./uniq_aa_results/NRC4b_to_others_f0.8_r50.txt
./uniq_aa.py NRC4other.afa ./uniq_aa_results/merged_fasta/excl_NRC4other.afa 0.8 50 > ./uniq_aa_results/NRC4other_to_others_f0.8_r50.txt
```

Merge all the results

```
paste ./uniq_aa_results/*_f0.8_r50.txt > ./uniq_aa_results/all_uniq_aa_results.txt
```

Results in all_uniq_aa_results.txt are used to make Excel sheet.