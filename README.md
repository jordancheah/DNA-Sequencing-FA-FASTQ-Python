# DNA Sequencing: Sampe Analysis with lambda_virus.fa and FASTQ file

### Genome Data

* lambda virus genome, at https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/lambda_virus.fa

* real DNA sequencing reads derived from a human https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/ERR037900_1.first1000.fastq

SOURCE: Ajay, S. S., Parker, S. C., Abaan, H. O., Fajardo, K. V. F., & Margulies, E. H. (2011). Accurate and comprehensive sequencing of personal genomes. Genome research, 21(9), 1498-1505. 

### Algorithms Implmented

* Naive exact matching algorithm that is strand-aware. That is, instead of looking only for occurrences of P in T, additionally look for occurrences of the reverse complement of P in T

* Naive matching algirithm that allows up to 2 mismatches per occurrence.
