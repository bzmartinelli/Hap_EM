# Hap_EM

Requirements: 

1. Python version 2.7.10 or above.
2. Unix based environment.


Input file
A tab separated file containing the coordinate of the CpG, their methylation state and the sequence read ID. Our code accepts 3 different notations representing the methylation status of the CpG:

   methylated   unmethylated
1)	m	u
2)	Z	z
3)	1	0


Example of input file: 

coordinate	methylation status	read sequence ID
166205		m		2a30292e95256700 
166205		m		bc082d70a7caec58
166205		u		d6a4eed8da88a9f9
166244		u		2010f4587de3ae8d
166244		m		caee1fa8acee2ad5
