-h, --help  show this help message and exit

# parameter
  1. input --string, Path of a input FASTA
  2. output --string, Enter the name of output csv, the output file will be created in 'ORF_output' folder behind this python.
  3. min_orf_length --integer, Minmum length(a.a.) of ORF (default = 20)
  4. seq_repeat --integer(1~4), How many times sequence repeat to find ORF (default = 4 times)

# notice
  1. Only 'ATCG' are avalible in input FASTA file sequence.
  2. The sequence will not be analized that have character except 'ATCG', and the id of sequence will be show in terminal with unavalible characters.
  3. Output file will be created in a folder named 'ORF_output' next to this python file, and the file will be named after string you entered.