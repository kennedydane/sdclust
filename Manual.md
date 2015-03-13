# sd\_clust version 0.1.0 pre-alpha #

Usage: sd\_clust `[`options`]` 

&lt;filename&gt;



> `[`options`]`
> > -v, --version        : show version number (stdout)
> > -h, --help           : show this help screen (stdout)
> > -p, --progress       : show percentage progress (stderr)
> > -f, --fasta          : show summary about the fasta sequences (stderr)


> -c, --show\_compact   : show the compact cluster table (stdout)
> -C, --show\_extended  : show the extended cluster table (stdout)
> -q, --single         : do single linkage clustering (default)
> -Q, --complete       : find all overlapping ESTs
> -n, --no\_rc          : don't evaluate reverse compliment
> -H, --heur\_only      : cluster based on heuristics only

> -S, --skip 

&lt;int&gt;

     : set number of words to skip (default = 10)
> -s, --sd 

&lt;float&gt;

     : set standard deviation threshold (default = 4)
> -t, --wc\_thresh 

&lt;int&gt;

: set word count threshold (default = 3)
> -w, --word\_len 

&lt;int&gt;

 : set word length (default = 14)

> -L, --d2\_l 

&lt;int&gt;

     : set d2 window length (default = 100)
> -T, --d2\_t 

&lt;int&gt;

     : set d2 threshold (default = 40)
> -W, --d2\_w 

&lt;int&gt;

     : set d2 word length (default = 6)

> -e, --compare        : Find d2 score of two sequences (stdout)
> -E, --compare\_all    : Find d2 score of every pair of sequences and
> > outputs successfull matches in -C format (stdout)

> -1, --1 

&lt;int&gt;

        : first sequence to be compared (default = 0)
> -2, --2 

&lt;int&gt;

        : second sequence to be compared (default = 1)

> 

&lt;filename&gt;

 is a fasta formated file of sequences

Notes:
  1. The format of the compact cluster table (-c) is simply that each line
> > contains all the sequences in the cluster. Sequences are numbered from 0
> > and each line is terminated with a ".".

> 2. The format of the extended cluster table (-C) is each line is a tab
> > separated list of:
> > > Sequence1        Sequence2       ReverseComplement       d2Score

> > where Sequence1 and Sequence2 are the numbers of the sequences in the
> > file (starting at 0), ReverseComplement indicates whether the match was
> > with the reverse complement (0=no, 1=yes) and d2Score indicates the d2
> > score between the sequences.

> 3. Each time sd\_clust finishes an analysis it appends the log file
> > (sd\_clust.log). The line it appends the following info:
> > > f\_name w\_thresh w\_size sd\_thresh skip\_val secs heur\_pass d1\_pass d2\_pass

> > where f\_name is the fasta file's name, w\_thresh is the heuristic's
> > word count threshold, w\_size is the heuristic's word size, sd\_thresh
> > is the standard deviation threshold, skip\_val specifies every
> > skip\_valth word to check, secs is the time to complete, heur\_pass
> > is the number of potential matches the heuristic found, d1 pass is
> > is the number of times the d1 heuristic passes and d2\_pass is the
> > number of matches that were found.

Examples:

> sd\_clust -c my\_seqences.fasta
> > This will use the default d2 and sd options and output the compact cluster
> > table to standard output.


> sd\_clust -C -Q mysequences.fasta
> > This will output the extended cluster table showing all overlaps between
> > ESTs.


> sd\_clust -c -t 20 -w 16 -s 2.0 -S 3 -n my\_sequences.fasta
> > This changes the heuristic to use a word length of 16, a word count
> > threshold of 20, a standard deviation threshold of 2 and  it only checks
> > every 3rd word, but doesn't check the reverse complement. Compact cluster
> > table is output.


> sd\_clust -e -1 34 -2 83 my\_sequences.fasta
> > This will calculate the d2 score between sequence #34 and sequence #83.
> > Note that sequences are numbered from 0.