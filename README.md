Mononucleotide repeat finder
============================

`mnrf` is a simple program to find mononucleotide repeats of specified length in a FASTA file with one or more tracks.

Prerequisites
-------------

* gcc >= 4.3


Compilation
-----------

```
./configure
make 
make install
```


Usage
-----

```
./mnrf

    Usage: mnrf <nrepeats> <input.fa> <output.bed>

```

Run the test case by:

```
./mnrf 5 test/test.fa output.bed
diff output.bed test/answer.bed
```

(No output indicates success.)

