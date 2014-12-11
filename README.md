Mononucleotide repeat finder
============================

`mnrf` is a simple program to find mononucleotide repeats of specified size
(number times of repeat) in a FASTA file with one or more tracks.

Author:    David J. H. Shih  <djh.shih@gmail.com>
Version:   0.1
Date:      2013-11-18
License:   GPL >= 3


Prerequisites
-------------

* gcc >= 4.3


Compilation
-----------

Edit the `Makefile` as necessary. Then compile.

$ make 


Usage
-----

$ ./mnrf
Usage: mnrf <nrepeats> <input.fa> <output.bed>

Run the test case by:

$ ./mnrf 5 test/test.fa output.tsv
$ diff output.tsv test/answer.tsv

(No output indicates success.)

