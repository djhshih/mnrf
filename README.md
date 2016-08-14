Mononucleotide repeat finder
============================

[![Build
Status](https://travis-ci.org/djhshih/mnrf.svg?branch=master)](https://travis-ci.org/djhshih/mnrf)
[![codecov](https://codecov.io/gh/djhshih/mnrf/branch/master/graph/badge.svg)](https://codecov.io/gh/djhshih/mnrf)

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

