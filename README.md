# RaPID-Query

Author: Yuan Wei (yuan.wei@Knights.ucf.edu)

Copyright (c) 2022 University of Central Florida and University of Texas Health Science Center at Houston

The tool is free for non-commercial usage. Please use at your own discretion. For commercial use please contact us.

## Run RaPID-Query Program

RaPID-Query program is compiled using GCC 9.1.0 with -Os optimization flag under a 64-bit Unix based operating system.

RaPID-Query program has below parameters:
- p, or panel `<INPUT PANEL FILE>`, where `<INPUT PANEL FILE>` is the input panel path and file name.
- q, or query `<INPUT QUERY FILE>`, where `<INPUT QUERY FILE>` is the input query path and file name.
- g, or genetic `<INPUT GENETIC MAPPING FILE>`, where `<INPUT GENETIC MAPPING FILE>` is the input genetic mapping path and file name.
- m, or match `<OUTPUT MATCH FILE>`, where `<OUTPUT MATCH FILE>` is the output match path and file name.
- lm, or length_marker `<MINIMUM IBD MARKERS IN NUMBER OF SITES>`, where `<MINIMUM IBD MARKERS IN NUMBER OF SITES>` is the minimum IBD markers in number of sites.
- lmh, or length_marker_high `<MINIMUM IBD MARKERS IN NUMBER OF SITES (HIGH)>`, where `<MINIMUM IBD MARKERS IN NUMBER OF SITES (HIGH)>` is the minimum IBD markers in number of sites for high resolution.
- d, or distance `<MINIMUM IBD LENGTH IN CM (LOW)>`, where `<MINIMUM IBD LENGTH IN CM (LOW)>` is the minimum IBD length in centiMorgans for low resolution.
- dh, or distance_high `<MINIMUM IBD LENGTH IN CM (HIGH)>`, where `<MINIMUM IBD LENGTH IN CM (HIGH)>` is the minimum IBD length in centiMorgans for high resolution.
- dg, or distance_gap `<MAXIMUM GAP BETWEEN IBDS IN CM>`, where `<MAXIMUM GAP BETWEEN IBDS IN CM>` is the maximum gap between IBDs in centiMorgans.
- w, or window `<WINDOW SIZE>`, where `<WINDOW SIZE>` is the panel window size for sub sampling.
- r, or run `<NUMBER OF RUNS>`, where `<NUMBER OF RUNS>` is the number of runs (sub-panels).
- c, or count `<MINIMUM NUMBER OF COUNT OF SUCCESSES>`, where `<MINIMUM NUMBER OF COUNT OF SUCCESSES>` is the minimum number of count of successes (hits) as a real match.

An example command of running RaPID-Query program:
```
./RaPID-Query_v1.0 -w 13 -r 5 -c 1 -d 7.0 -lm 700 -dh 1.0 -lmh 100 -dg 2.0 -m ./output_ibds.txt -p ./test_panel.vcf -q ./test_query.vcf -g ./test_map.txt 
```

The command to get the help of the program:
```
./RaPID-Query_v1.0 -h
```
## Input and Output Files
Three input files are required to run RaPID-Query program: the panel file in VCF format, the query file in VCF format, and the genetic map file in HapMap format. More than one individuals can be included in the query file. To perform the query against the panel, the values of *CHROM* field and *POS* field in both panel and query VCF files should match. The genetic map file uses the HapMap format, whose description should be found in the first line of the file. The format of each line starting with the second line contains four tab-delimited fields: *Chromosome*, *Position(bp)*, *Rate(cM/Mb)*, and *Map(cM)*. Note that the value of the *Rate(cM/Mb)* field is not used. If the genetic mapping of the physical position in VCF file is not found, interpolation is used to estimate the genetic distance of such position.

The output file contains the identity by descent (IBD) found by RaPID-Query program in a comma-delimited text format. The first line indicates the format of the file, and the identified IBD segments are listed starting from the second line. The output format contains nine fields: *Query individual id*, *Query individual haplotype id* (0 or 1), *Panel individual id*, *Panel individual haplotype id* (0 or 1), *Physical start position*, *Physical end position*, *Genetic length*, *Site start index* (in VCF file), and *Site end index* (in VCF file).
