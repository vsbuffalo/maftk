# Validation examples

This is just public note-taking to ensure this is working properly.

Region `chr21:5010905-5010955` has two blocks that overlap it.

```
$ cargo run  query --intersect-only chr21:5010905-5010955
    Finished `dev` profile [unoptimized + debuginfo] target(s) in 0.41s
     Running `target/debug/maftk query --intersect-only 'chr21:5010905-5010955'`
 INFO maftk::binary: Query elapsed 1.12428175s, 2 overlapping alignment blocks found.
[block 0 | score 6871]
                       0
hg38          5010904 + g
panTro4       5010904 - g
[block 1 | score 39400]
                       0         10        20        30        40
hg38          5010905 + tacatttcagggagatgcgaactgcaggtggaatcagaaaacagtacacg
panTro4       5010905 - tacatttcagggagatgcgaactgcaggtggaatcagaaaacagtacacg
calJac3       5010905 - tacatttcagagagacggaaatcgcagatagaatcagaaatcagtacGCG

Legend:
A Match
A Mismatch
- Indel
Coordinates shown as: species start_position [strand] (positions 0-based indexing)

```

The first block is the first basepair at this position in the human genome.
From BLAST:

```
browser new tab details YourSeq    51     1    51    51   100.0%  chr21  +     5010905   5010955     51

0000001 gtacatttcagggagatgcgaactgcaggtggaatcagaaaacagtacac 0000050
>>>>>>> |||||||||||||||||||||||||||||||||||||||||||||||||| >>>>>>>
5010905 gtacatttcagggagatgcgaactgcaggtggaatcagaaaacagtacac 5010954

0000051 g 0000051
>>>>>>> | >>>>>>>
5010955 g 5010955
```

which matches the coordinates (these are 1-based).

Looking at the alignment from maftk, there is 100% identity with panTro4 and
hg38, but 10 substitutions with calJac3. That's 10 / 50 basepairs or a
substitution rate of 20%. Now if we run both species pairs with `maftk stats`:

```
$ maftk stats --species 'hg38,panTro4' chr21:5010905-5010955
    Finished `dev` profile [unoptimized + debuginfo] target(s) in 0.33s
     Running `target/debug/maftk stats --species hg38,panTro4 'chr21:5010905-5010955'`
 INFO maftk::binary: Found 2 blocks in range chr21:5010904-5010955
 INFO maftk::binary: Total elapsed 1.12822675s
chrom	start	end	ref_aligned_start	ref_aligned_end	hg38_panTro4_subst_rate	hg38_panTro4_gap_rate	hg38_panTro4_num_subst	hg38_panTro4_num_single_gaps	hg38_panTro4_num_double_gaps	hg38_panTro4_valid_positions	hg38_panTro4_total_positions
chr21	5010904	5010905	5010828	5010905	0	0	0	0	0	1	1
chr21	5010905	5010955	5010905	5011097	0	0	0	0	0	50	50

$ maftk stats --species 'hg38,calJac3' chr21:5010905-5010955
    Finished `dev` profile [unoptimized + debuginfo] target(s) in 0.16s
     Running `target/debug/maftk stats --species hg38,calJac3 'chr21:5010905-5010955'`
 INFO maftk::binary: Found 2 blocks in range chr21:5010904-5010955
 INFO maftk::binary: Total elapsed 1.023960834s
chrom	start	end	ref_aligned_start	ref_aligned_end	calJac3_hg38_subst_rate	calJac3_hg38_gap_rate	calJac3_hg38_num_subst	calJac3_hg38_num_single_gaps	calJac3_hg38_num_double_gaps	calJac3_hg38_valid_positions	calJac3_hg38_total_positions
chr21	5010904	5010905	5010828	5010905	NA	NA	NA	NA	NA	NA	NA
chr21	5010905	5010955	5010905	5011097	0.2	0	10	0	0	50	50

```

In this test, the results are correct.

## Validation 2

This is a test with two blocks. Region is `chr21:5011558-5011561`.


```
$ matftk --verbosity info query --intersect-only chr21:5011558-5011561
 INFO maftk::binary: Query elapsed 1.136147s, 2 overlapping alignment blocks found.
 INFO maftk::binary: total query region length: 5
[block 0 | score 1826]
                          0
hg38          5011557 +   T
panTro4       5011557 -   T
ponAbe2       5011557 -   T
rheMac3       5011557 +   C
macFas5       5011557 +   C
chlSab2       5011557 -   C
calJac3       5011557 -   C
saiBol1       5011557 +   c
tupChi1       5011557 -   T
musFur1       5011557 -   T
ailMel1       5011557 +   T
odoRosDiv1    5011557 +   T
[block 1 | score 4428]
                          0
hg38          5011558 +   GCA
panTro4       5011558 -   GCG
ponAbe2       5011558 -   GCA
rheMac3       5011558 +   GCG
macFas5       5011558 +   GCG
chlSab2       5011558 -   GCG
calJac3       5011558 -   TCG
saiBol1       5011558 +   tcg
tupChi1       5011558 -   ACC
vicPac2       5011558 +   GCT
musFur1       5011558 -   GTA
ailMel1       5011558 +   GTA
odoRosDiv1    5011558 +   GTA

Legend:
A Match
A Mismatch
- Indel
Coordinates shown as: species start_position [strand] (positions 0-based indexing)
```

There is one substitution between hg38 and panTro4, in four overlapping basepairs.

Looking at the stats, we see both overlapping blocks (they have not been merged
since they are technically different blocks).

```
$ matftk  stats --species 'hg38,panTro4' chr21:5011558-5011561  | cut -f 1,2,3,6,9,10
    Finished `dev` profile [unoptimized + debuginfo] target(s) in 0.39s
     Running `target/debug/maftk stats --species hg38,panTro4 'chr21:5011558-5011561'`
 INFO maftk::binary: Found 2 blocks in range chr21:5011557-5011561
 INFO maftk::binary: Total elapsed 1.109409083s
chrom	start	end	hg38_panTro4_num_subst	hg38_panTro4_valid_positions	hg38_panTro4_total_positions
chr21	5011557	5011558	0	1	1
chr21	5011558	5011561	1	3	3

```

This confirms what we expect: first block (one basepair) has no differences,
second (three bases) has one difference.

The overlapping raw MAF entries:

```
$ zgrep  -A 10 5011557  chr21.maf.gz

s hg38.chr21          5011557 1 +  46709983 T
s panTro4.chr21       2548936 1 -  32799110 T
i panTro4.chr21       C 0 C 0
s ponAbe2.chr21       2639860 1 -  48394510 T
i ponAbe2.chr21       C 0 I 1
s rheMac3.chr3        2486774 1 + 198365852 C
i rheMac3.chr3        C 0 C 0
s macFas5.chr3        2509800 1 + 192294377 C
i macFas5.chr3        C 0 C 0
s chlSab2.chr2        2375720 1 -  90373283 C
i chlSab2.chr2        C 0 C 0

$ zgrep  -A 10 5011558 chr21.maf.gz
s hg38.chr21          5011558 3 +  46709983 GCA
s panTro4.chr21       2548937 3 -  32799110 GCG
i panTro4.chr21       C 0 C 0
s ponAbe2.chr21       2639862 3 -  48394510 GCA
i ponAbe2.chr21       I 1 C 0
s rheMac3.chr3        2486775 3 + 198365852 GCG
i rheMac3.chr3        C 0 C 0
s macFas5.chr3        2509801 3 + 192294377 GCG
i macFas5.chr3        C 0 C 0
s chlSab2.chr2        2375721 3 -  90373283 GCG
i chlSab2.chr2        C 0 C 0
```


## Validation 3

I've annotated the output below and split up stats columns onto multiple lines for clarity.

```
vsb@ponderosa /space/s1/vsb/projects/human_bgs/data/hg38/subrates main*
base ❯ maftk query chr1:10918-10986 --intersect-only
 INFO maftk::binary: Query elapsed 2.416644519s, 1 overlapping alignment blocks found.
 INFO maftk::binary: total *query* region length: 70
[block 0 | score 34237]
                       0         10        20        30        40        50        60
hg38            10917 + gagaggcgcaccgcgccggcgcaggcgcagagacacatgctagcgcgtccaggggtggaggcgtggcgc
panTro4         10917 - gagaggcgcaccgcgccggcgcag------agacacatactagcgcgtcctgggg-ggaggtgcggcgc
                                                              s           s          s s
                                                xxxxxx                         x
s: sub, x: gap

Legend:
A Match
A Mismatch
- Indel
Coordinates shown as: species start_position [strand] (positions 0-based indexing)

vsb@ponderosa /space/s1/vsb/projects/human_bgs/data/hg38/subrates main*
base ❯ maftk stats chr1:10918-10986 --species hg38,panTro4
 INFO maftk::binary: Found 1 blocks in range chr1:10917-10986
 INFO maftk::binary: Total elapsed 2.409640852s
chrom   start   end     
chr1    10917   10986

ref_aligned_start       ref_aligned_end 
10917                   11396   

hg38_panTro4_num_subst  hg38_panTro4_num_single_gaps    hg38_panTro4_num_double_gaps
4                       7                               0    


hg38_panTro4_valid_positions  hg38_panTro4_total_positions
62                            69
```

SO this looks right.
