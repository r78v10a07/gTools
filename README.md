# gTools

This code was developed by Roberto Vera Alvarez (r78v10a07@gmail.com) and 
Leonardo Mariño-Ramírez (l.marino.ramirez@gmail.com) to calculate 
the TPM values for introns and exons from a SAM files and a GTF file.

The program was tested using genes.gtf file from hg19 (UCSC).  

Papers used:

1. Wagner GP, Kin K, Lynch VJ. Measurement of mRNA abundance using RNA-seq data: 
RPKM measure is inconsistent among samples. Theory Biosci. 2012 Dec;131(4):281-5. 
doi: 10.1007/s12064-012-0162-3. https://www.ncbi.nlm.nih.gov/pubmed/22872506

2. The SAM/BAM Format Specication Working Group. Sequence Alignment/Map Format 
Specication. http://samtools.github.io/hts-specs/SAMv1.pdf

Install:

For GNU/Linux and MAC OSX with Makefile installed.

Type make for compilation.

Files:

    gTools/dist/Release/libgtools.a
    TPMCalculator/dist/Release/TPMCalculator
    TPMOutPutMerger/dist/Release/TPMOutPutMerger
 