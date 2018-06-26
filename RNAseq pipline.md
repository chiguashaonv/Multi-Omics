# RNAseq pepline

## 数据下载与处理
### 将数据从sra格式转化为fastq.gz格式
对于双端测序，用`fastq-dump`命令的`--split-3`选项
```powershell
fastq-dump --split-3 -A xxx.sra \
--gzip \
-O outpath 
```
  
### 检测数据质量QC
```powershell
fastqc -o outpath --noextarct xxx.fastq.gz
```

### 自动化处理数据
可同时处理双端数据
```powershell
fastp -i xxx.1.fastq.gz -o xxx.1.out.fastq.gz -I xxx.2.fastq.gz -O xxx.2.out.fastq.gz
```

## 序列比对

### 收集并了解不同的比对软件
[比较不同的比对软件的优势以及劣势](https://blog.csdn.net/qazplm12_3/article/details/76700045)。常用的比对软件包括：STAR, TopHat和HISAT2
>STAR具有最高比例的在基因组上有唯一比对位置的reads，尤其是对读长为300 nt的MCF7样品也有最高的比对率。与TopHat和HISAT2不同，STAR只保留双端reads都比对到基因组的序列，但对低质量的比对 (允许更多的错配碱基和soft-clip事件) 容忍度高。这一点在长reads >(MCF7-300)样品中的体现更为明显。TopHat则不允许soft-clip事件。

>在比对速度方面，HISAT2比STAR快2.5倍，比TopHat快大约100倍。

>soft-clip事件: 即reads末端存在低质量碱基或接头导致比对不上的, STAR会自动尝试截去未比对部分，只保留比对上的部分。

接下来，我们对这三种软件都安装并且进行尝试：[STAR安装](https://github.com/alexdobin/STAR)，[TopHat2安装](http://ccb.jhu.edu/software/tophat/index.shtml)，[HISAT2安装](https://ccb.jhu.edu/software/hisat2/index.shtml)

在了解这三种软件的创始人以及使用目的时，我发现TopHat2和HISAT2都是约翰霍普金斯大学计算生物学中心发表的软件，而且再TopHat2的首页上也已经提出了：
>Please note that TopHat has entered a low maintenance, low support stage as it is now largely superseded by HISAT2 which provides the same core functionality (i.e. spliced alignment of RNA-Seq reads), in a more accurate and much more efficient way.

因此，我们下面只尝试STAR和HISAT2这两种方法
**STAR**
[使用手册](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
[生信菜鸟团的使用说明](http://www.bio-info-trainee.com/727.html)
1. 生成genome Index
*--runThreadN* 线程数

*--runMode* 运行模式：genomeGenerate

*--genomeDir* /path/to/genomeDir

*--genomeFastaFiles* /path/to/genome/fasta1 /path/to/genome/fasta2 ...

--sjdbGTFfile /path/to/annotations.gtf

*--sjdbOverhang* ReadLength-1

*--runThreadN* option defines the number of threads to be used for genome generation, it has to be set to the number of available cores on the server node.

*--runMode* genomeGenerate option directs STAR to run genome indices generation job.

*--genomeDir* specifies path to the directory (henceforth called ”genome directory” where the
genome indices are stored. This directory has to be created (with mkdir) before STAR run
and needs to writing permissions. The filesystem needs to have at least 100GB of disk space
available for a typical mammalian genome. It is recommended to remove all files from the
genome directory before running the genome generation step. This directory path will have to
be supplied at the mapping step to identify the reference genome.

*--genomeFastaFiles* specified one or more FASTA files with the genome reference sequences.
Multiple reference sequences (henceforth called chromosomes) are allowed for each fasta file.
You can rename the chromosomes names in the chrName.txt keeping the order of the chromosomes
in the file: the names from this file will be used in all output alignment files (such as
.sam). The tabs are not allowed in chromosomes names, and spaces are not recommended.

*--sjdbGTFfile* specifies the path to the file with annotated transcripts in the standard GTF
format. STAR will extract splice junctions from this file and use them to greatly improve
accuracy of the mapping. While this is optional, and STAR can be run without annotations,
using annotations is highly recommended whenever they are available.

*--sjdbOverhang* specifies the length of the genomic sequence around the annotated junction
to be used in constructing the splice junctions database. Ideally, this length should be equal
to the ReadLength-1, where ReadLength is the length of the reads. For instance, for Illumina
2x100b paired-end reads, the ideal value is 100-1=99. In case of reads of varying length, the
ideal value is max(ReadLength)-1. In most cases, a generic value of 100 will work as
well as the ideal value.

```powershell
STAR \
 --runThreadN 30	\
 --runMode genomeGenerate	\
 --genomeDir /home/liuke/reference/genomeDir \
 --genomeFastaFiles /home/liuke/reference/hg19.fa	\
 --sjdbGTFfile /home/liuke/reference/hg19.gtf	\
 --sjdbOverhang 99
```

2. 进行序列比对
*--runThreadN* NumberOfThreads

*--genomeDir* /path/to/genomeDir

*--readFilesIn* /path/to/read1 [/path/to/read2 ]

*--runThreadN* option defines the number of threads to be used for genome generation, it has
to be set to the number of available cores on the server node.

*--genomeDir* specifies path to the genome directory where genome indices where generated.

*--readFilesIn* name(s) (with path) of the files containing the sequences to be mapped (e.g. RNA-seq FASTQ files). If using Illumina paired-end reads, the read1 and read2 files have to be supplied. STAR can process both FASTA and FASTQ files. Multi-line (i.e. sequence split in multiple lines) FASTA file are supported. If the read files are compressed, use the

*--readFilesCommand* UncompressionCommand option, where UncompressionCommand is the un-compression command that takes the file name as input parameter, and sends the uncompressed output to stdout. For example, for gzipped files (*.gz) use --readFilesCommand zcat OR --readFilesCommand gzip -c. For bzip2-compressed files, use --readFilesCommand bzip2 -c.

```powershell
STAR \
--runThreadN 30 \
--genomeDir genomeDir \
--readFilesIn fastq1 fastq2 \
--outFileNamePrefix ${mapping_dir}/${sample_id}. 
```
























