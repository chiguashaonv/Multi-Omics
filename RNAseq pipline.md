# RNAseq pepline

## 数据下载与处理
### 将数据从sra格式转化为fastq.gz格式
对于双端测序，用`fastq-dump`命令的`--split-3`选项
```
fastq-dump --split-3 -A xxx.sra \
--gzip \
-O outpath 
  ```
  
### 检测数据质量QC
```
fastqc -o outpath --noextarct xxx.fastq.gz
```

### 自动化处理数据
可同时处理双端数据
```
fastp -i xxx.1.fastq.gz -o xxx.1.out.fastq.gz -I xxx.2.fastq.gz -O xxx.2.out.fastq.gz
```
