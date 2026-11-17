# Reference-checked cohort alleles v1
离线科研数据整理者输入恰含fasta/vcf两个文本的JSON；FASTA和VCF各<=2MB，外层JSON<=8MB。
FASTA唯一简短contig名、非空大写ACGTN序列。VCF固定4.2、制表符分隔，1..200唯一样本，最多5000行、每行1..16个唯一非REF的ACGTN ALT；仅GT格式、QUAL/INFO为点、FILTER为PASS或点。不解释任意其他VCF字段或符号变异。
位置1起点，原始REF须逐字匹配参考；GT仅单倍体/二倍体，以/或|分隔，缺失为点，等位索引不得越界。部分缺失时AN仍计已知等位，样本剂量为null；全缺失频率为null。phasing不影响剂量。
每ALT独立输出AC/AN精确分数、按样本剂量；先裁公共后缀再前缀且REF/ALT各至少保留一位。不是重复序列左对齐，规范化后身份碰撞整批拒绝，避免擅自合并独立记录。
合成样例chr1:2 C>T/G：AN5，T AC2、AF2/5，G AC2、AF2/5，S3剂量null。空记录但有效header允许。非临床解读，不做变异调用、基因注释或大型基因组流式处理。
