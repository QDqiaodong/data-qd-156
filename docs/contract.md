# Reference-checked cohort alleles v1
离线科研数据整理者输入恰含fasta/vcf两个文本的JSON；FASTA和VCF各<=2MB，外层JSON<=8MB。
FASTA唯一简短contig名、非空大写ACGTN序列。VCF固定4.2、制表符分隔，1..200唯一样本，最多5000行、每行1..16个唯一非REF的ACGTN ALT；仅GT格式、QUAL/INFO为点、FILTER为PASS或点。不解释任意其他VCF字段或符号变异。
位置1起点，原始REF须逐字匹配参考；GT仅单倍体/二倍体，以/或|分隔，缺失为点，等位索引不得越界。部分缺失时AN仍计已知等位，样本剂量为null；全缺失频率为null。phasing不影响剂量。
每ALT独立输出AC/AN精确分数、按样本剂量；先裁公共后缀再前缀且REF/ALT各至少保留一位。不是重复序列左对齐，规范化后身份碰撞整批拒绝，避免擅自合并独立记录。
合成样例chr1:2 C>T/G：AN5，T AC2、AF2/5，G AC2、AF2/5，S3剂量null。空记录但有效header允许。非临床解读，不做变异调用、基因注释或大型基因组流式处理。

# 分片相位整理 v2（--stitch-phase）
带 `--stitch-phase` 时同一CLI改用v2契约，把同一样本、同一contig的多个VCF局部相位块按证据拼成确定性全局块；不带此开关的v1入口、输入格式、JSON字段与结果完全不变，两模式互不混用。

外层JSON恰含fasta与fragments，fragments为1..32项，每项恰含id与vcf。id为1..64位 `[A-Za-z0-9_.-]+` 且在整批内唯一。FASTA与每份VCF文本各<=2,000,000字符，外层原始JSON<=8,000,000字节。沿用v1的FASTA规则、VCF4.2固定列、参考核对、每行最多16个ALT及受限JSON规则（拒绝重复键、非有限数值、未定义字段）。全部分片合计最多5000数据行，样本（header名）并集最多200。每份VCF须有合法header，各含1..200个唯一样本，可只有header而无数据行；无观察的样本只进入samples，不制造块。

分片格式仅允许 `FORMAT=GT:PS`，每个样本单元恰含两个冒号字段。GT恰为两个由 `|` 连接的等位索引或点，索引须合法（0..ALT数）；拒绝 `/`、单倍体、缺PS、`GT` 后跟额外字段及PS前导零/非正整数。PS即使在全缺失调用中也必须提供，其身份仅在 (fragment_id, sample, contig) 内有效；不同分片的相同PS不得直接视为同一个块。数据行、样本列、ALT枚举及输入分片均无排序保证，算法只依赖规范化身份，不按行号、ALT序号或原始坐标配对。

每个ALT沿用v1“先裁后缀再前缀”的规范化身份 (contig, position, ref, alt)；同一分片内规范化身份碰撞整批拒绝。跨分片相同身份是相位重叠证据。该ALT对应的每个GT拷贝映为1（此ALT）、0（其他已知等位）或未知（点）：0/1是针对单个ALT的成员指示，其它ALT不当作缺失。

对同一样本、同一contig的任意两个局部块，考察“保持”与“交换”左右拷贝两种相对方向：一种方向仅当所有共享规范化身份中、所有两侧都已知的对应拷贝值相等时才兼容。只剩一种兼容方向才产生强制连接；两种都兼容则保持独立且不猜方向；无兼容方向整批拒绝。未知不是0也不是反对票，不能丢弃纯合或部分缺失调用；没有共享身份时两种方向都兼容。相容性须汇总两个块的全部共享身份后判断，不能因单个位点支持而忽略其他位点冲突。所有强制连接必须能同时满足：即使任意两块单独可兼容，只要多块闭环方向矛盾也整批拒绝。

输出块恰为强制连接的连通分量，不连接不同sample/contig，也不合并没有强制定向证据的块；全缺失的孤立块仍输出。每个输出块以成员中按 (fragment_id字典序, PS整数值) 最小者为anchor，其方向固定为不翻转；其余成员方向由全部强制连接共同决定。按该方向合并每个规范化身份的两份拷贝，已知值填补未知，均未知仍为null；任何已知冲突整批拒绝，不以输入先后覆盖。

报告顶层恰含 `schema_version=2`、`mode="phase_stitch"`、`source_sha256` 与 `result`；result恰含samples（样本并集字典序）、source_fragments（分片数）、source_records（数据行总数）与blocks。每个block恰含sample、contig、anchor、members与alleles；anchor恰含fragment_id与phase_set；members各项恰含fragment_id、phase_set与flipped布尔值；alleles各项恰含position、ref、alt与copies，copies为两个0/1/null。phase_set与position输出整数。members按成员键、alleles按(position,ref,alt)、blocks按(sample,contig,anchor成员键)排序。

对分片、数据行、样本列和ALT枚举同时置换并同步重映射GT后，除source_sha256外报告须完全一致。非法输入或任何冲突退出2并保留已有输出的原字节；输入与输出不得相同，拒绝符号链接输出，成功后原子替换。临床解读、变异调用或基因注释均不在范围内。合成样例见 tests/test_phase.py，覆盖链式连接含一次交换、两两兼容却闭环矛盾、同PS跨分片无重叠仍分离、多ALT换序及可裁剪表示对齐、部分缺失双向兼容不连接、全缺失块保留与置换不变性。
