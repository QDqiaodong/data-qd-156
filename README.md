# data-qd-156 · VCF群体变异数据整理

Reference-checked VCF cohort allele datasets。业务输入、计算和当前边界详见 [数据契约](docs/contract.md)。

## 安装与演示

Python >=3.11，运行和测试仅用标准库；先确认 `python3 --version`。
每仓建立独立环境，构建工具锁版本。首次安装无公共缓存时可能联网；业务运行不访问网络。

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-build.lock
.venv/bin/python -m pip install --no-build-isolation --no-deps .
.venv/bin/python -m pip check
.venv/bin/variant-cohort --source fixtures/input.json --output run/report.json
.venv/bin/variant-cohort --source fixtures/phase.json --output run/phase.json --stitch-phase
.venv/bin/variant-cohort --source fixtures/projection.json --output run/projection.json --project-assembly
.venv/bin/variant-cohort --source fixtures/pedigree.json --output run/pedigree.json --audit-pedigree
.venv/bin/variant-cohort --source fixtures/replicate.json --output run/replicate.json --audit-replicates
.venv/bin/variant-cohort --source fixtures/cohort_merge.json --output run/cohort_merge.json --merge-cohorts
.venv/bin/variant-cohort --source fixtures/qc_cohort.json --output run/qc_cohort.json --qc-cohort
.venv/bin/python -m unittest discover -s tests -v
```

本地已有wheel时，安装构建工具可加 `--no-index --find-links <公共wheel缓存目录>`。
修改源码后重新安装。本仓不共享其他仓环境，不依赖Docker、外部账号或持续服务。

输入JSON最多8MB，拒绝重复键、非有限数值和未定义字段；源字节SHA256随结果保留。
输出是可独立读取的JSON；成功后原子替换，输入/计算失败退出2并保留旧报告。
输入和输出不得相同，输出符号链接拒绝。单操作者批处理，不承诺并发发布或断电恢复。
同输入重跑不累加。测试使用临时目录，运行产物在忽略的run/。

src/variant_cohort/model.py 为v1业务模型，phase.py 为v2分片相位整理，project.py 为v3跨参考组装投影，pedigree.py 为v4家系最小屏蔽审计，replicate.py 为v5队列样本重复候选筛查，merge.py 为v6多批次联合基因型矩阵，vcf.py 共用受限FASTA/VCF读取；
qc_cohort.py 为v7队列质量筛选，common.py 提供受限JSON及报告发布，cli.py 衔接数据链路并以 `--stitch-phase`、`--project-assembly`、`--audit-pedigree`、`--audit-replicates`、`--merge-cohorts` 或 `--qc-cohort` 选择模式（六者互斥）。
fixtures 为明确自制样例，代码与样例采用MIT许可。

## 业务与样例

# Reference-checked cohort alleles v1
离线科研数据整理者输入恰含fasta/vcf两个文本的JSON；FASTA和VCF各<=2MB，外层JSON<=8MB。
FASTA唯一简短contig名、非空大写ACGTN序列。VCF固定4.2、制表符分隔，1..200唯一样本，最多5000行、每行1..16个唯一非REF的ACGTN ALT；仅GT格式、QUAL/INFO为点、FILTER为PASS或点。不解释任意其他VCF字段或符号变异。
位置1起点，原始REF须逐字匹配参考；GT仅单倍体/二倍体，以/或|分隔，缺失为点，等位索引不得越界。部分缺失时AN仍计已知等位，样本剂量为null；全缺失频率为null。phasing不影响剂量。
每ALT独立输出AC/AN精确分数、按样本剂量；先裁公共后缀再前缀且REF/ALT各至少保留一位。不是重复序列左对齐，规范化后身份碰撞整批拒绝，避免擅自合并独立记录。
合成样例chr1:2 C>T/G：AN5，T AC2、AF2/5，G AC2、AF2/5，S3剂量null。空记录但有效header允许。非临床解读，不做变异调用、基因注释或大型基因组流式处理。

# 分片相位整理 v2（可选 `--stitch-phase`）
不加开关时v1入口、格式与结果不变。加开关后输入改为 `fasta` + `fragments`（1..32个 `{id,vcf}`），
每份VCF为 `GT:PS` 分片相位调用；同一样本、同一contig的局部块按共享规范化等位的已知拷贝方向强制连接，
未知不投票，矛盾整批拒绝。输出 `schema_version=2`、`mode="phase_stitch"` 的确定性全局块，
每个块以最小 (fragment_id, PS) 成员为anchor，成员带 flipped 标志，等位copies为两个0/1/null。
全部分片合计<=5000行、样本并集<=200；置换分片/行/样本/ALT顺序结果不变（source_sha256除外）。
完整规则与验收覆盖见 [数据契约](docs/contract.md)，`fixtures/phase.json` 为链式连接样例。

# 跨参考组装投影 v3（可选 `--project-assembly`，与 `--stitch-phase` 互斥）
输入为 `source_fasta`、`target_fasta`、`vcf` 与 `chains`（1..128条显式链，每条1..5000个正向坐标block）。
VCF按 source_fasta 沿用v1核对；每个完整REF区间须由同一chain连续覆盖（跨block时两侧都须紧邻），
`-` 链整体反向互补（N仍为N），插入/删除仍是单跨度序列等位。目标REF须逐字等于转换后REF或某个ALT，
否则记为 reference_mismatch；无完整覆盖记为 unmapped；多个不同可用目标记为 ambiguous；恰一个才投影。
输出 `schema_version=3`、`mode="assembly_projection"`，含 samples、source_records、projected 与 rejected；
GT按序列身份重编号，忽略相位分隔、保留倍性与缺失。源记录/样本/ALT/chains置换后结果不变（source_sha256除外）。
完整规则与验收覆盖见 [数据契约](docs/contract.md)，`fixtures/projection.json` 为端到端样例。

# 家系最小屏蔽审计 v4（可选 `--audit-pedigree`，与其余开关互斥）
输入改为 `fasta` + `vcf` + `families` + `weights`：1..32个样本、最多200行、每行1..4个ALT，仅接受二倍体GT（允许部分或全部缺失）。
families声明双亲与子女，样本/家庭二部图必须为森林；weights为每样本1..1,000,000的整数屏蔽代价。
每条记录独立求解：在满足全部家系传递约束的完整GT赋值中最小化屏蔽代价之和，精确统计最优赋值数
（同一GT赋值的多个传递见证只计一次），并给出每样本全局GT支持与必屏蔽/可屏蔽集合。
输出 `schema_version=4`、`mode="pedigree_audit"`；行/样本/ALT/家庭置换后结果不变（source_sha256除外）。
完整规则与验收覆盖见 [数据契约](docs/contract.md)，`fixtures/pedigree.json` 为端到端样例。

# 队列样本重复候选筛查 v5（可选 `--audit-replicates`，与其余开关互斥）
输入改为 `fasta` + `vcf` + `minimum_informative_sites` + `max_discordance`：2..100个样本、最多5000行、每行1..16个ALT，仅接受二倍体GT（允许部分或全部缺失，缺失拷贝所在记录不进入该样本对的任一分母）。
每条GT映射为REF/ALT序列的无序二倍体多重集，`/` 与 `|` 等价、ALT枚举顺序不影响身份；记录身份 `(contig,position,ref,排序后alts)` 重复整批拒绝。
逐对统计 complete/informative/concordant/discordant 位点，informative 仅计双方完整且至少一方携带非REF等位的记录；比例判定用精确整数交叉相乘，边界等号属于 `candidate`，证据不足为 `insufficient_evidence`。
输出 `schema_version=5`、`mode="replicate_audit"`，每对含至多8个按规范身份排序的错配例子与省略计数；不生成聚类、不做传递闭包、不推断身份或亲缘。样本列/行/ALT/分隔符置换后结果不变（source_sha256除外）。
完整规则与验收覆盖见 [数据契约](docs/contract.md)，`fixtures/replicate.json` 为端到端样例。

# 多批次联合基因型矩阵 v6（可选 `--merge-cohorts`，与其余开关互斥）
输入改为 `fasta` + `cohorts`（1..32个 `{id,vcf}`，id为1..64位 `[A-Za-z0-9_.-]+` 且全局唯一）：
每份VCF含1..100个唯一样本、最多5000行、每行1..16个ALT，GT允许单倍体或二倍体（`/` 或 `|`，点缺失）；
全部批次样本名互不重复且并集<=200，数据行合计<=5000，允许只有合法header的空VCF批次。
批次以精确 `(contig,position)` 站点身份对齐，同站点REF须逐字相同（冲突整批失败），批次内重复站点整批拒绝，不trim、不左对齐、不合坐标不同的等价indel。
联合 `alts` 取各批次ALT序列的字典序集合并集，每个GT拷贝按碱基序列身份重映射（REF仍为0、点仍为null），保留拷贝顺序、倍性、缺失及 `/` 与 `|` 相位；
样本所在批次未观测该站点时输出null，显式 `0/0` 才是观测参考纯合。
输出 `schema_version=6`、`mode="cohort_merge"`，含 samples、source_cohorts、source_records 与 sites；
批次/VCF行/样本列/ALT枚举置换后结果不变（source_sha256除外），失败退出2并保留旧输出原字节。
完整规则与验收覆盖见 [数据契约](docs/contract.md)，`fixtures/cohort_merge.json` 为端到端样例。

# 队列质量筛选 v7（可选 `--qc-cohort`，与其余开关互斥）
输入为 `fasta` + `vcf` + `policy`：2..100个样本、最多5000行、每行1..16个ALT，仅接受二倍体GT（`/` 或 `|`，允许部分或全部缺失）。
policy 恰含 `min_sample_call_rate`、`min_site_call_rate`、`max_sample_heterozygosity`（各为 `{numerator,denominator}`）与 `min_site_alt_count`（0..100000整数）。
先按规范化等位身份重映射GT并拒绝身份碰撞；样本阶段以全部数据行为分母，只有两拷贝均已知的GT计为called，`called=0` 时杂合率为null且不能通过，
阈值用精确整数交叉相乘、边界等号保留。位点阶段只在保留样本上重新计算 called_samples/call_rate/alt_count/alt_carriers，
先样本后位点、不迭代回补；双失败原因固定取字典序较小者（`high_heterozygosity`/`low_call_rate` 与 `low_alt_count`/`low_call_rate`）。
输出 `schema_version=7`、`mode="cohort_qc"`，含 policy、source_records、samples、sample_metrics、excluded_samples、sites 与 excluded_sites；
比值输出为最简整数分数、零分母为null，GT保留相位分隔符与拷贝顺序。行/样本列/ALT枚举/等价编码置换后结果不变（source_sha256除外）。
完整规则与验收覆盖见 [数据契约](docs/contract.md)，`fixtures/qc_cohort.json` 为端到端样例。
