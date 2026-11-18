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

# 跨参考组装投影 v3（--project-assembly）
带 `--project-assembly` 时同一CLI改用v3契约，把一份已按source_fasta核对过的VCF记录经显式chains投影到另一份参考组装，产出可审计的等位数据集；不带开关的v1入口与 `--stitch-phase` 的v2入口、输入格式、字段与结果完全不变。该开关与 `--stitch-phase` 互斥，由argparse互斥组拒绝同时出现。本功能不推断相位、不拼接PS块、不做变异调用、基因注释或临床解读。

外层JSON恰含 `source_fasta`、`target_fasta`、`vcf` 与 `chains`，沿用受限JSON规则（拒绝重复键、非有限数值、未定义字段）。`source_fasta` 与 `target_fasta` 各按v1、v2的FASTA规则校验：唯一简短contig名、非空大写ACGTN序列、每份至少一个非空contig，但二者是独立参考、可有不同contig集合。两份FASTA各<=2,000,000字符，外层原始JSON<=8,000,000字节。`vcf` 按 `source_fasta` 沿用v1读取：固定4.2、制表符分隔、1..200唯一样本、最多5000行、每行1..16个唯一非REF的ACGTN ALT，GT保留单倍体/二倍体、`/`或`|` 分隔、点缺失、索引不得越界。每个数据行须给出所有header样本的GT单元，即使全部缺失，仍参与候选重编号与投影收敛的一致性比较。允许VCF只有合法header：`source_records` 为0，`projected`/`rejected` 为空，`samples` 仍列出header中全部样本。

`chains` 为1..128项，每项恰含 `id`、`source_contig`、`target_contig`、`strand`、`blocks`。`id` 为1..64位 `[A-Za-z0-9_.-]+` 且整批唯一；`source_contig`/`target_contig` 必须在对应FASTA中存在；`strand` 仅 `+` 或 `-`。不要求chain覆盖每个contig。

每个chain的 `blocks` 为1..5000项，整批最多20000项；每项恰含 `source_start`、`target_start`、`length`，均为正整数且不得为bool。两个start均为相应FASTA正向序列的1起点左边界；`length` 在两个参考上相同，且区间不得越界。blocks按 `source_start` 严格递增且不重叠；目标区间在 `+` 链严格递增不重叠，在 `-` 链严格递减不重叠，未满足即整批输入非法。不同chain可以重叠。

一行VCF的完整原始REF区间必须由同一个chain连续覆盖；跨多个block时，覆盖区间在source和target上都须紧邻且方向一致，不能跨任一侧缺口、拼接不同chain或只映POS。`+` 链按正向序列投影，`-` 链反向映射坐标并对REF及每个ALT整体做ACGTN反向互补（N仍为N）；插入和删除均仍为单个REF跨度上的序列等位，不拼接目标缺口。

每个完整覆盖候选先取 `target_fasta` 上的映射REF跨度。只有与固定映射跨度长度相同、且逐字等于转换后REF或任一转换后ALT时才可用；不同长度的ALT可被保留为ALT，但不能凭附近序列强行作为目标REF。允许目标参考恰是原ALT，此时该等位成为新REF、原REF须成为新ALT，其他等位不得丢失（若序列重复导致新集合不一致则整批拒绝）。新ALT按字典序排列，并按序列身份重编号每个样本GT；相位分隔符忽略，已知索引升序且null排最后，保留原有倍性与缺失。`genotypes` 每个元素为整数索引或null组成的数组，长度等于原GT倍性，输出不保留 `/` 或 `|` 分隔符，无缺失时也不得改为字符串。

对每个可用候选同时裁剪REF与全部ALT的共同后缀、再裁共同前缀，每条序列至少保留一位；裁前缀时递增目标position；不做重复序列左对齐。候选等价仅当规范化目标 `(contig,position,ref,alts)` 与所有样本的重编号GT都相同；不能仅按chain数或POS判歧义。因此 `chains` 只含通过该来源唯一候选等价判断的可用chain id，不纳入失败或其他歧义候选。

每条源记录没有完整连续覆盖候选时标为 `unmapped`；有完整覆盖但没有目标参考兼容候选时标为 `reference_mismatch`；去除等价候选后多于一个时标为 `ambiguous`；恰一个时才进入projected。未投影记录保留来源身份与原因；不能因一条失败丢弃整批。非法输入、任何内部不一致或下述合并冲突才整批退出2并保留旧输出原字节。

多条源记录收敛到同一规范化目标 `(contig,position,ref,alts)` 时，只有各样本重编号GT完全相同才合并为一条并合并来源；任一样本冲突须整批拒绝，不得取首条、按频率投票或覆盖。来源身份为原始 `(contig,position,ref,按字典序排列的alts)`，输入中来源身份重复须整批拒绝，即使原ALT枚举顺序不同。

报告顶层恰含 `schema_version=3`、`mode="assembly_projection"`、`source_sha256` 与 `result`；result恰含 `samples`、`source_records`、`projected`、`rejected`。`samples` 按名字排序；`source_records` 为原数据行数。`projected` 各项恰含 `contig`、`position`、`ref`、`alts`、`genotypes`、`sources`、`chains`：`genotypes` 是按samples对齐的GT数组，`sources` 为来源身份对象数组，`chains` 为所有等价可用候选的chain id去重排序并集。`rejected` 各项恰含 `source` 与 `reason`；来源对象恰含 `contig`、`position`、`ref`、`alts`。来源与投影均按 `(contig,position,ref,alts)` 排序，`rejected` 按来源排序，所有position输出整数。

对源记录、样本列、ALT枚举（同步重映射GT）及chains列表置换后，除 `source_sha256` 外报告须完全一致。非法输入、任何冲突退出2并保留已有输出的原字节；输入与输出不得相同，拒绝符号链接输出，成功后原子替换。运行保持离线且无第三方依赖，不下载参考。合成样例见 `fixtures/projection.json` 与 `tests/test_projection.py`，覆盖负链多ALT含缺失的GT重编号、目标参考为原ALT的REF/ALT交换、REF跨紧邻block成功而任一侧有gap则unmapped、两个chain规范化后等价不误报ambiguous、多个不同可用目标保留ambiguous、参考不兼容与无覆盖原因区分、两源收敛的一致合并及冲突保旧、置换不变性，并保留v1/v2回归。

# 家系最小屏蔽审计 v4（--audit-pedigree）
带 `--audit-pedigree` 时同一CLI改用v4契约，对家系森林做遗传一致性最小屏蔽审计；不带开关的v1入口与 `--stitch-phase`、`--project-assembly` 入口、输入格式、字段与结果完全不变。该开关与另两个开关互斥，由argparse互斥组拒绝同时出现。本功能不推断相位、不做变异调用、基因注释或临床解读。

外层JSON恰含 `fasta`、`vcf`、`families`、`weights`，沿用受限JSON规则（拒绝重复键、非有限数值、未定义字段）。`fasta` 按v1的FASTA规则校验；`vcf` 沿用v1的VCF4.2固定列、参考逐字核对与FORMAT仅GT规则，但限制1..32个样本、最多200条数据行、每行1..4个ALT；FASTA与VCF文本各<=2,000,000字符，外层原始JSON<=8,000,000字节。每个样本单元必须给恰好两个拷贝，以 `/` 或 `|` 分隔，每拷贝为合法等位索引或点；忽略相位分隔，拒绝单倍体（包括单独的点），允许部分或全部缺失。允许只有合法header而无数据行：`source_records` 为0、`records` 为空、`samples` 仍列出全部header样本。

`weights` 是恰好覆盖全部header样本的对象，各值为1..1,000,000的整数，不接受bool。`families` 为0..32个对象，每项恰含 `parents` 与 `children`；`parents` 恰为两个不同样本名，顺序无意义；`children` 为非空且无重复的样本名数组，不含该家庭任何亲本。所有成员须在VCF header中；每个样本最多作为一个家庭的子女，但可再成为其他家庭的亲本；同一无序亲本对只能声明一个家庭，其同胞须列在同一 `children` 中。家庭节点与样本节点是不同类型的身份：对每项建立一个家庭节点，与两亲及各子女各连一条边，所得无向二部图必须为森林（有环即整批拒绝），允许同胞、跨代和无关联样本；不能只检查亲到子的有向无环。

每条记录独立推导全部样本的完整无序二倍体GT：候选等位只来自该行REF及全部ALT，一个完整GT为两个有重数的等位（0/0含两份0）。每个子女的两个等位须能分别取自两亲各一个；同一样本跨家庭使用同一GT，不能逐家庭单独修正。已知调用按多重集匹配：完整GT须含全部已知拷贝及其重数，未知不约束（`0/.` 兼容0/0或0/1，`1/1` 不兼容0/1，`./.` 兼容任意完整GT）。对不兼容GT只移除该样本该行的观察约束并收取一次 `weights[sample]`，其完整GT仍须满足家系传递；权重不按拷贝数、家庭数或传递次数累乘。各记录独立，不使用跨位点相位、频率先验或性别推断。

在满足全部家庭传递约束的赋值中求全体样本屏蔽代价之和的全局最小值，精确统计所有最优完整GT赋值，并给出每个样本在所有最优解中出现的GT集合。不同赋值按全体样本GT区分；亲本传递路径和同一赋值的重复见证不重复计数（如双亲0/1×0/1、子女0/1只计1个赋值）。`required_masks` 为所有最优赋值都不匹配原调用的样本，`possible_masks` 为至少一个最优赋值不匹配的样本（包含 `required_masks`）；GT支持须为全局可达到的集合，不是各家庭局部候选的并集。不把两个不同的全体GT赋值因屏蔽集合相同而折叠。

报告顶层恰含 `schema_version=4`、`mode="pedigree_audit"`、`source_sha256` 与 `result`；result恰含 `samples`、`source_records`、`records`。`samples` 为全部样本名字典序，`source_records` 为VCF数据行数。`records` 每项恰含 `contig`、`position`、`ref`、`alts`、`minimum_cost`、`optimal_assignments`、`required_masks`、`possible_masks`、`genotypes`。位点身份保留原始contig/POS/REF且将ALT按字典序排列，不裁剪或左对齐；同一 `(contig,position,ref,排序后alts)` 重复整批拒绝，其余重叠记录各自分析。输入GT按序列身份重编号到输出ALT顺序。`minimum_cost` 和 `optimal_assignments` 为精确整数，后者必须为正，不得取模、近似或把大整数转浮点。两个mask数组按样本名排序；`genotypes` 为与 `samples` 对齐的数组，每项是该样本的GT支持列表，各GT为两个按数值升序的整数索引，支持列表去重后按二元组字典序排序。`records` 按 `(contig,position,ref,alts)` 排序，position为整数。

对VCF行、样本列、ALT枚举（同步重映射GT）、families列表及亲本/子女顺序置换后，除 `source_sha256` 外报告须完全一致。非法输入退出2并保留已有输出的原字节；输入与输出不得相同，拒绝符号链接输出，成功后原子替换。运行保持离线且无第三方依赖。合成样例见 `fixtures/pedigree.json` 与 `tests/test_pedigree.py`，覆盖规范并列例（双亲0/0、子女0/1：代价1、赋值5、三者皆可屏蔽）、传递见证去重（0/1×0/1→0/1计1）、跨代共享缺失样本使两个局部可满足家庭整体非零、同胞联合约束、代价并列的全部支持与精确计数、ALT换序重编号、部分缺失重数、多等位、孤立缺失样本组合计数、多独立家庭计数乘积、32样本长链与宽同胞家族、小例独立穷举核对、非法家庭环/重复亲本对/多重子女身份/非法权重/单倍体及各项规模越界，并保留v1/v2/v3回归。

