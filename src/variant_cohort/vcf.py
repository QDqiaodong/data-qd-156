"""Strict bounded FASTA and biallelic/multiallelic VCF 4.2 reader."""
import re

_SEQ = re.compile('[ACGTN]+\\Z')
_POS = re.compile('[1-9][0-9]*\\Z')
_SAMPLE = re.compile('[A-Za-z0-9_.-]+\\Z')
_COLUMNS = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

def reference(text):
    if not isinstance(text,str) or len(text)>2_000_000:raise ValueError('FASTA exceeds limit')
    result={};name=None
    for line in text.splitlines():
        if line.startswith('>'):
            name=line[1:]
            if not re.fullmatch(r'[A-Za-z0-9_.-]+',name) or name in result:raise ValueError('invalid/duplicate FASTA contig')
            result[name]=''
        else:
            if name is None or not re.fullmatch('[ACGTN]+',line):raise ValueError('expected uppercase FASTA sequence')
            result[name]+=line
    if not result or any(not s for s in result.values()):raise ValueError('empty reference')
    return result

def _rows(text,ref,fmt,parse_cells,max_samples=200,max_rows=5000,max_alts=16):
    """Walk one bounded VCF, validating every shared structural rule.

    parse_cells(cell, alt_count) parses one sample cell for the declared format.
    """
    if not isinstance(text,str) or len(text)>2_000_000:raise ValueError('VCF exceeds limit')
    lines=text.splitlines()
    if not lines or lines[0]!='##fileformat=VCFv4.2':raise ValueError('VCFv4.2 header required')
    header=None;samples=[];out=[];seen=set()
    for line in lines[1:]:
        if line.startswith('##'):
            if header is not None:raise ValueError('metadata after column header')
            continue
        if line.startswith('#'):
            h=line.split('\t')
            if header is not None or h[:9]!=_COLUMNS:raise ValueError('invalid VCF columns')
            samples=h[9:]
            if not 1<=len(samples)<=max_samples or len(set(samples))!=len(samples) or any(not _SAMPLE.fullmatch(s) for s in samples):raise ValueError('invalid samples')
            header=h;continue
        if header is None:raise ValueError('missing VCF header')
        row=line.split('\t')
        if len(row)!=len(header) or len(out)>=max_rows:raise ValueError('VCF row width/count')
        chrom,pos,rid,base,alt,qual,filt,info,declared=row[:9]
        if chrom not in ref or not _POS.fullmatch(pos) or not _SEQ.fullmatch(base):raise ValueError('invalid locus/ref')
        pos=int(pos);alts=alt.split(',')
        if len(alts)>max_alts or len(set(alts))!=len(alts) or any(not _SEQ.fullmatch(a) or a==base for a in alts):raise ValueError('only distinct sequence ALT supported')
        if ref[chrom][pos-1:pos-1+len(base)]!=base:raise ValueError('reference mismatch')
        if filt not in ('PASS','.') or declared!=fmt or qual!='.' or info!='.':raise ValueError('unsupported FILTER/FORMAT/QUAL/INFO')
        key=(chrom,pos,base,tuple(alts))
        if key in seen:raise ValueError('duplicate variant record')
        seen.add(key);calls=[]
        for cell in row[9:]:
            calls.append(parse_cells(cell,len(alts)))
        out.append((chrom,pos,base,alts,calls))
    if header is None:raise ValueError('VCF header missing')
    return samples,out

def _gt_cell(cell,alt_count):
    if not re.fullmatch(r'(?:[0-9]+|\.)(?:[/|](?:[0-9]+|\.))?',cell):raise ValueError('haploid/diploid GT required')
    parts=re.split(r'[/|]',cell);alleles=[None if x=='.' else int(x) for x in parts]
    if any(a is not None and a>alt_count for a in alleles):raise ValueError('GT allele out of range')
    return alleles

def records(text,ref):
    return _rows(text,ref,'GT',_gt_cell)

def _gt_cell_diploid(cell,alt_count):
    if not re.fullmatch(r'(?:[0-9]+|\.)(?:[/|])(?:[0-9]+|\.)',cell):raise ValueError('diploid GT required')
    copies=tuple(None if x=='.' else int(x) for x in re.split(r'[/|]',cell))
    if any(x is not None and x>alt_count for x in copies):raise ValueError('GT allele out of range')
    return copies

def diploid_records(text,ref):
    """Bounded GT-only reader: every sample cell is exactly two allele copies."""
    return _rows(text,ref,'GT',_gt_cell_diploid,max_samples=32,max_rows=200,max_alts=4)

def phased_records(text,ref):
    """Bounded GT:PS reader: every cell is a phased diploid GT plus phase set."""
    def cell(value,alt_count):
        fields=value.split(':')
        if len(fields)!=2:raise ValueError('sample cell must contain exactly GT:PS')
        gt,ps=fields
        if not re.fullmatch(r'(?:[0-9]+|\.)(?:\|)(?:[0-9]+|\.)',gt):raise ValueError('diploid phased GT with | required')
        copies=tuple(None if x=='.' else int(x) for x in gt.split('|'))
        if any(x is not None and x>alt_count for x in copies):raise ValueError('GT allele out of range')
        if not _POS.fullmatch(ps):raise ValueError('PS must be a positive integer without leading zeros')
        return copies,int(ps)
    return _rows(text,ref,'GT:PS',cell)
