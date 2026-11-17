"""Strict bounded FASTA and biallelic/multiallelic VCF 4.2 reader."""
import re

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

def records(text,ref):
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
            if header is not None or h[:9]!=['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']:raise ValueError('invalid VCF columns')
            samples=h[9:]
            if not 1<=len(samples)<=200 or len(set(samples))!=len(samples) or any(not re.fullmatch(r'[A-Za-z0-9_.-]+',s) for s in samples):raise ValueError('invalid samples')
            header=h;continue
        if header is None:raise ValueError('missing VCF header')
        row=line.split('\t')
        if len(row)!=len(header) or len(out)>=5000:raise ValueError('VCF row width/count')
        chrom,pos,rid,base,alt,qual,filt,info,fmt=row[:9]
        if chrom not in ref or not re.fullmatch('[1-9][0-9]*',pos) or not re.fullmatch('[ACGTN]+',base):raise ValueError('invalid locus/ref')
        pos=int(pos);alts=alt.split(',')
        if len(alts)>16 or len(set(alts))!=len(alts) or any(not re.fullmatch('[ACGTN]+',a) or a==base for a in alts):raise ValueError('only distinct sequence ALT supported')
        if ref[chrom][pos-1:pos-1+len(base)]!=base:raise ValueError('reference mismatch')
        if filt not in ('PASS','.') or fmt!='GT' or qual!='.' or info!='.':raise ValueError('unsupported FILTER/FORMAT/QUAL/INFO')
        key=(chrom,pos,base,tuple(alts))
        if key in seen:raise ValueError('duplicate variant record')
        seen.add(key);calls=[]
        for gt in row[9:]:
            if not re.fullmatch(r'(?:[0-9]+|\.)(?:[/|](?:[0-9]+|\.))?',gt):raise ValueError('haploid/diploid GT required')
            parts=re.split(r'[/|]',gt);alleles=[None if x=='.' else int(x) for x in parts]
            if any(a is not None and a>len(alts) for a in alleles):raise ValueError('GT allele out of range')
            calls.append(alleles)
        out.append((chrom,pos,base,alts,calls))
    if header is None:raise ValueError('VCF header missing')
    return samples,out
