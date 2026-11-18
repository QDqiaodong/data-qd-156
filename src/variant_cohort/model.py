from fractions import Fraction
from .common import require
from .vcf import reference,records
from .merge import merge_cohorts

def normalize(pos,ref,alt):
    while len(ref)>1 and len(alt)>1 and ref[-1]==alt[-1]:ref,alt=ref[:-1],alt[:-1]
    while len(ref)>1 and len(alt)>1 and ref[0]==alt[0]:pos+=1;ref,alt=ref[1:],alt[1:]
    return pos,ref,alt

def analyze(data):
    require(data,'fasta vcf');ref=reference(data['fasta']);samples,variants=records(data['vcf'],ref);out=[];seen=set()
    for chrom,pos,base,alts,calls in variants:
        an=sum(a is not None for gt in calls for a in gt)
        for index,alt in enumerate(alts,1):
            p,r,a=normalize(pos,base,alt);key=(chrom,p,r,a)
            if key in seen:raise ValueError('normalized allele collision; reconcile source first')
            seen.add(key);ac=sum(x==index for gt in calls for x in gt)
            out.append({'contig':chrom,'position':p,'ref':r,'alt':a,'allele_count':ac,'allele_number':an,
                        'frequency':str(Fraction(ac,an)) if an else None,'dosages':{s:None if None in gt else gt.count(index) for s,gt in zip(samples,calls)}})
    return {'samples':samples,'source_records':len(variants),'alleles':sorted(out,key=lambda x:(x['contig'],x['position'],x['ref'],x['alt']))}
