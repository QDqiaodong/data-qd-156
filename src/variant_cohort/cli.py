DESCRIPTION='Reference-checked VCF cohort allele datasets'
import argparse, json, sys
from pathlib import Path
from .common import read, save
from .model import analyze
from .pedigree import audit
from .phase import stitch
from .project import project
from .replicate import audit_replicates as run_replicate_audit

def run(source,stitch_phase=False,project_assembly=False,audit_pedigree=False,audit_replicates=False):
    data,digest=read(source)
    if audit_replicates:return {'schema_version':5,'mode':'replicate_audit','source_sha256':digest,'result':run_replicate_audit(data)}
    if audit_pedigree:return {'schema_version':4,'mode':'pedigree_audit','source_sha256':digest,'result':audit(data)}
    if project_assembly:return {'schema_version':3,'mode':'assembly_projection','source_sha256':digest,'result':project(data)}
    if stitch_phase:return {'schema_version':2,'mode':'phase_stitch','source_sha256':digest,'result':stitch(data)}
    return {'schema_version':1,'source_sha256':digest,'result':analyze(data)}

def main():
    parser=argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    modes=parser.add_mutually_exclusive_group()
    modes.add_argument('--stitch-phase',action='store_true',help='stitch fragment-local phased blocks into global blocks')
    modes.add_argument('--project-assembly',action='store_true',help='project reference-checked records onto another assembly')
    modes.add_argument('--audit-pedigree',action='store_true',help='audit pedigree genotype consistency with minimal masking')
    modes.add_argument('--audit-replicates',action='store_true',help='offline pairwise replicate-candidate data-quality screen')
    args=parser.parse_args()
    try:
        if args.source.resolve()==args.output.resolve():raise ValueError('input and output must differ')
        report=run(args.source,args.stitch_phase,args.project_assembly,args.audit_pedigree,args.audit_replicates);save(args.output,report)
    except (ValueError,TypeError,KeyError,OSError,OverflowError,RecursionError) as error:
        print(str(error),file=sys.stderr);return 2
    print(json.dumps(report,sort_keys=True,ensure_ascii=False,allow_nan=False));return 0

if __name__=='__main__':raise SystemExit(main())
