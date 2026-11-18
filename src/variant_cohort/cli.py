DESCRIPTION='Reference-checked VCF cohort allele datasets'
import argparse, json, sys
from pathlib import Path
from .common import read, save
from .model import analyze
from .phase import stitch

def run(source,stitch_phase=False):
    data,digest=read(source)
    if stitch_phase:return {'schema_version':2,'mode':'phase_stitch','source_sha256':digest,'result':stitch(data)}
    return {'schema_version':1,'source_sha256':digest,'result':analyze(data)}

def main():
    parser=argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--stitch-phase',action='store_true',help='stitch fragment-local phased blocks into global blocks')
    args=parser.parse_args()
    try:
        if args.source.resolve()==args.output.resolve():raise ValueError('input and output must differ')
        report=run(args.source,args.stitch_phase);save(args.output,report)
    except (ValueError,TypeError,KeyError,OSError,OverflowError,RecursionError) as error:
        print(str(error),file=sys.stderr);return 2
    print(json.dumps(report,sort_keys=True,ensure_ascii=False,allow_nan=False));return 0

if __name__=='__main__':raise SystemExit(main())
