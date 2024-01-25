from Scripts.Phasing import get_args,PhasingRun
from pathlib import Path


def main():
    args = get_args()
    sampledir = args.sampledir
    P = Path(sampledir)
    ref = args.ref
    bed = args.bed
    PhasingRun(P,ref,bed)
    
if __name__ == "__main__":
    main()