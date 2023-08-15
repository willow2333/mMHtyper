#!/anaconda3/envs/FEALPy/bin python3.8
# -*- coding: utf-8 -*-
# ---
# @Software: PyCharm
# @File: ForensicPhasing.py
# @Author: willow
# @Site: 
# @Time: 4月 28, 2023
# ---


import argparse
from pathlib import Path
from run import Run
import os

NanoFilt=r'/path/software/Anaconda3/envs/Nano/bin/NanoFilt'
minimap2 = '/path/software/Anaconda3/envs/Forensic/bin/minimap2'
samtools = '/path/software/Anaconda3/envs/Forensic/bin/samtools'

def PhasingRun(P, ref,bed):
    namelist = [i for i in os.listdir(P) if i.endswith('.fq.gz')]
    for name in namelist:
        j = name[0:-6]
        newpath = P / j
        tmp = newpath/'tmp'
        if not os.path.exists(newpath):
            os.mkdir(newpath)
        if not os.path.exists(tmp):
            os.mkdir(tmp)
        print('Filter Starting !')
        os.system('gunzip -c {0} |{2}  -q 10 -l 300 --maxlength 1500|gzip > {1}'.format(P / name,newpath / '{}_filter.fq.gz'.format(j),NanoFilt))
        print('Alignment !')
        os.system('{3} -ax map-ont -t 20 --secondary=no {0} {1} -o {2} 2> /dev/null'.format(ref,newpath / '{}_filter.fq.gz'.format(j),newpath / '{}.sam'.format(j),minimap2))
        os.system('{2} view  -q 50 -Sb {0} >{1}'.format(newpath / '{}.sam'.format(j),newpath / '{}.bam'.format(j),samtools))
        os.system('{2} sort -@6 -O bam -o {0} {1}'.format(newpath / '{}.sort.bam'.format(j),newpath / '{}.bam'.format(j),samtools))
        os.system('{1} index {0}'.format(newpath / '{}.sort.bam'.format(j),samtools))
        os.system("{2} view -H  {0} >{1}".format(newpath/'{}.sam'.format(j),tmp/'{}_uniq.sam'.format(j),samtools))
        os.system("{2} view {0}|awk -F'\t' '{{if($2==0) print $0}}' >>{1}".format(newpath / '{}.sort.bam'.format(j),tmp/'{}_uniq.sam'.format(j),samtools))
        os.system("{2} view {0}|awk -F'\t' '{{if($2==16) print $0}}' >>{1}".format(newpath / '{}.sort.bam'.format(j),tmp / '{}_uniq.sam'.format(j),samtools))
        os.system('{2} view  -Sb {0} >{1}'.format(tmp / '{}_uniq.sam'.format(j),tmp / '{}_uniq.bam'.format(j),samtools))
        os.system('{2} sort -@6 -O bam -o {0} {1}'.format(tmp / '{}_uniq.sort.bam'.format(j),tmp / '{}_uniq.bam'.format(j),samtools))
        os.system('{1} index {0}'.format(tmp / '{}_uniq.sort.bam'.format(j),samtools))
        R = Run()
        af = 0.8
        R.run(tmp,'{}_uniq.sort.bam'.format(j),bed,thread=20, af=af)
        os.system('mv {} {}'.format(tmp/'Phasing.txt' ,newpath))
        os.system('rm -r {}'.format(tmp))
        print('Phasing Successfully !')



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--sampledir",help='The absolute path of fastq files, the fastq files need be gzip and named *.fq.gz')
    parser.add_argument("--ref",help='The absolute path of reference file')
    parser.add_argument("--bed",help='The absolute path of bed file')
    args = parser.parse_args()
    sampledir = args.sampledir
    P = Path(sampledir)
    ref = args.ref
    bed = args.bed
    PhasingRun(P,ref,bed)

