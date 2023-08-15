#!/anaconda3/envs/FEALPy/bin python3.8
# -*- coding: utf-8 -*-
# ---
# @Software: PyCharm
# @File: multiSplit.py
# @Author: willow
# @Site: 
# @Time: 6月 23, 2022
# ---

import os
from multiprocessing import Pool

import pandas as pd
import subprocess

'''
According the bed file spilt the big bam file !!!
'''

samtools = '/path/software/Anaconda3/envs/Forensic/bin/samtools'
CONSENT = '/path/software/Anaconda3/envs/Forensic/bin/CONSENT-correct'
seqkit = '/path/software/Anaconda3/envs/Forensic/bin/seqkit'
ref = '/path/Reference/reference_genome/hg38_puire/hg.38.fa'
minimap2 = '/path/software/Anaconda3/envs/Forensic/bin/minimap2'


class MultiSplit:
    def __init__(self, bed):
        print('bed file load !')
        bedfiles = bed
        # bedfiles = r'/data/qinliu/Sever/Forensic_multiSNP/multiSNP_WF.bed'
        self.beddf = pd.read_csv(bedfiles, sep='\t')
        self.groupdf = self.beddf.groupby(by='MiniHap name')
        self.sitesbed = {}
        for k, v in self.groupdf:
            self.sitesbed[k] = v['chromosome'].tolist()[0] + ':' + '-'.join(
                [v['Position (GRCh38)'].astype('str').tolist()[0], v['Position (GRCh38)'].astype('str').tolist()[-1]])


    ####### consider correction the raw reads #############
    def split(self, tmp):
        P, sortbamfile, newsitesbed = tmp
        print('split bam file creation!')
        splitbamlist = []
        for k, v in newsitesbed.items():
            bamfilepath = P / '{}.bam'.format(k)
            sortbamfilepath = P / '{}.sort.bam'.format(k)
            fastqfile =  P / '{}.fastq'.format(k)
            fastafile = P / '{}.fasta'.format(k)
            os.system('{0} view -hb {1} {2} > {3}'.format(samtools, P/sortbamfile, v, bamfilepath))
            print('{} size is :{}'.format(k,os.path.getsize(bamfilepath)))
            if os.path.getsize(bamfilepath) > 1000:
                if k == 'mh01KK-172' or k =='mh21KK-324' or k == 'mh02KK-031':
                    os.system('{0} sort -@6 -O bam -o {1} {2}'.format(samtools, sortbamfilepath, bamfilepath))
                    os.system('rm -rf {0}'.format(bamfilepath))
                    os.system('{0} index {1}'.format(samtools, sortbamfilepath))
                    splitbamlist.append((sortbamfilepath, k))
                else:
                    os.system('{0} fastq {1} |seqkit sample -n 1000 > {2}'.format(samtools, bamfilepath,fastqfile))
                    os.system('{0} --in {1} --out {2} --type ONT -j 15'.format(CONSENT, fastqfile, fastafile))
                    os.system('{0} -ax map-ont -t 32 {1} {2} | samtools view -hb -F 0x904 > {3} 2> /dev/null'.format(minimap2,ref,fastafile,bamfilepath))
                    os.system('{0} sort -@6 -O bam -o {1} {2}'.format(samtools, sortbamfilepath, bamfilepath))
                    os.system('rm -rf {0}'.format(bamfilepath))
                    os.system('{0} index {1}'.format(samtools, sortbamfilepath))
                    splitbamlist.append((sortbamfilepath, k))
            else:
                os.system('{0} sort -@6 -O bam -o {1} {2}'.format(samtools, sortbamfilepath, bamfilepath))
                os.system('rm -rf {0}'.format(bamfilepath))
                os.system('{0} index {1}'.format(samtools, sortbamfilepath))
                splitbamlist.append((sortbamfilepath, k))
                print('{}:{}'.format(k,v))

            ####
            # splitbamlist.append((sortbamfilepath, k))
        # print(splitbamlist)
        return splitbamlist


    def multi(self, P, sortbamfile):
        thread = 4
        finalsplitbamlist = []
        sitesbedsplit = [dict(zip(list(self.sitesbed.keys())[i:i + 8],
                                  list(self.sitesbed.values())[i:i + 8])) for i in range(0, len(self.sitesbed),
                                                                                             8)]
        # print(sitesbedsplit,len(sitesbedsplit))
        pool = Pool(thread)
        splitbamlist = pool.map(self.split, [[P, sortbamfile,sitesbedsplit[i] ] for i in range(0,len(sitesbedsplit))])
        pool.close()
        pool.join()
        finalsplitbamlist += splitbamlist
        result = [j for i in finalsplitbamlist for j in i]
        return result


