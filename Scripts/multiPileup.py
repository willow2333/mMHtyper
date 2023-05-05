#!/anaconda3/envs/FEALPy/bin python3.8
# -*- coding: utf-8 -*-
# ---
# @Software: PyCharm
# @File: multiPileupPhasing.py
# @Author: willow
# @Site: 
# @Time: 6月 20, 2022
# ---

import pandas as pd

import pysam

import time

cigar_oper = {0: 'M', 1: 'I', 2: 'D', 3: 'N', 4: 'S', 5: 'H', 6: 'P', 7: '=', 8: 'X', 9: 'B'}


class MultiPileup:
    '''Create the IGV-likely files,if the depth of site over 10000 will be select randomly'''
    def bamextract(self, bedlist):
        print(bedlist)
        n1 = 0
        pmlist = bedlist[0]
        region = bedlist[1]
        P = bedlist[2]
        exapndname = bedlist[3]
        print('{} bam start: {}'.format(pmlist, time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))))
        bamfile = pmlist
        samfile = pysam.AlignmentFile(bamfile, 'r')
        id_seq_q = pd.DataFrame()
        bedregion = region.split(':')[1].split('-')
        chromosome = region.split(':')[0]
        start = int(bedregion[0])
        end = int(bedregion[-1])
        id_seq_q['chromosome'] = [chromosome for i in range(start, end + 1)]
        # id_seq_q['db SNP'] = []
        id_seq_q['pos'] = [i for i in range(start, end + 1)]

        for x in samfile.fetch(chromosome, start, end):
            # print(x.reference_name,chromosome,x.reference_start,start,end)
            # if x.reference_name == chromosome and  (start < x.reference_start + 1  < end):
            if n1 < 10000:
                expandseq = {}
                name = x.query_name
                rawcigarlist = x.cigar
                sequence = x.query_sequence
                rstart = x.reference_start + 1
                if rstart > start:
                    expandseq.update(dict([(i, 'na') for i in range(start, rstart)]))
                n = rstart
                m = 0
                if rawcigarlist[0][0] == 4:
                    newsequence = sequence[rawcigarlist[0][1]:]
                    cigarlist = rawcigarlist[1:]
                elif rawcigarlist[0][0] == 5:
                    newsequence = sequence
                    cigarlist = rawcigarlist[1:]
                else:
                    newsequence = sequence
                    cigarlist = rawcigarlist

                for i in range(len(cigarlist)):
                    if n <= end:
                        if cigarlist[i][0] == 0:
                            refpos = range(n, n + cigarlist[i][1])
                            qrpos = range(m, m + cigarlist[i][1])
                            if refpos[-1] >= start:
                                for j in range(len(refpos)):
                                    if refpos[j] <= end and refpos[j] >= start:
                                        expandseq[refpos[j]] = newsequence[qrpos[j]]
                                    else:
                                        continue
                            n += cigarlist[i][1]
                            m += cigarlist[i][1]
                        elif cigarlist[i][0] == 1:
                            if n - 1 >= start:
                                expandseq[n - 1] = ','.join([expandseq[n - 1], newsequence[m:m + cigarlist[i][1]]])
                            m += cigarlist[i][1]

                        elif cigarlist[i][0] == 2:
                            refpos = range(n, n + cigarlist[i][1])
                            for j in range(len(refpos)):
                                if refpos[j] <= end and refpos[j] >= start:
                                    expandseq[refpos[j]] = '*'
                                else:
                                    continue
                            n += cigarlist[i][1]
                        elif cigarlist[i][0] == 4 or cigarlist[i][0] == 5:
                            refpos = range(n, n + cigarlist[i][1])
                            for j in range(len(refpos)):
                                if refpos[j] <= end and refpos[j] >= start:
                                    expandseq[refpos[j]] = 'na'
                                else:
                                    continue
                            n += cigarlist[i][1]
                        else:
                            print(name, x.cigarstring)

                if len(expandseq) < (end - start) + 1:
                    for a in range(n, end + 1):
                        expandseq[a] = 'na'
                    # print(pmlist)
                    # print(len(id_seq_q))
                    # print(len(list(expandseq.values())))
                    id_seq_q['{}_base'.format(name)] = list(expandseq.values())
                else:
                    finalseq = {k: v for k, v in expandseq.items() if k <= end}
                    id_seq_q['{}_base'.format(name)] = list(finalseq.values())
                n1 += 1
            else:
                print('This {} site down sampling !'.format(pmlist))
                break
        print('bam end: {}'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))))
        id_seq_q.to_csv(P / '{}_ExpandSeqQual.txt'.format(exapndname), sep='\t', index=False)
