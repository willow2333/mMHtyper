#!/anaconda3/envs/FEALPy/bin python3.8
# -*- coding: utf-8 -*-
# ---
# @Software: PyCharm
# @File: run.py
# @Author: willow
# @Site: 
# @Time: 6月 23, 2022
# ---
import argparse
from multiprocessing import Pool
from multiPileup import *
from multiPhasing import *
# from newmultiPhasing import *
import multiSplit
from pathlib import Path
import pandas as pd

import warnings

warnings.filterwarnings('ignore')

class Run():
    def Splitbam(self,threads, P, bamfilelist, beddict):
        thread = threads
        newbamlist = [bamfilelist[i:i + thread] for i in range(0, len(bamfilelist), thread)]
        for bam in newbamlist:
            pool = Pool(len(bam))
            multi = MultiPileup()
            pool.map(multi.bamextract, [[i[0], beddict[i[1]], P, i[1]] for i in bam])
            pool.close()
            pool.join()


    def Callvaraint(self,P, beddf, namelist,af):
        sitesref = dict(zip(beddf['chromosome'] + ':' + beddf['Position (GRCh38)'].astype('str'), beddf['REF']))
        sitessample = dict(zip(beddf['chromosome'] + ':' + beddf['Position (GRCh38)'].astype('str'), beddf['MiniHap name']))
        sitesdbsnp = dict(zip(beddf['chromosome'] + ':' + beddf['Position (GRCh38)'].astype('str'), beddf['db SNP']))
        variantmergedf = pd.DataFrame()
        for name in namelist:
            exapandfile = P / '{}_ExpandSeqQual.txt'.format(name)
            variantdf = MultiPhasing().Genotypes(exapandfile)
            if len(variantdf) > 0:
                variantmergedf = pd.concat([variantmergedf, variantdf], axis=0)
        variantmergedf = variantmergedf[(variantmergedf['chr'].isin(beddf['chromosome'].tolist())) & (
            variantmergedf['pos'].isin(beddf['Position (GRCh38)'].tolist())) & (variantmergedf['depth'] > 100)]

        variantmergedf['REF'] = [
            sitesref[variantmergedf['chr'].tolist()[i] + ':' + variantmergedf['pos'].astype('str').tolist()[i]] for i in
            range(len(variantmergedf))]
        variantmergedf['db SNP'] = [
            sitesdbsnp[variantmergedf['chr'].tolist()[i] + ':' + variantmergedf['pos'].astype('str').tolist()[i]] for i in
            range(len(variantmergedf))]

        variantmergedf['sample'] = [
            sitessample[variantmergedf['chr'].tolist()[i] + ':' + variantmergedf['pos'].astype('str').tolist()[i]] for i in
            range(len(variantmergedf))]

        # variantmergedf.to_csv(P / 'sites_pileup.txt', sep='\t', index=False)

        callvariantdf = pd.DataFrame()
        Alt = []
        a1_ratio = []
        a2_ratio = []
        # Alt2 = []
        for i in range(len(variantmergedf)):
            varaintdict = {'A': variantmergedf['A'].tolist()[i],
                           'T': variantmergedf['T'].tolist()[i],
                           'C': variantmergedf['C'].tolist()[i],
                           'G': variantmergedf['G'].tolist()[i],
                           '*': variantmergedf['*'].tolist()[i]}
            varaintdictsort = dict(sorted(varaintdict.items(), key=lambda dict: dict[1], reverse=True))
            maxvariant = list(varaintdictsort.keys())[0]

            if maxvariant != variantmergedf['REF'].tolist()[i]:
                if varaintdictsort[maxvariant] > af:
                    Alt.append(','.join([maxvariant, maxvariant]))
                    a1_ratio.append(varaintdictsort[maxvariant])
                    a2_ratio.append(varaintdictsort[maxvariant])
                else:
                    try:
                        alratio = varaintdictsort[list(varaintdictsort.keys())[0]] / sum(list(varaintdictsort.values()))
                        if  alratio > af:
                            Alt.append(','.join(list([maxvariant, maxvariant])))
                            a1_ratio.append(alratio)
                            a2_ratio.append(alratio)
                        else:
                            Alt.append(','.join([maxvariant, list(varaintdictsort.keys())[1]]))
                            a1_ratio.append(varaintdictsort[maxvariant])
                            a2_ratio.append(varaintdictsort[list(varaintdictsort.keys())[1]])
                    except:
                        Alt.append(','.join(list([maxvariant, maxvariant])))
                        a1_ratio.append(varaintdictsort[maxvariant])
                        a2_ratio.append(varaintdictsort[maxvariant])
            else:
                try:
                    alratio = varaintdictsort[list(varaintdictsort.keys())[0]] / sum(list(varaintdictsort.values()))
                    # print(alratio)
                    if alratio >af:
                        Alt.append(','.join([maxvariant, maxvariant]))
                        a1_ratio.append(varaintdictsort[maxvariant])
                        a2_ratio.append(varaintdictsort[maxvariant])
                    else:
                        Alt.append(','.join([maxvariant, list(varaintdictsort.keys())[1]]))
                        a1_ratio.append(varaintdictsort[maxvariant])
                        a2_ratio.append(varaintdictsort[list(varaintdictsort.keys())[1]])
                except:
                    Alt.append(','.join(list([maxvariant, maxvariant])))
                    a1_ratio.append(varaintdictsort[maxvariant])
                    a2_ratio.append(varaintdictsort[maxvariant])

        callvariantdf['chr'] = variantmergedf['chr']
        callvariantdf['pos'] = variantmergedf['pos']
        callvariantdf['REF'] = variantmergedf['REF']
        callvariantdf['db SNP'] = variantmergedf['db SNP']
        callvariantdf['depth'] = variantmergedf['depth']
        callvariantdf['Alt'] = Alt
        callvariantdf['a1_ratio'] = a1_ratio
        callvariantdf['a2_ratio'] = a2_ratio

        callvariantdf.to_csv(P / 'variants.txt', sep='\t', index=False)
        return callvariantdf


    def CallPhasing(self,threads, P, groupbeddf, beddf, callvariantdf):
        sitessample = dict(zip(beddf['chromosome'] + ':' + beddf['Position (GRCh38)'].astype('str'), beddf['MiniHap name']))
        sitesdbsnp = dict(zip(beddf['chromosome'] + ':' + beddf['Position (GRCh38)'].astype('str'), beddf['db SNP']))
        sitesref = dict(zip(beddf['chromosome'] + ':' + beddf['Position (GRCh38)'].astype('str'), beddf['REF']))
        finalphasingdf = pd.DataFrame()
        thread = threads
        name_expand_list = [(k, v) for k, v in groupbeddf]
        expandlist = [name_expand_list[i:i + thread] for i in range(0, len(name_expand_list), thread)]
        for expand in expandlist:
            pool = Pool(len(expand))
            Phase = MultiPhasing()
            parameter = []
            for i in expand:
                expandfiles = P / '{}_ExpandSeqQual.txt'.format(i[0])
                regiondf = i[1]
                parameter.append([i[0], regiondf, expandfiles, sitessample, sitesdbsnp, sitesref, callvariantdf])
            phasingdf = pool.map(Phase.Phasing, parameter)
            pool.close()
            pool.join()
            for result in phasingdf:
                finalphasingdf = pd.concat([finalphasingdf, result], axis=0)
        finalphasingdf.to_csv(P / 'Phasing.txt', sep='\t', index=False)


    def run(self,P, bamfile, bed,thread,af):
        Split = multiSplit.MultiSplit(bed)
        beddict = Split.sitesbed
        beddf = Split.beddf
        groupbeddf = Split.groupdf
        bamlist = Split.multi(P,bamfile)
        namelist = [i[1] for i in bamlist]
        self.Splitbam(thread, P, bamlist, beddict)
        callvariantdf = self.Callvaraint(P, beddf, namelist,af)
        self.CallPhasing(thread, P, groupbeddf, beddf, callvariantdf)


