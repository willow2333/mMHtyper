#!/anaconda3/envs/FEALPy/bin python3.8
# -*- coding: utf-8 -*-
# ---
# @Software: PyCharm
# @File: multiPhasing.py
# @Author: willow
# @Site: 
# @Time: 6月 23, 2022
# ---


import pandas as pd

from tqdm import tqdm, trange
import time
from multiprocessing import Pool
#from spoa import poa

class MultiPhasing:

    def Genotypes(self, expandfiles):
        pileupdf = self.Depth(expandfiles)
        return pileupdf

    def Phasing(self, parameter):
        name, regiondf, expandfiles, sitessample, sitesdbsnp, sitesref,callvariantdf = parameter
        altdict = dict(zip(callvariantdf['chr'] + ':' + callvariantdf['pos'].astype('str'), callvariantdf['Alt']))
        depthdict = dict(zip(callvariantdf['chr'] + ':' + callvariantdf['pos'].astype('str'), callvariantdf['depth']))
        expandf = pd.read_csv(expandfiles, sep='\t')
        phaselist = []
        rawphasedict = {}
        # print('name:{}'.format(name))
        if len(expandf.columns) > 2:
            phasingdf = pd.DataFrame()
            newexpandf = expandf[
                (expandf['chromosome'].isin(regiondf['chromosome'].tolist())) & (
                    expandf['pos'].isin(regiondf['Position (GRCh38)'].tolist()))]

            cols1 = [x for x in newexpandf.columns if len(newexpandf[newexpandf[x].isin(['na'])]) != 0]

            newexpandf = newexpandf.drop(cols1, axis=1)
            # print(len(cols1),len(newexpandf.columns))
            if len(newexpandf) == len(regiondf) and len(newexpandf.columns) >102:
                phasingdf['chr'] = newexpandf['chromosome']
                phasingdf['pos'] = newexpandf['pos']
                phasingdf['REF'] = [
                    sitesref[
                        phasingdf['chr'].tolist()[i] + ':' + phasingdf['pos'].astype('str').tolist()[i]] for i
                    in
                    range(len(phasingdf))]
                phasingdf['sample'] = [
                    sitessample[
                        phasingdf['chr'].tolist()[i] + ':' + phasingdf['pos'].astype('str').tolist()[i]] for i
                    in
                    range(len(phasingdf))]
                phasingdf['db SNP'] = [
                    sitesdbsnp[
                        phasingdf['chr'].tolist()[i] + ':' + phasingdf['pos'].astype('str').tolist()[i]] for i
                    in
                    range(len(phasingdf))]
                phasingdf['genotypes'] =[
                    altdict[
                        phasingdf['chr'].tolist()[i] + ':' + phasingdf['pos'].astype('str').tolist()[i]] for i
                    in
                    range(len(phasingdf))]
                phasingdf['depth'] = [
                    depthdict[
                        phasingdf['chr'].tolist()[i] + ':' + phasingdf['pos'].astype('str').tolist()[i]] for i
                    in
                    range(len(phasingdf))]

                for i in newexpandf.columns :
                    if i not in ['chromosome','pos']:
                        phaselist.append(''.join([i[0] for i in newexpandf[i].tolist()]))

                # rawphasedict = {i:phaselist.count(i) for i in list(set(phaselist)) if '*' not in i}
                rawphasedict = {i: phaselist.count(i) for i in list(set(phaselist))}
                rawphasedictsorted = sorted(rawphasedict.items(), key=lambda dict: dict[1], reverse=True)


                for j1 in range(len(rawphasedictsorted)):
                    phase1 = rawphasedictsorted[j1][0]
                    correct1 = ['yes' for i in range(len(phase1)) if phase1[i] in list(phasingdf['genotypes'].tolist()[i])]
                    if correct1.count('yes') == len(phasingdf):
                        phasingdf['phase1'] = [i for i in phase1]
                        del rawphasedictsorted[j1]
                        break
                if 'phase1' in phasingdf.columns:
                    genotypes = [i.replace(' ','').split(',') for i in phasingdf['genotypes'].tolist()]
                    phasing1 = phasingdf['phase1'].tolist()
                    # print([i[0] for i in rawphasedictsorted])
                    phase2 = []
                    for j2 in range(len(phasingdf)):
                        p1 = phasing1[j2]
                        geno = genotypes[j2]
                        # if p1 in geno:
                        if p1 == geno[0]:
                            phase2.append(geno[1])
                        else:
                            phase2.append(geno[0])
                    phasingdf['phase2'] = phase2
                    # print(''.join(phase2))
                    # if ''.join(phase2) in [i[0] for i in rawphasedictsorted]:
                    #     phasingdf['phase2'] = phase2
                    # else:
                    #     'the sites genotypes wrong!'
                    #     phasingdf['phase2'] = [i for i in rawphasedictsorted[0][0]]

                    return phasingdf
                else:
                    print('phase1 is wrong !')
            else:
                print('The depth of {} sites are not complete !'.format(name))
        else:
            print('{} sites no reads !'.format(name))



    def phasingcount(self,n,new_total_phasinglist):
        beforemergephasing = []
        for p1 in new_total_phasinglist:
            if len(p1[0][0]) - 2 == 0:
                beforemergephasing.append([(i[0],i[1]) for i in p1])
            else:
                phasing = []
                for i in range(len(p1)):
                    if i != len(p1)-1:
                        phasing.append((p1[i][0][:2],p1[i][1][:2]))
                    else:
                        tailphase = p1[i]
                        for j in range(len(tailphase[0])):
                            if j < len(tailphase[0])-1:
                                # print(j)
                                phasing.append((tailphase[0][j:j+2],tailphase[1][j:j+2]))
                beforemergephasing.append(phasing)
        # print(beforemergephasing)

        mergephasing = []
        for n1 in range(n - 1):
            microphsing = []
            for p1 in beforemergephasing:
                microphsing.append(p1[n1])
            mergephasing.append(microphsing)
        # print(mergephasing)

        genotypes = []
        for i in range(len(mergephasing)):
            if i ==0:
                tmp0 = [j[0][0] for j in mergephasing[i]] +[j[1][0] for j in mergephasing[i]]
                tmp0dict = dict([(j,tmp0.count(j)) for j in set(tmp0)])
                sorttmp0dict = sorted(tmp0dict.items(), key=lambda dict: dict[1], reverse=True)
                if len(list(sorttmp0dict)) >1:
                    if sorttmp0dict[0][1]/sorttmp0dict[1][1] > 4:
                        genotypes.append(sorttmp0dict[0][0])
                    else:
                        genotypes.append(''.join([sorttmp0dict[0][0],sorttmp0dict[1][0]]))
                else:
                    genotypes.append(sorttmp0dict[0][0])

                tmp1 = [j[0][1] for j in mergephasing[i]]  + [j[1][1] for j in mergephasing[i]]+ [j[0][0] for j in mergephasing[i+1]] + [j[1][0] for j in mergephasing[i+1]]
                tmp1dict = dict([(j,tmp1.count(j)) for j in set(tmp1)])
                sorttmp1dict = sorted(tmp1dict.items(), key=lambda dict: dict[1], reverse=True)
                if len(list(sorttmp1dict)) > 1:
                    if sorttmp1dict[0][1] / sorttmp1dict[1][1] >= 4:
                        genotypes.append(sorttmp1dict[0][0])
                    else:
                        genotypes.append(''.join([sorttmp1dict[0][0],sorttmp1dict[1][0]]))
                else:
                    genotypes.append(sorttmp1dict[0][0])
            elif i == len(mergephasing)-1:

                tmp1 = [j[0][1] for j in mergephasing[i]] + [j[1][1] for j in mergephasing[i]]
                tmp1dict = dict([(j, tmp1.count(j)) for j in set(tmp1)])
                sorttmp1dict = sorted(tmp1dict.items(), key=lambda dict: dict[1], reverse=True)
                if len(list(sorttmp1dict)) > 1:
                    if sorttmp1dict[0][1] / sorttmp1dict[1][1] > 4:
                        genotypes.append(sorttmp1dict[0][0])
                    else:
                        genotypes.append(''.join([sorttmp1dict[0][0], sorttmp1dict[1][0]]))
                else:
                    genotypes.append(sorttmp1dict[0][0])

            else:
                tmp1 = [j[0][1] for j in mergephasing[i]]  + [j[1][1] for j in mergephasing[i]]+ [j[0][0] for j in mergephasing[i+1]] + [j[1][0] for j in mergephasing[i+1]]
                tmp1dict = dict([(j,tmp1.count(j)) for j in set(tmp1)])
                sorttmp1dict = sorted(tmp1dict.items(), key=lambda dict: dict[1], reverse=True)
                if len(list(sorttmp1dict)) > 1:
                    if sorttmp1dict[0][1] / sorttmp1dict[1][1] >= 4 :
                        genotypes.append(sorttmp1dict[0][0])
                    else:
                        genotypes.append(''.join([sorttmp1dict[0][0],sorttmp1dict[1][0]]))
                else:
                    genotypes.append(sorttmp1dict[0][0])

        # print(genotypes)
        phase1 = []
        phase2 = []
        # print(mergephasing)
        for a in range(len(mergephasing)):
            phase = mergephasing[a]
            newphase = [i[0] for i in phase] + [i[1] for i in phase]
            newphase1dict = dict([(i, newphase.count(i)) for i in newphase])
            newphase1dictsort = [i[0] for i in
                                 sorted(newphase1dict.items(), key=lambda dict: dict[1], reverse=True)]
            if len(newphase1dictsort)>1:
                number = 0
                for n1 in range(len(newphase1dictsort) - 1):
                    site1 = newphase1dictsort[n1]
                    for n2 in range(n1, len(newphase1dictsort)):
                        site2 = newphase1dictsort[n2]
                        if (sorted(set(list(site1[0] + site2[0]))) == sorted(list(genotypes[a]))) and (sorted(
                                set(list(site1[1] + site2[1]))) == sorted(list(genotypes[a + 1]))):
                            phase1.append(site1)
                            phase2.append(site2)
                            number += 1
                            break
                if number == 0:
                    phase1.append(newphase1dictsort[0])
                    phase2.append(newphase1dictsort[1])

            else:
                phase1.append(newphase1dictsort[0])
                phase2.append(newphase1dictsort[0])

        # print(phase1, phase2)
        newphase1 = [phase1[0][0],phase1[0][1]]+ [phase1[i][1] for i in range(1,len(phase1))]
        newphase2 = [phase2[0][0],phase2[0][1]]+ [phase2[i][1] for i in range(1,len(phase2))]
        # print(newphase1,newphase2)
        return newphase1,newphase2





    def Depth(self, pileupmatrix):
        print('{}:csv load start!'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))))
        newdf = pd.DataFrame()
        pos = []
        depth = []
        chrs = []
        A = []
        T = []
        C = []
        G = []
        deletion = []
        with open(pileupmatrix, 'r') as f:
            lines = f.readlines()
            for i in tqdm(lines):
                if not i.startswith('chromosome'):
                    tmp = i.strip().split('\t')
                    align = [tmp[i] for i in range(2, len(tmp)) if tmp[i] != 'na']
                    if len(align) > 0:
                        aligns = self.Count(len(align), align)
                        A.append(aligns['A'])
                        T.append(aligns['T'])
                        C.append(aligns['C'])
                        G.append(aligns['G'])
                        deletion.append(aligns['*'])
                        depth.append(len(align))
                        pos.append(int(tmp[1]))
                        chrs.append(tmp[0])
                    else:
                        continue
        print('{}:csv load successfully!'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))))
        newdf['chr'] = chrs
        newdf['pos'] = pos
        newdf['depth'] = depth
        newdf['A'] = A
        newdf['T'] = T
        newdf['C'] = C
        newdf['G'] = G
        newdf['*'] = deletion
        return newdf

    def Count(self, dep, baselist):
        base_count = {}

        for i in range(len(baselist)):
            b = baselist[i]
            if len(b) == 1:
                if b == '*':
                    if b not in list(base_count.keys()):
                        base_count[b] = 1

                    else:
                        base_count[b] += 1

                else:
                    if b not in list(base_count.keys()):
                        base_count[b] = 1

                    else:
                        base_count[b] += 1

            else:
                newb = ''.join(b.replace("'", '').replace('[', '').replace(']', '').split(',')).replace(' ', '')

                if newb not in list(base_count.keys()):
                    base_count[newb] = 1

                else:
                    base_count[newb] += 1

        newbase_count = {'A': 0, 'T': 0, 'C': 0, 'G': 0, '*': 0}
        for k, v in base_count.items():
            newbase_count[k[0]] += v
        for key in list(newbase_count.keys()):
            if newbase_count['*'] / dep > 0.75:
                newbase_count[key] = float('{:.4f}'.format(newbase_count[key] / dep))
            else:
                if key == '*':
                    newbase_count[key] = 0
                else:
                    newbase_count[key] = float('{:.4f}'.format(newbase_count[key] / (dep - newbase_count['*'])))
        return newbase_count
