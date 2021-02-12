# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 12:13:57 2020

@author: James

converts epitope file, b and t cell predictions into binned csv format 
compatible with the BepiTBR prediction model
"""
#from matplotlib import pyplot as plt
from collections import Counter
import os
from multiprocessing import Pool
import math
from time import time
import argparse


def bEpiLength(fn):
    '''
    get the B cell epitope length and type from file

    Parameters
    ----------
    fn : str
        path to file of epitopes.

    Returns
    -------
    d : dict
        dict of epitope lengths: {<epi_header>: {type: <epi_type>, length: <epi_length>}}.

    '''
    with open(fn) as f:
        l = f.readlines()
    l = [x[:-1] if x[-1] == '\n' else x for x in l]
    d = {l[x]:{'type':l[x].split(sep='_')[0].replace('>', ''), 'length':len(l[x + 1])} for x in range(0, len(l), 2)}
    
    return d

# def plotDir(p):
#     fns = os.listdir(p)
#     fns = [x for x in fns if ('.txt' in x) and (not 'log' in x) and (not 'protein' in x)]
#     for x in fns:
#         print(x)
#         d = bEpiLength(os.path.join(p, x))
#         c = Counter([x['length'] for x in d.values()])
#         plt.figure()
#         plt.bar([y for y in c], [c[y] for y in c])
#         plt.xlabel('epi length')
#         plt.ylabel('number of epitopes')
#         plt.title(x.replace('.txt', ''))

def pickCols(l, colNames):
    #returns only columns included in colNames from a list with headers as first element
    header = l[0]
    inds = [header.index(x) for x in colNames]
    l = [[y[x] for x in inds] for y in l]
    return l

def epiBinPos(headers, values, bEpiPos, tEpiLenngth, bEpiLength, distanceLimit, binSize, cutoff, method):
    '''
    dependency for readEpi
    find the bin the t epitope is located in with repect to b epitope

    Parameters
    ----------
    headers : list
        header line from MHCII prediction results in the form if list.
    values : list of lists
        values in MHCII results, excludes first column (sequences).
    bEpiPos : int
        starting position of b cell epitope.
    tEpiLenngth : int
        length of t cell epitope in peptides.
    bEpiLength : int
        length of b cell epitope in peptides.
    distanceLimit : int
        maximum distance allowed between t and b epitope
    binSize : int
        bin width.
    cutoff : float
        maximum percentile rank of a t epitope to count as a binder.
    method : int
        method to calculate bin number:
            0: simple difference in start position between t epitopes to 
                b epitope divided by bin width
            1: middle to middle distance divided by bin width
            2: minimum distance between t and b epitope divided by bin width,
                t epitopes that overlaps b epitopes are given bin number of
                distanceLimit/binSize + 1
            all bin numbers start at 1            

    Returns
    -------
    outDict : dict
        {<MHCII_allele>: {<bin#>: <#binders>}}.

    '''
    nBins = distanceLimit/binSize
    headerDict = {x: headers[x] for x in range(len(headers))}
    outDict = {}
    for x in values:
        tEpiPos = int(x[0])
        if method == 0:     #simple start to start calculation of distance
            epiDist = tEpiPos - bEpiPos
            if abs(epiDist) > distanceLimit:
                continue
            else:
                binNumber = int(math.floor(epiDist / binSize) + nBins + 1)
                hits = [f'{headerDict[y]}' for y in range(1, len(x)) if x[y] <= cutoff]               
            for allele in hits:
                if allele in outDict:
                    if binNumber in outDict[allele]:
                        outDict[allele][binNumber] += 1
                    else:
                        outDict[allele][binNumber] = 1
                else:
                    outDict[allele] = {binNumber: 1}
        elif method == 1:   #middle to middle distance
            epiDist = tEpiPos + tEpiLenngth/2 - bEpiPos - bEpiLength/2
            if abs(epiDist) > distanceLimit:
                continue
            else:
                binNumber = int(math.floor(epiDist / binSize) + nBins + 1)
                hits = [f'{headerDict[y]}' for y in range(1, len(x)) if x[y] <= cutoff]               
            for allele in hits:
                if allele in outDict:
                    if binNumber in outDict[allele]:
                        outDict[allele][binNumber] += 1
                    else:
                        outDict[allele][binNumber] = 1
                else:
                    outDict[allele] = {binNumber: 1}
        elif method == 2:   #shortest distance between epitopes
            #if t epitope is before b epitope
            if (tEpiPos + tEpiLenngth) < (bEpiPos):
                epiDist = tEpiPos + tEpiLenngth - bEpiPos
            #if t epitope is after b epitope
            elif tEpiPos > (bEpiPos + bEpiLength):
                epiDist = tEpiPos - bEpiPos - bEpiLength
            #if t and b epitopes overlap
            else:
                epiDist = 0
            if abs(epiDist) > distanceLimit:
                continue
            else:
                binNumber = int(math.floor(epiDist / binSize) + nBins + 1)
                hits = [f'{headerDict[y]}' for y in range(1, len(x)) if x[y] <= cutoff]               
            for allele in hits:
                if allele in outDict:
                    if binNumber in outDict[allele]:
                        outDict[allele][binNumber] += 1
                    else:
                        outDict[allele][binNumber] = 1
                else:
                    outDict[allele] = {binNumber: 1}
    return outDict
                    
            
def readEpi(t_data_dir, b_data_dir, epi_header, allelesTypes, tSoftwareList, bSoftwareList, boundary, bEpiLen, binSize, cutoff, method, refHeader=None):
    '''
    read t and b cell predictions for a given epitope and convert the data to a list

    Parameters
    ----------
    t_data_dir : str
        path to directory of MHCII prediction results.
    b_data_dir : str
        path to directory of b cell prediction results.
    epi_header : str
        header line from data file.
    allelesTypes : list of str
        list of types of alleles to aggregate hits.
    tSoftwareList : list of str
        list of MHCII binding prediction software.
    bSoftwareList : list of str
        list of b cell epitope prediction software.
    boundary : int
        maximum distance allowed between t and b epitope.
    bEpiLen : int
        length of b cell epitope in peptides.
    binSize : int
        bin width.
    cutoff : float
        maximum percentile rank of a t epitope to count as a binder.
    method : int
        method to calculate bin number:
            0: simple difference in start position between t epitopes to 
                b epitope divided by bin width
            1: middle to middle distance divided by bin width
            2: minimum distance between t and b epitope divided by bin width,
                t epitopes that overlaps b epitopes are given bin number of
                distanceLimit/binSize + 1
            all bin numbers start at 1            
    refHeader : list of str, optional
        expected header for the output data. The default is None.

    Returns
    -------
    lOut: list
        data in the format:
            [<b_cell_software1_result>, <b_cell_software2_result> ... <b_cell_software#n_result>,
            <allele1_MHCII_software1_bin1_result>, <allele1_MHCII_software1_bin2_result> ... <allele1_MHCII_software1_bin#n_result> ...
            <allele#n_MHCII_software1_bin1_result> ... <allele#n_MHCII_software1_bin#n_result> ...
            <allele1_MHCII_software2_bin1_result> ... <allele#n_MHCII_software2_bin#n_result> ...
            ... <allele#n_MHCII_software#n_bin#n_result>].

    '''
    tDirPath = os.path.join(t_data_dir, epi_header[1:])
    bDirPath = os.path.join(b_data_dir, epi_header[1:])
    #check if prediction data folders exist
    if not os.path.isdir(tDirPath):
        e = f't cell data directory for "{epi_header}" not found in {t_data_dir}'
        print(e)
        return 1
    if not os.path.isdir(bDirPath):
        e = f'b cell data directory for "{epi_header}" not found in {b_data_dir}'
        print(e)
        return 1
    bPaths = [os.path.join(bDirPath, f'{x}.out') for x in bSoftwareList]
    if refHeader:
        outHeader = ['Epitope'] + [x for x in bSoftwareList]
    lOut = [epi_header[1:]]
    try:
        for p in bPaths:
            with open(p) as f:
                bOut = f.readline().replace('\n', '')
                if not bOut:
                    e = f'b cell prediction {p} is empty; skipping epitope'
                    print(e)
                    return 1
                lOut.append(bOut)
        infoPath = os.path.join(tDirPath, 'information.txt')
        with open(infoPath) as f:
            l = f.readlines()
        for x in l:
            if 'Bepi_pos' in x:
                bEpiPos = int(x.split(sep='\t')[-1].replace('\n', ''))
                break
            else:
                bEpiPos = 'na'
    except:
        print(f'data file(s) not found for {epi_header}')
        return 1
    if bEpiPos == 'na':
        e = f'b cell epitope position for "{epi_header}" not found in {tDirPath}'
        print(e)
        return 1
    for software in tSoftwareList:
        dataPath = os.path.join(tDirPath, software + '_out.txt')
        if not os.path.isfile(dataPath):
            e = f' {software} data file for "{epi_header}" not found in {tDirPath}'
            print(e)
            return 1
        else:
            with open(dataPath) as f:
                l = f.readlines()
            l = [x.split(sep=' ') for x in l]
            if len(l) > 1: 
                tEpiLength = len(l[1][0])
                headers = l[0]
                at = allelesTypes + ['Pos']
                alleles = [x for x in headers if [y for y in at if y in x]]
                l = pickCols(l, alleles)
                headers = l[0]
                values = [[float(y) for y in x] for x in l[1:]]
                d = epiBinPos(headers, values, bEpiPos, tEpiLength, bEpiLen, boundary, binSize, cutoff, method)
            else:
                d = []
            totalBins = int(boundary/binSize) * 2
            if method == 2:
                totalBins += 1
            for t in allelesTypes:
                if d:
                    ct = [Counter(d[allele]) for allele in d if t in allele]
                    if ct:
                        cOut = ct[0]
                        for c in ct[1:]:
                            cOut.update(c)
                    else:
                        cOut = {}
                else:
                    cOut = {}
                dTmp = {f'{t}_{software}_{binNumber}': cOut[binNumber] if binNumber in cOut else 0 for binNumber in range(1, (totalBins + 1))}
                lOut += [dTmp[f'{t}_{software}_{binNumber}'] for binNumber in range(1, (totalBins + 1))]
                if refHeader:
                    outHeader += [f'{t}_{software}_{binNumber}' for binNumber in range(1, (totalBins + 1))]
    if refHeader:
        if outHeader == refHeader:
            return lOut
        else:
            e = f'finished header for {epi_header} ({outHeader}) is not the same as reference header ({refHeader})'
            print(e)
            return 1
    return lOut

def str2list(s):
    if ' ' in s:
        s.replace(' ', '')
    l = s.split(sep=',')
    return l
    
if __name__ == '__main__':

    # sourceDir = 'sample_processed_files/sample_source'
    # outDir = 'sampleOut'
    # t_data_dir = 'sample_processed_files'
    # b_data_dir = 'sample_processed_files/bcell'
    # epiLengthLimit = 30
    # allelesTypes = ['D', 'DRB', 'DP', 'DQ']
    # tSoftwareList = ["MixMHC2pred","netMHCIIpan"]
    # bSoftwareList = ['bepipred1.0', 'bepipred2.0', 'LBEEP']
    # boundary = 100
    # binSize = 10
    # cutoff = 20
    # method = 2
    #plotDir(p)  
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', type=str, 
                        help='epitope file')
    parser.add_argument('-o', type=str, 
                        help='output directory')
    parser.add_argument('-t', type=str, 
                        help='directory of t cell predictions')
    parser.add_argument('-b', type=str, 
                        help='directory of b cell predictions')
    parser.add_argument('--l', type=int, default=30,
                        help='maximum peptide length of b epitopes; default: 30')
    parser.add_argument('--a', type=str, default='D,DRB,DP,DQ',
                        help='types of MHCII alleles to include in output, comma seperated; default: "D, DRB, DP, DQ"')
    parser.add_argument('--tl', type=str, default="MixMHC2pred,netMHCIIpan",
                        help='names of MHCII binidng prediction software, comma seperated; default: "MixMHC2pred, netMHCIIpan"')
    parser.add_argument('--bl', type=str, default='bepipred1.0,bepipred2.0,LBEEP',
                        help='names of b cell epitope prediction software, comma seperated; default: "bepipred1.0, bepipred2.0, LBEEP"')    
    parser.add_argument('--d', type=int, default=180,
                        help='maximum distance allowed between t and b epitope; default: 180')
    parser.add_argument('--bs', type=int, default=20,
                        help='bin width; default: 20')
    parser.add_argument('--c', type=float, default=1.5,
                        help='maximum percentile rank of a t epitope to count as a binder; default: 1')
    parser.add_argument('--m', type=int, default=1,
                        help='''method to calculate bin number:\n\t
            0: simple difference in start position between t epitopes to 
                b epitope divided by bin width\n\t
            1: middle to middle distance divided by bin width\n\t
            2: minimum distance between t and b epitope divided by bin width,\n\t\t
                t epitopes that overlaps b epitopes are given bin number of\n\t\t
                distanceLimit/binSize + 1\n\t
            all bin numbers start at 1\n\n\t
            default: 1''')
    args = parser.parse_args()
    sourceFn = args.s
    outDir = args.o
    t_data_dir = args.t
    b_data_dir = args.b
    epiLengthLimit = args.l
    allelesTypes = str2list(args.a)
    tSoftwareList = str2list(args.tl)
    bSoftwareList = str2list(args.bl)
    boundary = args.d
    binSize = args.bs
    cutoff = args.c
    method = args.m
    
    time0 = time()
    time1 = time()
    print('reading epitope files...')
    dAll = bEpiLength(sourceFn)
    dAll = {x: dAll[x] for x in dAll if dAll[x]['length'] < epiLengthLimit}
    time2 = time()
    tSeg = round(time2 - time1, 3)
    tTot = round(time2 - time0, 3)
    print(f'finished reading epitope files, {len(dAll)} epitopes total, {tSeg} s, {tTot} s total')
    time1 = time()
    print('reading t and b cell prediction data...')
    totalBins = int(math.ceil(boundary/binSize)) * 2
    if method == 2:
        totalBins += 1
    refHeader = ['Epitope'] + [x for x in bSoftwareList]
    for software in tSoftwareList:
        for t in allelesTypes:
            refHeader += [f'{t}_{software}_{binNumber}' for binNumber in range(1, (totalBins + 1))]
    if not os.path.isdir(outDir):
        os.mkdir(outDir)
    inputList = [tuple([t_data_dir, b_data_dir, 
                 x, allelesTypes, tSoftwareList, bSoftwareList, boundary, 
                 dAll[x]['length'], binSize, cutoff, method, refHeader]) 
                 for x in dAll]
    with Pool() as pool:
        l = pool.starmap(readEpi, inputList)
    results = [refHeader] + [x for x in l if type(x) != int]
    results = [[str(x) for x in y] for y in results]
    percentSuccess = round((len(results) - 1)/len(l)*100, 2)
    time2 = time()
    tSeg = round(time2 - time1, 3)
    tTot = round(time2 - time0, 3)
    print(f'finished reading, {(len(results) - 1)} out of {len(l)} found ({percentSuccess}%), {tSeg} s, {tTot} s total')
    time1 = time()
    origFn = os.path.split(sourceFn)[-1]
    fn = os.path.join(outDir, f'{origFn}.csv')
    with open(fn, 'w') as f:
        s = '\n'.join([','.join(x) for x in results])
        f.write(s)