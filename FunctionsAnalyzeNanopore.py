import gzip
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import math
import csv
import multiprocessing as multiprocessing
from Bio import SeqIO
from Bio import AlignIO
import os
import re
from scipy import stats
from Bio import pairwise2
from Bio.pairwise2 import format_alignment 
from Bio.Seq import Seq 
from multiprocessing import Pool
sns.set_style("whitegrid")

def ProcessFile(file,GFPseqF,Dir,DirOut,Dictionaries,Parameters,AlignmentParameters,IndexPositions):
    DictTranslate,DictTranslateREV,DictComp,DictMultiplexing = Dictionaries
    MaxLength,LibraryLength,QualityThres,AlignmentThres,PrimerFWD,PrimerREV = Parameters
    ScoreMatch,ScoreMisMatch,GapStart,GapCont = AlignmentParameters
    # filter for reads that are at least 80% the lenght of the reference sequence
    if len(GFPseqF) > 2:
        MinLength = round(0.9*len(GFPseqF))
    else:
        MinLength = round(0.9*len(GFPseqF[0]))
    print(MinLength)
#     if re.match(r"K.*R1_001\.fastq\.gz$", file):
#         GFPseq = str(Seq(GFPseqR).reversecomplement())
#     else:
#         GFPseq = GFPseqF.upper()
    with gzip.open(Dir+file, "rt", encoding="utf8") as handle:
        i = 0
        for record in SeqIO.parse(handle, "fastq"):
            i+=1
        imax = i    
#         imax = 100
    Matrix = np.zeros((imax+1,MaxLength+1),int)
    MatrixQ = np.zeros((imax+1,MaxLength+1),int)
    MatrixO = np.zeros((imax+1,MaxLength+1),int)
    MatrixQO = np.zeros((imax+1,MaxLength+1),int)
    MatrixA = np.zeros((imax+1,LibraryLength),int)
    MatrixIn = np.zeros((imax+1,LibraryLength),int)
    Score = np.zeros((imax+1,1),int)
    Indexes = ['NaN']*(imax+1) 
    Experiments = ['NaN']*(imax+1)   
    if IndexPositions[0] == 0:
            Split = file.split('_merged.fastq.gz')
            File = Split[0]
            for key in DictMultiplexing.keys():
#                 print(key)
                if DictMultiplexing[key]['Genotype'] == File:
                    NewIndex = DictMultiplexing[key]['Sample']
                    NewExp = DictMultiplexing[key]['Experiment']

    with gzip.open(Dir+file, "rt", encoding="utf8") as handle:
        i = 0
        if len(GFPseqF) > 2:
            GFPseqT = GFPseqF.upper().translate(DictTranslate)
        else:
            GFPseqT = GFPseqF[0].upper().translate(DictTranslate)
        Matrix[i,0:len(GFPseqT)] = np.array([list(map(int, GFPseqT))])
        MatrixA[i,0:len(GFPseqT)] = np.array([list(map(int, GFPseqT))])
        for record in SeqIO.parse(handle, "fastq"):
            i+=1
            if i < imax:
                SeqO = str(record.seq)
                Quality = record.letter_annotations["phred_quality"]
                SeqR = SeqO.translate(DictTranslate)
                if len(Quality) > MaxLength:
                    Quality = Quality[0:MaxLength]
                    SeqR = SeqR[0:MaxLength]
                MatrixQO[i,0:len(Quality)] = np.array([list(map(int, Quality))])
                MatrixO[i,0:len(SeqR)] = np.array([list(map(int, SeqR))])
#                 print((sum(Quality)/len(Quality)))
                if (sum(Quality)/len(Quality) > QualityThres) & (len(Quality) >= MinLength):
                    MatrixQ[i,0:len(Quality)] = np.array([list(map(int, Quality))])
                    Matrix[i,0:len(SeqR)] = np.array([list(map(int, SeqR))])
                    Reverse = 0
                    if len(GFPseqF) > 2:
                        seq1 = Seq(GFPseqF.upper()) 
                        seq2 = Seq(SeqO)
                        seq3 = Seq(SeqO).reverse_complement()
                        # RUN ALIGMENT IN FWD AND REV AND PICK BEST ONE
                        alignmentsFWD = pairwise2.align.localms(seq1, seq2, ScoreMatch,ScoreMisMatch,GapStart,GapCont,one_alignment_only = True) 
                        alignmentsREV = pairwise2.align.localms(seq1, seq3, ScoreMatch,ScoreMisMatch,GapStart,GapCont,one_alignment_only = True)
                        if alignmentsFWD[0].score > alignmentsREV[0].score:
                            alignments = alignmentsFWD
                        else:
                            alignments = alignmentsREV
                            seq2 = seq3
                            Reverse = 0
                    else:
                        seq1 = Seq(GFPseqF[0]) 
                        seq1R = Seq(GFPseqF[1])
                        seq2 = Seq(SeqO)
                        seq3 = Seq(SeqO).reverse_complement()
                        # RUN ALIGMENT IN FWD AND REV AND PICK BEST ONE
                        alignmentsFWD1 = pairwise2.align.localms(seq1, seq2, ScoreMatch,ScoreMisMatch,GapStart,GapCont,one_alignment_only = True) 
                        alignmentsREV1 = pairwise2.align.localms(seq1, seq3, ScoreMatch,ScoreMisMatch,GapStart,GapCont,one_alignment_only = True)
                        alignmentsFWD2 = pairwise2.align.localms(seq1R, seq2, ScoreMatch,ScoreMisMatch,GapStart,GapCont,one_alignment_only = True) 
                        alignmentsREV2 = pairwise2.align.localms(seq1R, seq3, ScoreMatch,ScoreMisMatch,GapStart,GapCont,one_alignment_only = True)
#                         print([alignmentsFWD1[0].score,alignmentsREV1[0].score,alignmentsFWD2[0].score,alignmentsREV2[0].score])
                        
                        if alignmentsFWD1[0].score > alignmentsREV1[0].score:
                            alignments1 = alignmentsFWD1
                        else:
                            alignments1 = alignmentsREV1
                            seq2 = seq3
                        if alignmentsFWD2[0].score > alignmentsREV2[0].score:
                            alignments2 = alignmentsFWD2
                        else:
                            alignments2 = alignmentsREV2
                            seq2 = seq3
                        if alignments1[0].score > alignments2[0].score:
                            alignments = alignments1
                            Reverse = 1
                        else:
                            alignments = alignments2
                            Reverse = 2
#                         print(seq2)    
#                         print(Reverse)

                    for alignment in alignments: 
                        if alignment.score > AlignmentThres:
                            #print(format_alignment(*alignment))
                            Score[i] = alignment.score
                            SeqA = list(alignment.seqA)
                            SeqBO = list(alignment.seqB)
                            
                            for j in range(0,len(SeqA)):
                                if SeqA[j] == "-": 
                                    SeqBO[j] = "I"
                            SeqBO = ''.join(SeqBO)
                            #save insertions in a new matrix before deleting with DictTranslate
                            SeqBI = SeqBO
                            SeqB = SeqBO.translate(DictTranslate)
                            SeqBnoI = SeqB.translate(DictTranslateREV)
                            Index = SeqBnoI[IndexPositions[0]:IndexPositions[1]]
                          
                            SeqBIn = ['0']*len(SeqB)
                            #reduce multiple III to I to save individual ocurrences of insertion
                            while 'II' in SeqBI:
                                 SeqBI = SeqBI.replace("II", "I")
                            #annotate where that indel happened. save coordinate and delete to annotate next
                            while 'I' in SeqBI:
                                Pos = SeqBI.find('I')
                                SeqBIn[Pos-1] = '1'
                                SeqBI = SeqBI.replace("I", "",1)
                            if Reverse == 2:
                                #SeqBOREV = SeqBO[0:len(PrimerFWD)]+SeqBO[-len(PrimerREV)-1:len(PrimerFWD)+1:-1].translate(DictComp)+SeqBO[-len(PrimerREV):-1]
                                SeqBOREV = SeqBnoI[0:len(PrimerFWD)]+SeqBnoI[-len(PrimerREV)-1:len(PrimerFWD)+1:-1].translate(DictComp)+SeqBnoI[-len(PrimerREV):-1:1]
#                                 SeqBOREV = SeqBO.translate(DictComp)
                                SeqB = SeqBOREV.translate(DictTranslate)
                                SeqBIn = SeqBIn[0:len(PrimerFWD)]+SeqBIn[-len(PrimerREV)-1:len(PrimerFWD)+1:-1]+SeqBIn[-len(PrimerREV):-1:1]
#                                 print('here')
                            else:
                                SeqBIn = ''.join(SeqBIn)
         
                            MatrixA[i,0:len(SeqB)] = np.array([list(map(int, SeqB))])
                            MatrixIn[i,0:len(SeqBIn)] = np.array([list(map(int, SeqBIn))])
                              
                            #print(Index)
                            found = False
                            if IndexPositions[0] == 0:
                                Indexes[i] = NewIndex
                                Experiments[i] = NewExp
                            else:
                                for key in DictMultiplexing.keys():
                                    if DictMultiplexing[key]['Index'] == Index:
                                        Indexes[i] = DictMultiplexing[key]['Sample']
                                        Experiments[i] = DictMultiplexing[key]['Experiment']
                                        #Descriptions[i] = DictMultiplexing[key]['Experiment']
                                        found = True
                                if not found:
                                    Indexes[i] = 'Unassigned'
                                    Experiments[i] = 'None'
                                


        #there is a bug in np.delete than when there is nothing to delete it reorders input
        #Sum = np.sum(Matrix,axis=1)
        #Matrix = np.delete(Matrix, (Sum == 0), axis=0)
#         hm = sns.clustermap(Matrix,cmap="viridis",vmin=0, vmax=5,col_cluster=False,row_cluster=False)
#         hm.savefig(DirOut+file+"_heatmapFilt.png")
#         hm.fig.clf()
        
#         hmO = sns.clustermap(MatrixO,cmap="viridis",vmin=0, vmax=5,col_cluster=False,row_cluster=False)
#         hmO.savefig(DirOut+file+"_heatmapO.png")
#         hmO.fig.clf()
        
        #Sum = np.sum(MatrixQ,axis=1)
        #MatrixQ = np.delete(MatrixQ, (Sum == 0), axis=0)
#         hmQ = sns.clustermap(MatrixQ,cmap="viridis",col_cluster=False,row_cluster=False)
#         hmQ.savefig(DirOut+file+"_heatmapQFilt.png")
#         hmQ.fig.clf()
        
#         hmQO = sns.clustermap(MatrixQO,cmap="viridis",col_cluster=False,row_cluster=False)
#         hmQO.savefig(DirOut+file+"_heatmapQO.png")
#         hmQO.fig.clf()
                    
        if len(GFPseqF) > 2:
            Score[0] = len(GFPseqF)*ScoreMatch
        else:
            Score[0] = len(GFPseqF[0])*ScoreMatch
        Indexes[0] = 'REF'
        Experiments[0] = 'REF'
        Score = Score[:,0]
  
        Order = np.argsort(Score)
        Score = Score[::-1]
        Order = Order[::-1]
        Indexes = [Indexes[i] for i in Order]
        Experiments = [Experiments[i] for i in Order]
        MatrixA = MatrixA[Order,:]
        MatrixIn = MatrixIn[Order,:]
        Sum = np.sum(MatrixA,axis=1)
        MatrixA = np.delete(MatrixA, (Sum == 0), axis=0)
        MatrixIn = np.delete(MatrixIn, (Sum == 0), axis=0)
        Indexes = [x for x in Indexes if x != 'NaN']
        Experiments = [x for x in Experiments if x != 'NaN']  

        Score = np.delete(Score, (Sum == 0), axis=0)
#         ScoreDF = pd.DataFrame(Score)
#         row_colors = ScoreDF[0].map(lut).to_numpy()
#         hmA = sns.clustermap(MatrixA,cmap="viridis",vmin=0,vmax=5,col_cluster=False,row_cluster=False,row_colors=row_colors)
        DF = pd.DataFrame(MatrixA)
        DF['Score'] = Score
        DF['Sample'] = Indexes
        DF['Experiment'] = Experiments
        

        DFIn = pd.DataFrame(MatrixIn)
        DFIn['Score'] = Score
        DFIn['Sample'] = Indexes
        DFIn['Experiment'] = Experiments
#        Ceil = math.ceil(LibraryLength/6)
#        Replace = [0,1,2,3,4,5]*Ceil
#        DF.iloc[-1,0:LibraryLength] = Replace[0:LibraryLength]
        lut1 = dict(zip(set(list(range(1,len(Order),1))), sns.color_palette("viridis",max(Score))))
        lut2 = dict(zip(set(DF['Sample']), sns.color_palette("viridis",len(DF['Sample'].unique()))))
        lut3 = dict(zip(set(DF['Experiment']), sns.color_palette("inferno",len(DF['Experiment'].unique()))))
        row_colors1 = DF['Score'].map(lut1)
        row_colors2 = DF['Sample'].map(lut2)
        row_colors3 = DF['Experiment'].map(lut3)
        row_colors = pd.concat([row_colors1,row_colors2,row_colors3],axis=1)
        Cbar_kws={'ticks':[],'label':'N      A      C      T      G      D','boundaries':[0,1,2,3,4,5,6]}
        hmA = sns.clustermap(DF.iloc[:,0:LibraryLength],cmap="viridis",vmin=0, vmax=5,figsize=(8.15, 35),cbar_pos=(0.1, 0.5, 0.05, 0.12),col_cluster=False,row_cluster=False,row_colors=row_colors,cbar_kws=Cbar_kws)
        hmA.ax_heatmap.axhline(1, color='white', lw=1)
        hmA.ax_heatmap.text(-18, len(Score)*0.99,str(max(Score)))
        hmA.ax_heatmap.text(-18, len(Score)*0.02,str(min(Score)),color='k')
        hmA.savefig(DirOut+file+"_heatmapA_SortedScore.png")
        hmA.fig.clf()
        plt.close()
        
        DF = DF.sort_values('Sample', ascending=True)
        DFIn = DFIn.sort_values('Sample', ascending=True)
        lut1 = dict(zip(set(list(range(1,len(Order),1))), sns.color_palette("viridis",max(Score))))
        lut2 = dict(zip(set(DF['Sample']), sns.color_palette("viridis",len(DF['Sample'].unique()))))
        lut3 = dict(zip(set(DF['Experiment']), sns.color_palette("inferno",len(DF['Experiment'].unique()))))
        row_colors1 = DF['Score'].map(lut1)
        row_colors2 = DF['Sample'].map(lut2)
        row_colors3 = DF['Experiment'].map(lut3)
        row_colors = pd.concat([row_colors1,row_colors2,row_colors3],axis=1)
        Cbar_kws={'ticks':[],'label':'N      A      C      T      G      D','boundaries':[0,1,2,3,4,5,6]}
        hmA = sns.clustermap(DF.iloc[:,0:LibraryLength],cmap="viridis",vmin=0, vmax=5,figsize=(8.15, 35),cbar_pos=(0.1, 0.5, 0.05, 0.12),col_cluster=False,row_cluster=False,row_colors=row_colors,cbar_kws=Cbar_kws)
        hmA.ax_heatmap.axhline(1, color='white', lw=1)
        hmA.ax_heatmap.text(-18, len(Score)*0.99,str(max(Score)))
        hmA.ax_heatmap.text(-18, len(Score)*0.02,str(min(Score)),color='k')
        hmA.savefig(DirOut+file+"_heatmapA_SortedIndex.png")
        hmA.fig.clf()
        plt.close()
        
        sns.set(rc={'figure.figsize':(10,10)})
    
        hist = sns.histplot(Score)
        hist.set(title=file+" distribution of aligment scores to GFP")
        fig = hist.get_figure()
        fig.savefig(DirOut+file+"_ScoreHist.pdf")
        fig.clf()
        plt.close()
        
        #MatrixA[MatrixA==0]=np.nan
        #MatrixA[MatrixA==5]=np.nan
#         m = stats.mode(MatrixA)
#         Mismatch = []
#         for i in range(0,len(MatrixA[0])):
#             k = 0
#             for j in MatrixA[:,i]:
#                 if j!=m[0][0][i]: 
#                     k+=1
#             Mismatch.append(k/len(MatrixA[:,i]))
#         plot = sns.lineplot(x=range(0,len(MatrixA[0])),y=Mismatch)
#         plot.set(title=file)
#         fig = plot.get_figure()
#         fig.savefig(DirOut+file+"_mismatched.pdf")
#         fig.clf()
#         plt.close()
        
        DF.to_csv(DirOut+file+"_DF.csv")
        DFIn.to_csv(DirOut+file+"_DFInsertions.csv")
        return
                                     
        #separate by sample and calculate mutation profile for each
#         MatrixMM = np.zeros((len(DF['Sample'].unique()),LibraryLength),float)
#         MatrixIndels = np.zeros((len(DF['Sample'].unique()),LibraryLength),float)
#         REF = DF[DF['Sample'] == 'REF']
#         REFgRNAs = REF.copy()
#         REFgRNAs.loc[:,111:131] = 0
#         REFgRNAs.loc[:,162:182] = 0
#         REFgRNAs.loc[:,177:200] = 0
#         REFgRNAs.loc[:,207:228] = 0
#         DF[DF['Sample'] == 'REF'] = REFgRNAs
#         samplei=0
#         fig, axes = plt.subplots(5, 6, sharex=True, sharey=True, figsize=(20,10))
#         for sample in DF['Sample'].unique():
#             SubDF = DF[DF['Sample'] == sample]
#             SubDFIn = DFIn[DFIn['Sample'] == sample]
#             Mismatch = []
#             Indels = []
#             for i in range(0,LibraryLength):
#                 mm = 0
#                 indel = 0
#                 for j in range(0,len(SubDF.iloc[:,i])):
#                     if (SubDF.iloc[j,i]!=REF[i].item()) & (SubDF.iloc[j,i]!=5): 
#                         mm+=1
#                     if (SubDF.iloc[j,i]==5) | (SubDFIn.iloc[j,i]==1):
#                         indel+=1
#                 Mismatch.append(mm/len(SubDF.iloc[:,i]))
#                 Indels.append(indel/len(SubDF.iloc[:,i]))
#             MatrixMM[samplei,:] = Mismatch
#             MatrixIndels[samplei,:] = Indels
#             samplei+=1
#             plot = sns.lineplot(ax=axes[math.ceil(samplei/6)-1,samplei-((math.ceil(samplei/6))-1)*6-1],x=range(0,len(Mismatch)),y=Mismatch)
#             plot = sns.lineplot(ax=axes[math.ceil(samplei/6)-1,samplei-((math.ceil(samplei/6))-1)*6-1],x=range(0,len(Indels)),y=Indels)
#             plot.set(title=sample)
# #             fig = plot.get_figure()
#         fig.savefig(DirOut+file+"_mismatched_samples.pdf")
#         fig.clf()
#         plt.close()
        
#         hmM = sns.clustermap(MatrixMM,col_cluster=False,row_cluster=False,vmin=0,vmax = 0.4)
#         hmM.savefig(DirOut+file+"_heatmapMismatched.png")
#         hmM.fig.clf()
#         plt.close()
def CalculateMismatchIndels(DFall,DFInall,sample,LibraryLength,REF):
    SubDF = DFall[DFall['Sample'] == sample]
    SubDFIn = DFInall[DFInall['Sample'] == sample]
    Neach = len(SubDF.iloc[:,1])
    SubDF = np.array(SubDF.iloc[:,1:LibraryLength+1])
    SubDFIn = np.array(SubDFIn.iloc[:,1:LibraryLength+1])
    REFlong = np.array(np.tile(REF.iloc[0,1:LibraryLength+1],(Neach,1)))
    # element-wise compare read to refence, discount mismatched counted by deletions
    Comp = np.not_equal(SubDF, REFlong)
    Mismatch = np.sum(Comp ^ (SubDF==5),axis=0)
    Indels = np.sum(SubDF==5+SubDFIn,axis=0)
    AC = np.sum(Comp & (REFlong==1) & (SubDF==2),axis=0)
    AT = np.sum(Comp & (REFlong==1) & (SubDF==3),axis=0)
    AG = np.sum(Comp & (REFlong==1) & (SubDF==4),axis=0)
    CA = np.sum(Comp & (REFlong==2) & (SubDF==1),axis=0)
    CT = np.sum(Comp & (REFlong==2) & (SubDF==3),axis=0)
    CG = np.sum(Comp & (REFlong==2) & (SubDF==4),axis=0)
    TA = np.sum(Comp & (REFlong==3) & (SubDF==1),axis=0)
    TC = np.sum(Comp & (REFlong==3) & (SubDF==2),axis=0)
    TG = np.sum(Comp & (REFlong==3) & (SubDF==4),axis=0)
    GA = np.sum(Comp & (REFlong==4) & (SubDF==1),axis=0)
    GC = np.sum(Comp & (REFlong==4) & (SubDF==2),axis=0)
    GT = np.sum(Comp & (REFlong==4) & (SubDF==3),axis=0)
    DFMutations = pd.DataFrame([AC,AT,AG,CA,CT,CG,TA,TC,TG,GA,GC,GT])
    DFMutations = DFMutations.T
    DFMutations.columns = ['AC','AT','AG','CA','CT','CG','TA','TC','TG','GA','GC','GT']
    AllDiff = (Mismatch+Indels)/Neach
    Mismatch = Mismatch/Neach
    Indels = Indels/Neach
    return(Mismatch,Indels,AllDiff,DFMutations,Neach,sample)

def ConcatenateDataFrame(FileSet,DirOut):
#     DFall = pd.DataFrame()
#     for f in FileSet:
#         DFeach = pd.read_csv(DirOut+f)
#         DFall = pd.concat([DFall, DFeach])
#     DFall = pd.concat((pd.read_csv(DirOut+) for f in FileSet), ignore_index=True)
#     DFall = pd.DataFrame()
# sooooo much faster to make a list (1min for 50 %data )and concatenate all at once than append or concatenate  v
    List = [0]*len(FileSet)
    i = 0
    for f in FileSet:
        #print(f)
        List[i] = pd.read_csv(DirOut+f)
        #List[i] = pd.read_csv(DirOut+f, encoding='iso-8859-1')
        i+=1
    DFall = pd.concat(List)
    return(DFall)

def CountSelectedMismatches(Samples,sample,DFall,REF,LibraryLength,MutationsSelected,GFPseqF,PrimerFWD,PrimerREV):
    Which = Samples.index(sample)
    SubDF = DFall[DFall['Sample'] == sample]
    Neach = len(SubDF.iloc[:,1])
    SubDF = np.array(SubDF.iloc[:,1:LibraryLength+1])
    REFlong = np.array(np.tile(REF.iloc[0,1:LibraryLength+1],(Neach,1)))
    # element-wise compare read to refence, discount mismatched counted by deletions
    Comp = np.not_equal(SubDF, REFlong)
    CT = Comp & (REFlong==2) & (SubDF==3)
    GA = Comp & (REFlong==4) & (SubDF==1) 
#     MismatchSelected = MutationsSelected & (Comp ^ (SubDF==5))
    MismatchSelected = MutationsSelected & (CT | GA)
    MismatchSelected = MismatchSelected[:,len(PrimerFWD)+15:len(GFPseqF)-len(PrimerREV)-75]
    NumberMutations = np.sum(MismatchSelected,axis=1)
    return(NumberMutations,Neach)

def CountAllMismatches(Samples,sample,DFall,REF,LibraryLength,GFPseqF,PrimerFWD,PrimerREV):
    Which = Samples.index(sample)
    SubDF = DFall[DFall['Sample'] == sample]
    Neach = len(SubDF.iloc[:,1])
    SubDF = np.array(SubDF.iloc[:,1:LibraryLength+1])
    REFlong = np.array(np.tile(REF.iloc[0,1:LibraryLength+1],(Neach,1)))
    # element-wise compare read to refence, discount mismatched counted by deletions
    Comp = np.not_equal(SubDF, REFlong)
#     MismatchSelected = MutationsSelected & Comp
    MismatchSelected = Comp[:,len(PrimerFWD)+15:len(GFPseqF)-len(PrimerREV)-75]
    NumberMutations = np.sum(MismatchSelected,axis=1)
    return(NumberMutations,Neach)
