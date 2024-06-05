#!/usr/bin/env python -i

print("Welcome to HQAlign_DT Version June 4 2024")

## 10/14/18  Fixed a bunch of bugs related to fenceposting of output data for tabular files, and for including base+length counts for antisense as well as sense strand.  These are essential fixes
## 10/14/18  Also fixed problems in importing excel files, so HQAlign should now work well with XLS or XLSX files.  Note that a major problem was from Microsoft switching boolean values from "True"/"False" to '1'/'1' in Excel.  XLSX import is still experimental
## 06/04/24  Fixed file open to avoid deprecated 'rU' option

from shutil import copy2
from array import *
from sys import *
from os.path import *
import datetime, gzip
from time import gmtime, strftime, asctime, localtime, time
from re import split as resplit
from subprocess import check_output,check_call
import platform
import ctypes
from os import chdir,getcwd,stat,walk
import string

## the following section has a number of interface enhancements that can be modified as described for "safe" compiling"
## for safe compiling comment out the next fourteen lines
## Note a bunch of this is so that HQAlign works in Python 3.xx.  But 3.xx is currently somehwat not optimal in that the default string type for 3.xx is unicode, which wastes both memory and clock cycles.
try:
    import cPickle
except:
    import pickle as cPickle
try:
    from Tkinter import *
    from tkFileDialog import asksaveasfilename,askopenfilename
except:
    from tkinter import *
    from tkinter.filedialog import askopenfilename,asksaveasfilename
try:
    import statvfs
except:
    pass
    
def GoodChars(s0):
    '''Filter out everything but 0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~ \t\n\r\x0b\x0c'''
    return ''.join(list(filter(lambda x: x in string.printable, s0)))

def UserGetFileNameRead1(msg1):
    FilePath1=''

## for safer compiling comment out the next four lines
    root = Tk()
    root.withdraw()
    FilePath1=askopenfilename(title=msg1,initialdir=getcwd(),parent=root)
    root.quit()

    if not(isfile(FilePath1)):
        print("No file selected: Using Default File HQSheet.txt")
        FilePath1='HQSheet.txt'
    return FilePath1

def UserGetFileNameWrite1(msg1):
    FilePath1=''

## for safer compiling comment out the next four lines
    root = Tk()
    root.withdraw()
    FilePath1=asksaveasfilename(title=msg1,initialdir=getcwd(),parent=root)
    root.quit()

    if not(isfile(FilePath1)):
        print("No file selected: Stopping HQAlign")
        exit()
    return FilePath1
    
def ImportXLS1(BatchFilePath):
    Setup11=''
    ## comment out the next 11 lines for "safe" compiling (also this won't let you use Excel files for input)"
    try:
        import xlrd
    except:
        print("Excel File Converter 'xlrd' not installed on this system, try saving spreadsheet as a 'txt' or 'csv' file and rerunning")
        exit()
    wb1 = xlrd.open_workbook(BatchFilePath)
    sh1 = wb1.sheet_by_index(0)
    for rownum in range(sh1.nrows):
        for colnum in range(sh1.ncols):
            Setup11+=str(sh1.cell(rownum,colnum)).split(':')[-1].lstrip('u').strip("'")+'\t'
        Setup11+='\r'
    if Setup11=='':
        print("Excel File Converter 'xlrd' not installed on this system, try saving spreadsheet as a 'txt' or 'csv' file and rerunning")
        exit()
    return Setup11


def arguments1():
    ## comment out the following seven lines for safe compile 
    try:
        args1=[]
        for arg11 in argv:
            args1.append(arg11)
        return args1
    except:
        return ''

    return ''

## from Paul Wicks, Frankovskyi Bogdan, TriPhoenix, JF Sebastian, http://stackoverflow.com/questions/51658/cross-platform-space-remaining-on-volume-using-python

def get_free_space(folder):
    """ Return folder/drive free space (in bytes)"""
    ## comment out the following nine lines for safe compile 
    try:
        if platform.system() == 'Windows':
            free_bytes = ctypes.c_ulonglong(0)
            ctypes.windll.kernel32.GetDiskFreeSpaceExW(ctypes.c_wchar_p(folder), None, None, ctypes.pointer(free_bytes))
            return free_bytes.value
        else:
            return statvfs(folder).f_bfree*statvfs(folder).f_frsize
    except:
        return 0
    return 0

gzipimported1=False ## Set to true for the first time gzip is imported

fastqDumpPaths1=[getcwd(),'~/Downloads','C:\Downloads']
fastqDumpPath1=join(getcwd(),'fastq-dump')
def sraTranslate(sraPath,fastqDumpPath):
    '''translates sra file sraPath into .fasta format using fastq-dump program at fastqDumpPath.  Returns fastApath,fastqDumpPath)'''
    Qprint(time0,LogFileName,'Trying .sra to .fasta translation using fastq-dump, with default location '+fastqDumpPath) 
    translated=False
    try:
        j=check_call([fastqDumpPath,'--fasta',sraPath])
        if j==0:
            translated=True
    except:
        pass
    if not(translated):
        for r1 in fastqDumpPaths1:
            if not(isdir(r1)):
                continue
            T0=time()
            for curdir, dirs, files in walk(r1):
                T1=time()-T0
                if T1>15:
                    break
                for f in files:
                    if f.startswith('fastq-dump'):
                        fastqDumpPath=os.path.join(r1,curdir,f)
                        j=check_call([fastqDumpPath,'--fasta',sraPath])
                        if j==0:
                            translated=True
                    if translated:
                        break
                if translated:
                    break
            if translated:
                break
    if not(translated):
        Qprint(time0,LogFileName,'fastqDump not found at the locations that were checked within time limit.. asking user to locate '+fastqDumpPath) 

        fastqDumpPath=UserGetFileNameRead1('sra to fasta conversion failed, please locate fastq-dump')
        j=check_call([fastqDumpPath,'--fasta',sraPath])
        if j==0:
            translated=True
    if not(translated):
        Qprint(time0,LogFileName,'Unable to translated .sra file, try putting sratools into the directory with the script')
    return (sraPath[:-3]+'fasta',fastqDumpPath)        
    

#################################
## Start of user-set variables ##
#################################
## Debug options in this search can be located by doing a text search for "Debug"
## Option 0: Run the script as a __main__ program (for an interactive debug after running the program).  Uncomment the line ## if true below
## Option 1: Cprofile, comment the line "##if true..." and uncomment ## import Cprofile and Cprofile.run('HQMain()'))
## Option 2: is a dump of local variables at the end of the text

def Log(LogFileNameA,sLog):
    LogTime=strftime("%d %b %Y %H:%M:%S", localtime())
    LogFile=open(LogFileNameA,'a')
    sss=sLog.split('\r')
    if len(sss)==1:
        LogFile.write('!Log ' + LogTime + '   ' + sLog.strip('\r')+'\r')
    else:
        LogFile.write('!Log ' + LogTime)
        for sLog1 in sss:
            LogFile.write('!Log ' + LogTime + '   ' + sLog1.strip('\r')+'\r')
    LogFile.close()
def Qprint(timeA0,LogFileNameQ,qprint1):
    curtime=time()-timeA0
    print(qprint1+' | t=%.3f sec.'%curtime)
    Log(LogFileNameQ,qprint1+' | t=%.3f sec.'%curtime)

def SuperParse(bbs1,Dlocals):
    '''parse a string bbs1 into a series of safe assignment statements
    some conventions
    #-to-end of line or delimiter is comment
    commas, tabbs, and spaces outside of quotes separate elements
    leading and trailing spaces are deleted
    equals outside of a quote is king'''
    bbs2=bbs1+''
    bbs3=[] ## an array of entries, each of which is a 2-tuple with a value and a type, t=(quoted) text, u=(unquoted)
    bbs4=''
    bt1=''
    BatchList=[]
    while bbs2:
        bbs2=bbs2.strip()
        if bbs2.startswith('#StartBatchTable'):
            if bbs2.find('#EndBatchTable')>=0:
                bt1=bbs2[:bbs2.find('#EndBatchTable')]
                bbs2=bbs2[bbs2.find('#EndBatchTable'):]
            else:
                bt1=bbs2+''
                bbs2=''
                break
        bbc1=bbs2[0]
        if not(bbs2):
            break
        if bbc1=="'":
            if bbs2[1:].split("'")[0]:
                bba1=bbs2[1:].split("'")[0].strip()
                if bba1 and bba1[0]!='#':
                    bbs3.append((bba1,'t'))
            bbs2=bbs2[1:].split("'",1)[-1]
            continue
        if bbc1=='"':
            if bbs2[1:].split('"')[0]:
                bba1=bbs2[1:].split('"')[0].strip()
                if bba1 and bba1[0]!='#':
                    bbs3.append((bba1,'t'))
            bbs2=bbs2[1:].split('"',1)[-1]
            continue
        if bbc1=='#':
            bbs5=[]
            if bbs2.find('\t')>=0:
                bbs5.append(bbs2.find('\t'))
            if bbs2.find('\r')>=0:
                bbs5.append(bbs2.find('\r'))
            if bbs2.find('\n')>=0:
                bbs5.append(bbs2.find('\n'))
            if not(bbs5):break
            bbs2=bbs2[min(bbs5):]
            bbs2=bbs2.strip()
            if not(bbs2):
                break
            continue
        bbs5=[]
        if bbs2.find('\t')>=0:
            bbs5.append(bbs2.find('\t'))
        if bbs2.find('\r')>=0:
            bbs5.append(bbs2.find('\r'))
        if bbs2.find('\n')>=0:
            bbs5.append(bbs2.find('\n'))
        if bbs2 and not(bbs5):
            bbs3.append((bbs2.strip(),'u'))
            break
        bbs3.append((bbs2[:min(bbs5)].strip(),'u'))
        bbs2=bbs2[min(bbs5):].strip()
    bbs6=[]
    bbs8=[]
    for bbs7 in bbs3:
        if bbs7[1]=='t':
            bbs6.append(bbs7[0])
        else:
            bbs9=resplit('[\t\r\n, =]',bbs7[0])
            for bbs10 in bbs9:
                bbs10=bbs10.strip().strip("'"+'"').strip()
                bbs6.append(bbs10)
    bbs11=[]
    for bbj1 in bbs6:
        bbk1=bbj1.strip('=').strip()
        if bbk1 in Dlocals:
            bbt1=type(Dlocals[bbk1])
            bbs11.append([bbk1,''])
            continue
        if not(bbs11):
            continue
        if bbj1=='' or bbs11[-1][1]!='':continue
        if bbt1==bool:
            bbj1=bbj1[0].upper()+bbj1[1:].lower()
            if 'f' in bbj1.lower() or '0' in bbj1:  #oops, python does bool('False')=True, need to correct for this.  So any word with an 'f' (upper or lower case) will be interpreted as 'false' for a boolean variable.  So keep your language clean!
                bbj1=''
        try:
            bbs11[-1][1]=bbt1(bbj1)
        except:
            try:
                bbs11[-1][1]=bbt1(float(bbj1))
            except:
                bbs11[-1][1]=bbj1
    if bt1:
        bt1s=bt1.splitlines()
        HeaderRowNow=False
        HeaderRow=[]
        currow1=1
        for btl1 in bt1s:
            blli1=btl1.split('\t')
            for blli11 in blli1:
                blli11=blli11.strip('=').strip()
                if blli11 in Dlocals:
                    HeaderRowNow=True
            if HeaderRowNow:
                HeaderRow=[blli11.strip('=').strip() for blli11 in blli1]
                break
            currow1+=1
        if not(HeaderRowNow):
            currow1=1
        for btl1 in bt1s[currow1:]:
            BatchList.append([])
            blli1=btl1.split('\t')
            blci1=0
            for blli11 in blli1:
                blli11=blli11.strip("'"+'"\t\r\n =')
                if "=" in blli11:
                    BatchList[-1].append([blli11.split('=')[0].strip("'"+'"\t\r\n '),blli11.split('=')[1].strip("'"+'"\t\r\n ')])
                else:
                    if blci1<len(HeaderRow) and HeaderRow[blci1]:
                        BatchList[-1].append([HeaderRow[blci1],blli11.strip("'"+'"\t\r\n =')])
                blci1+=1
        bbs11.append(['BatchTable',BatchList])
    return bbs11



def AntiSense(s,filterminus):  
    return s.translate(filterminus)[::-1]

def SeqStrip(s,filterplus):
    return s.translate(filterplus).translate(filterplus,' ')
## the two stage translation here is to allow the first stage to convert all stray characters to ' ', then nuke all of these characters looking for a single character

def sdi1(X):
        X2=list(X.keys())
        if X2: X2.sort()
        X4='('
        for X3 in X2: X4=X4+str(X3)+':'+str(X[X3])+','
        if X4[-1]==',':X4=X4[:-1]
        X4+=')'
        return X4

def fillerchar(nf1):
    if nf1%10==0:
        return '|'
    if nf1%10==5:
        return ','
    return '.'
def fillerspace(a,b):
    fs1=''
    for indef1 in range(a,b+1):
        fs1+=fillerchar(indef1)
    return fs1
        
def Ruler(a,b,nm1):
    R1='#'
    R1+=str(a)
    R1+=fillerchar(len(R1)+1)    
    R2i=a+len(R1)
    R1add=''    
    while R2i+len(R1add)<b:
        R1+=R1add
        R2i=a+len(R1)
        Next20=(R2i//20)*20+20
        if Next20>=b:
            break
        L20=len(str(Next20))
        SpaceToFill=Next20-L20-R2i
        if SpaceToFill<1:
            R1add=fillerspace(R2i,b-1)
        else:
            R1add=fillerspace(R2i,R2i+SpaceToFill-1)+str(Next20)+fillerchar(Next20)
        if Next20%80==0 and Next20+1+len(nm1)<b:
            R1add+=nm1
    R1+=fillerspace(R2i,b-1)+'#'+str(b)
    return R1

def Transpose(infileT,outfileT,Memb):
    fT=open(infileT,'rt')
    fU=open(outfileT,'w')
    lines=1
    for L2 in fT:
        if lines==1:cols=len(L2.split('\t'))
        lines+=1
    fT.seek(0)
    CurFea=0
    lineM=Memb//lines
    NotDone=True
    ElementsInChunk=9999999999999999
    LL=['moose']*lines
    CurPos=0
    BreakIndex=0
    while NotDone:
        linex=0
        for L2 in fT:
            L2=L2.strip('\r\n')
            L3='\t'.join(L2.split('\t')[CurPos:])
            if len(L3)>lineM:
                L3=L3[:lineM]
                L3=L3.split('\t')[:-1]
            else:
                L3=L3.split('\t')
            if len(L3)<ElementsInChunk:
                ElementsInChunk=len(L3)
            LL[linex]=L3[:ElementsInChunk]
            linex+=1
        fT.seek(0)
        if ElementsInChunk==0:
            return "Insufficient memory for File Transposition"
            break
        for EIC in range(ElementsInChunk):
            for linex in range(lines-1):
                fU.write(LL[linex][EIC])
                if linex<lines-2:
                    fU.write('\t')
            if EIC<ElementsInChunk-1: fU.write('\r')
        CurPos=CurPos+ElementsInChunk
        if CurPos>=cols: break
        fU.write('\r')
        ElementsInChunk=9999999999999999
    fT.close()
    fU.close()
    return "File Transposition completed successfully"

def MainHQ():
    pass

## Debug option 1... to turn this script into a function (allowing Cprofile among other options), comment out the next line

if True:
    print("Main Loop")
    ## Character of overall Experiment
    ExptDesc='ExptDesc-NotSet'   ## Description for the entire experiment
    ReferenceFileName='ReferenceFileName-Default_Not_Set' ##'/Users/firelab08/Desktop/JuliaTryptic.txt'       ## file with reference sequences, 'ask' to ask the user

    ## Description of an individual samples that will be matched to Reference Files (the Sequencer output and what it is)
    SampleDesc='FirstSample'## Text Describing the individual sample
    ReadFileName='ReadFileName-Default_Not_Set'       ## 'ask' asks the user using a dialog box
    StartBarcode=''         ## Barcode at start of sequence (e.g. 'GCAC'). '' means no barcode.  Use a series of k "N"'s before the barcode (or just a series of k "N"'s and no barcode) to cut off the first k bases of the sequence before analysis 
    EndLinker=''            ## Linker at end of sequence ('GCACGGACCAGATTGA'). '' means no linker or not in sequence, can specify just a few bases at the beginning if these are unique in the linker, and this will require that only these match the linker
    DataUpper=True  ## True=All data is in upper case (also all nonstandard characters will be treated as "A".  Setting this to false performs an upper casing and filtering of data for AacCGgTt

    ## BatchTable provides A tab or comma-delimited list format for batch processing of samples
    MultipleDataSets=False       ## Set MultipleDataSets=True to use batch format for entry of sample id's.  Other variable settings can be any command that should be run
                            ## at the time the individual alignments for a given sample are to be executed (e.g. you can have two lines with the
                            ## same files and barcodes but with different Sample Descriptions and number of mismatches).

    BatchTable="""BatchTable-Default_Not_Set"""

    ## Character of Alignment
    MinRead=19              ## Minimal read length for analysis (Linker and Barcode lengths will be added to this value by HQAlign).  Sequences shorter than this will not be tested for alignment
    MaxRead=999             ## Maximal analyzed length (Linker and Barcode lengths will be added to this value by HQAlign).  HQAlign Ignores all bases beyond this point.
    MaxMisMatch=0           ## Maximum allowed mismatches.  In batch mode, put the largest value used here.  
    MinMatch=19             ## Minimum number of (non-linker non-barcode) matches, starting from the first base and going forward that must match, with fewer than "MM" mismatches interspersed
    ExtensionOK=True        ## True: Keep all reads meeting initial match criteria, allowing any number of downstream mismatches.  False: count mismatches through the entire read (To require perfect match from start to end of read set ExtensionOK=False and MaxMisMatch=0)
    RequireEndLinker=True  ## True- Require the specified bases of the 3' linker to be present of the tag is ignored
    Seed=12                  ## Seed length... this is the segment size that the program looks for as an initial perfect match.  Odd things can happen if Seed*MM>MinRead
    FivePrimeTrim=0
    Multiplicity1=1
    ## How to treat Multiple hits from the same tag
    MultipleMatchMode='AllChampions'

    #    Settings for MultipleMatchMode:
    #    FirstMatch: Return 1st match in reference file that meets match criteria
    #    FirstChampion: Return 1st match in reference file that (i) meets match criteria, and (ii) has the best match score
    #    AllMatches: Return all matches in reference file that meet match criteria
    #    AllChampions: Return all matches in reference file that (i) meets match criteria, and (ii) tie for best match score
    #    UniqueMatchesOnly: Return only matches in the reference file that (i) meet the match criteria and (ii) where no other sites match criteria
    #    UniqueChampionsOnly: Return only matches in the reference file that (i) have best match scores and (iii) where no other target sites have that score
    #    SuperHits_n; Return only matches in the reference file that (i) meet match criteria and (ii) where any other match to the query has at least n additional mismatches'''



    ## How to treat Multiple hits to the same sequence
    StartCollapse=False     ## StartCollapse keeps only one instance on each strand for any given base (equivalent to start occurences rather than instances)
    StartLenCollapse=False  ## StartCollapse keeps only one instance on each strand for any given base (equivalent to start+end occurences rather than instances)
    DataPreCollapsed=False   ## Data file contains read sequences that are already collapsed.  Tag names in FastA file are of the format >nnnnn-xxx where xxx is the multiplicity of each tag

    ## some Primitive ways to filter the input file (e.g. for 26G RNAs)
    TargetStart='None'           ## if defined (A,G,C, or T) then only sequence reads starting with that base will be considered
    TargetLength=0         ## if >0 (an integer) then only sequence reads of that length (following barcode and linker trimming) will be considered.  Zero to ignore

    ## Character of output data gathering- What "bins" and groups of sequences to report on 
    TotalBin=True          ## one bin that will be a summary of the total file
    GeneBin=True           ## do gene by gene summaries for each gene? Can be done as total, sense or antisense.
    BinBin=True            ## do segment by segment summaries for each gene?
    BaseBin=True           ## provide a base-by-base summary
    Granularity=30          ## how fine a granularity; 1 means that every base is recorded as a bin, 100=bins of 100

    ## Character of output formatting- Info about Reference Sequences that will be reported
    Bin_Name=True           ## Output a column with the gene name for each bin
    Bin_Size=True           ## Output a column with the size of each bin
    Bin_Start=True          ## Output a column with the start of each bin
    Bin_End=True            ## Output a column with the end of each bin
    Bin_Sequence=True       ## Output a column with the sequence of each bin
    Bin_Composition=True   ## Output 4 columns with the sequence composition

    ## Character of output formatting- What info to provide for each bin
    TabularOutput=True         ## True to provide a detailed output file (non-exclusive alternative is SimpleOutput or Pileup, but this is required for all other output)
    TabularOutputFileName='auto'      ## Output file name... auto is an automatically generated time-specific name
    TransposeFile=True      ## provides the output file in transposed format (horizontal to vertical shift)

    StartSense=True         ## Report the start of each sense match?
    StartAntiSense=True     ## Report the start of each AntiSense match?
    StartTotal=True         ## Report the start of every match?

    DyadSense=True         ## Report a hypothetical Dyad distribution, with Dyad being DyadOffset away from each start (sense reads)
    DyadAntiSense=True     ## Report a hypothetical Dyad distribution, with Dyad being DyadOffset away from each start (antisense reads)
    DyadTotal=True          ## Report a hypothetical Dyad distribution, with Dyad being DyadOffset away from each start
    DyadOffset=73            ## This is the distance between the beginning of a read and the relevant feature, e.g. 73 for a Nucleosome dyad

    EndSense=True           ## Report the end of each sense match?
    EndAntiSense=True       ## Report the end of each antisense match?
    EndTotal=True           ## Report the end of each match?

    CoverageSense=True      ## Report coverage by sense reads
    CoverageAntiSense=True  ## Report coverage by AntiSense reads
    CoverageTotal=True      ## Report coverage by total reads
    VirtualSegmentLen=147  ## This means the length of the insert will be used as the segment length for coverage (e.g. 147 for a nucleosome)

    SizeHistogramSense=True     ## Report a sense size histogram for each bin
    SizeHistogramAntiSense=True ## Report an antisense size histogram for each bin
    SizeHistogramTotal=True     ## Report an "all-hit" size histogram for each bin
    Base1SizeHistogram=True    ## Augments the size histogram to provide both size and first base (e.g. 22G versus 22A)

    BaseMatchesSense=True        ## Report a summary of matching and mismatching bases (sense hits) for each bin
    BaseMatchesAntiSense=True    ## Report a summary of matching and mismatching bases (antisense hits) for each bin
    BaseMatchesTotal=True        ## Report a summary of matching and mismatching bases (total hits) for each bin

    ReadMatchesSense=True     ## Report a summary of matching and mismatching reads for each bin (sense)-- Gets used for Pileup also
    ReadMatchesAntiSense=True ## Report a summary of matching and mismatching reads for each bin (antisense)-- Gets used for Pileup also
    ReadMatchesTotal=True   ## Report a summary of matching and mismatching reads for each bin (antisense)-- Gets used for Pileup also

    CompositionMatrix=True      ## prepare a composition list for sites around each hit.  Can be indexed (for each position) or non-indexed (one list)
    CompositionCenterEnd=True  ## False Centers on the end of each read, true on the beginning of each read.
    CompositionIndexed=True     ## True gives a composition for each position in the segment that is being analyzed.  False gives a single composition for the whole segment.
    CompositionTupleLen=2    ## >1 for longer oligonucleotide queries (tuples)
    Composition_Start=-25       ## how far upstream of a given match to start the composition window
    Composition_End=+5         ## how far downstream of a given match to extend the composition window

    SimpleOutput=True           ## True to provide a simple output file with the locations of each hit
    SimpleOutputFileName='auto' ## FileName for simple output file "ask" asks in the program, and "auto" gives a name automatically
    SimpleOutputFormat='hqall'   ## File formats for simple output
                                ## 'hqall' gives a match by match output
                                ## 'hqsum' aggregates matches for each position (only for short reference sequences below ~1MB)
                                ## 'bed' UCSC-browser compatible bed file
                                ## 'psl' Standard psl alignment output 
                                ## 'psls'.  Standard psl format with sequence appended as an extra field

    PileUp=True                 ## Generate a pileup file
    PileUpFileName='auto'       ## PileUp file name... auto is an automatically generated time-specific name
    PileSep=-1                  ## minimum separation between entries in the pileup view of the alignment.  -1 for a full pileup
    KeepFullQuery=True         ## True: Adds the full read for each match to the target (False allows the program to give just the segment that matches the target with <MaxMisMatch mismatches, i.e. the reads will effectively be truncated and ignored after that number of mismatches)

    ## How to store pre-compiled reference indexes
    StoreIndex=True         ##Keep an index based on the reference sequence?  Speeds up subsequent runs with same refereence for reference sequence sets >5MB.  Requires up to (or more than) 10x more disk space than the original file
    IndexFilePath='auto'    ##Where to store the index.  'auto' will store the index in the same directory as the original reference sequence file (original file name+'.Seed-nn_MM-nn.hqref')
    MinFreeGB=5             ##To keep HQAlign from filling your hard drive.  If the program thinks less than this amount of space will be available on the relevant storage drive after storing the index, it will avoid storing

    ## How often to report "progress on the terminal 
    ReGran=1000000          ## How frequently to report the indexing process
    ReportInterval=100000   ## How often to report progress as the program is matching queries
    StopAfterReference=0    ## For trial runs, stop after a limited number of bases in the reference sequence file; =0 for "Don't stop until finished"
    StopAfterQuery=0        ## For trial runs, stop after a limited number of queries (tags); =0 for "Don't stop until finished" 



    HardLocals={'RequireEndLinker':RequireEndLinker,
    'ReGran':ReGran,
    'ExtensionOK':ExtensionOK,
    'BaseMatchesAntiSense':BaseMatchesAntiSense,
    'Bin_Sequence':Bin_Sequence,
    'EndTotal':EndTotal,
    'SampleDesc':SampleDesc,
    'MultipleMatchMode':MultipleMatchMode,
    'ReportInterval':ReportInterval,
    'Bin_Size':Bin_Size,
    'Bin_Start':Bin_Start,
    'StopAfterQuery':StopAfterQuery,
    'SimpleOutput':SimpleOutput,
    'TabularOutputFileName':TabularOutputFileName,
    'Bin_Composition':Bin_Composition,
    'CompositionCenterEnd':CompositionCenterEnd,
    'PileUp':PileUp,
    'EndSense':EndSense,
    'CoverageAntiSense':CoverageAntiSense,
    'Composition_Start':Composition_Start,
    'MinMatch':MinMatch,
    'StopAfterReference':StopAfterReference,
    'DyadAntiSense':DyadAntiSense,
    'CompositionIndexed':CompositionIndexed,
    'Composition_End':Composition_End,
    'CoverageTotal':CoverageTotal,
    'MaxMisMatch':MaxMisMatch,
    'DyadOffset':DyadOffset,
    'CoverageSense':CoverageSense,
    'Bin_End':Bin_End,
    'PileUpFileName':PileUpFileName,
    'GeneBin':GeneBin,
    'SizeHistogramAntiSense':SizeHistogramAntiSense,
    'Granularity':Granularity,
    'CompositionTupleLen':CompositionTupleLen,
    'TotalBin':TotalBin,
    'KeepFullQuery':KeepFullQuery,
    'DyadSense':DyadSense,
    'Base1SizeHistogram':Base1SizeHistogram,
    'TargetLength':TargetLength,
    'StartBarcode':StartBarcode,
    'EndLinker':EndLinker,
    'BinBin':BinBin,
    'VirtualSegmentLen':VirtualSegmentLen,
    'PileSep':PileSep,
    'DyadTotal':DyadTotal,
    'SimpleOutputFileName':SimpleOutputFileName,
    'MultipleDataSets':MultipleDataSets,
    'SimpleOutputFormat':SimpleOutputFormat,
    'ReadMatchesTotal':ReadMatchesTotal,
    'ReadMatchesSense':ReadMatchesSense,
    'EndAntiSense':EndAntiSense,
    'MaxRead':MaxRead,
    'BaseMatchesTotal':BaseMatchesTotal,
    'ReferenceFileName':ReferenceFileName,
    'TargetStart':TargetStart,
    'ReadMatchesAntiSense':ReadMatchesAntiSense,
    'CompositionMatrix':CompositionMatrix,
    'DataUpper':DataUpper,
    'TransposeFile':TransposeFile,
    'SizeHistogramSense':SizeHistogramSense,
    'StartLenCollapse':StartLenCollapse,
    'DataPreCollapsed':DataPreCollapsed,
    'StartSense':StartSense,
    'Bin_Name':Bin_Name,
    'ExptDesc':ExptDesc,
    'BatchTable':BatchTable,
    'SizeHistogramTotal':SizeHistogramTotal,
    'StartTotal':StartTotal,
    'TabularOutput':TabularOutput,
    'Seed':Seed,
    'BaseMatchesSense':BaseMatchesSense,
    'MinRead':MinRead,
    'StartAntiSense':StartAntiSense,
    'StartCollapse':StartCollapse,
    'BaseBin':BaseBin,
    'ReadFileName':ReadFileName,
    'StoreIndex':StoreIndex,
    'IndexFilePath':IndexFilePath,
    'MinFreeGB':MinFreeGB}

    ###############################
    ## End of user-set variables ##
    ###############################
    ## start program and import from various other places

    time0=time()


    ## For eventual python 3 compatibility
    if version[0]!='2':
        xrange=range

    ASB11="AaCc-*-gGtT"  ## * serves as a barrier between different sequence entries in the full sequence lists (T will be all sense sequences catenated, U all antisense sequences
    ASB8="AaCcgGtT"
    filterplusupper=''
    for i in range(256):
        if chr(i) in ASB8:
            filterplusupper+=chr(i).upper()
        else:
            filterplusupper+=' '
    filterplus=''
    for i in range(256):
        if chr(i) in ASB8:
            filterplus+=chr(i)
        else:
            filterplus+=' '
    filterminusupper=''
    for i in range(256):
        if chr(i) in ASB11:
            filterminusupper+=ASB11[10-ASB11.find(chr(i))].upper()
        else:
            filterminusupper+=' '
    filterminus=''
    for i in range(256):
        if chr(i) in ASB11:
            filterminus+=ASB11[10-ASB11.find(chr(i))]
        else:
            filterminus+=' '

        
    BatchFilePath=''
 
    if len(arguments1())>1 and isfile(arguments1()[1]):
        BatchFilePath=arguments1()[1]


                
        

    ## Generate a version of the date and time that can be put into a file name
    now1=strftime('_Date-%m-%d-%y_Time-%H-%M-%S_',localtime())




    if not(BatchFilePath):
        BatchFilePath=UserGetFileNameRead1("File with instructions for HQ Align (e.g., 'HQSheet.xls') ")
##    if not(BatchFilePath):
##        ReadFileDir=getcwd()           ## '' default directory when ASK-ing manually for the reference file ''=directory with this program
##        ReferenceFileDir=getcwd()           ## '' default directory when ASK-ing manually for the experimental file ''=directory with this program 
##        SimpleOutputFileDir=getcwd()  ## File Directory for simple output file
##        TabularOutFileDir=getcwd()           ## Output file directory... '' is the directory with the reference file
##        PileUpFileDir=getcwd()        ## PileUp file directory... '' is the directory with the reference file
##        LogFileDir=getcwd()        ## PileUp file directory... '' is the directory with the reference file
##        LogFileName=join(LogFileDir,'HQAlignLog.txt')
    TempBatch=open(BatchFilePath,mode='r')
    BatchFileDir=dirname(abspath(TempBatch.name))
    BatchFileName=basename(TempBatch.name)
    TempBatch.close()
    chdir(BatchFileDir)
    BatchFileType=BatchFileName.rsplit('.',1)[-1]
    if BatchFileType.lower() in ('xls','xlsx'):
        Setup1=ImportXLS1(BatchFilePath)
    else:
        SetupF1=open(BatchFilePath,'rb').read()
        Setup1=''
        for i in SetupF1:
            if type(i)==str:
                i=ord(i)
            if i<128:
                Setup1+=chr(i)
        if Setup1.count('\t')<Setup1.count(','):
            Setup1=Setup1.replace(',','\t')
    Setup1=GoodChars(Setup1)
    DefaultBatchFile1=open('HQSheet.txt',mode='w')
    DefaultBatchFile1.write(Setup1)
    DefaultBatchFile1.close()
            

    ReadFileDir=BatchFileDir+''          ## '' default directory when ASK-ing manually for the reference file ''=directory with this program
    ReferenceFileDir=BatchFileDir+''           ## '' default directory when ASK-ing manually for the experimental file ''=directory with this program 
    SimpleOutputFileDir=BatchFileDir+''  ## File Directory for simple output file
    TabularOutFileDir=BatchFileDir+''          ## Output file directory.
    PileUpFileDir=BatchFileDir+''        ## PileUp file directory... '' is the directory with the reference file
    LogFileDir=BatchFileDir+''        ## PileUp file directory... '' is the directory with the reference file
    LogFileName=join(LogFileDir,'HQAlignLog.txt')
    Qprint(time0,LogFileName,'BatchFile='+TempBatch.name)
    SettingsList1=SuperParse(Setup1,HardLocals)
    for J1 in SettingsList1:
        if type(J1[1])==str:
            J1[1] = J1[1].strip().strip('"'+"'")
        if J1[0]=='RequireEndLinker':RequireEndLinker=J1[1]
        if J1[0]=='ExtensionOK':ExtensionOK=J1[1]
        if J1[0]=='ReGran':ReGran=J1[1]
        if J1[0]=='BaseMatchesAntiSense':BaseMatchesAntiSense=J1[1]
        if J1[0]=='Bin_Sequence':Bin_Sequence=J1[1]
        if J1[0]=='EndTotal':EndTotal=J1[1]
        if J1[0]=='SampleDesc':SampleDesc=J1[1]
        if J1[0]=='MultipleMatchMode':MultipleMatchMode=J1[1]
        if J1[0]=='ReportInterval':ReportInterval=J1[1]
        if J1[0]=='Bin_Size':Bin_Size=J1[1]
        if J1[0]=='Bin_Start':Bin_Start=J1[1]
        if J1[0]=='StopAfterQuery':StopAfterQuery=J1[1]
        if J1[0]=='SimpleOutput':SimpleOutput=J1[1]
        if J1[0]=='TabularOutputFileName':TabularOutputFileName=J1[1]
        if J1[0]=='Bin_Composition':Bin_Composition=J1[1]
        if J1[0]=='CompositionCenterEnd':CompositionCenterEnd=J1[1]
        if J1[0]=='PileUp':PileUp=J1[1]
        if J1[0]=='EndSense':EndSense=J1[1]
        if J1[0]=='CoverageAntiSense':CoverageAntiSense=J1[1]
        if J1[0]=='Composition_Start':Composition_Start=J1[1]
        if J1[0]=='MinMatch':MinMatch=J1[1]
        if J1[0]=='StopAfterReference':StopAfterReference=J1[1]
        if J1[0]=='DyadAntiSense':DyadAntiSense=J1[1]
        if J1[0]=='CompositionIndexed':CompositionIndexed=J1[1]
        if J1[0]=='Composition_End':Composition_End=J1[1]
        if J1[0]=='CoverageTotal':CoverageTotal=J1[1]
        if J1[0]=='MaxMisMatch':MaxMisMatch=J1[1]
        if J1[0]=='DyadOffset':DyadOffset=J1[1]
        if J1[0]=='CoverageSense':CoverageSense=J1[1]
        if J1[0]=='Bin_End':Bin_End=J1[1]
        if J1[0]=='PileUpFileName':PileUpFileName=J1[1]
        if J1[0]=='GeneBin':GeneBin=J1[1]
        if J1[0]=='SizeHistogramAntiSense':SizeHistogramAntiSense=J1[1]
        if J1[0]=='Granularity':Granularity=J1[1]
        if J1[0]=='CompositionTupleLen':CompositionTupleLen=J1[1]
        if J1[0]=='TotalBin':TotalBin=J1[1]
        if J1[0]=='KeepFullQuery':KeepFullQuery=J1[1]
        if J1[0]=='DyadSense':DyadSense=J1[1]
        if J1[0]=='Base1SizeHistogram':Base1SizeHistogram=J1[1]
        if J1[0]=='TargetLength':TargetLength=J1[1]
        if J1[0]=='StartBarcode':
            StartBarcode=J1[1]
            FivePrimeTrim=0
            while StartBarcode and StartBarcode[0].upper()=='N':
                FivePrimeTrim+=1
                StartBarcode=StartBarcode[1:]
            if StartBarcode:
                sbc=StartBarcode[0]                
        if J1[0]=='EndLinker':EndLinker=J1[1]
        if J1[0]=='BinBin':BinBin=J1[1]
        if J1[0]=='VirtualSegmentLen':VirtualSegmentLen=J1[1]
        if J1[0]=='PileSep':PileSep=J1[1]
        if J1[0]=='DyadTotal':DyadTotal=J1[1]
        if J1[0]=='SimpleOutputFileName':SimpleOutputFileName=J1[1]
        if J1[0]=='MultipleDataSets':MultipleDataSets=J1[1]
        if J1[0]=='SimpleOutputFormat':SimpleOutputFormat=J1[1]
        if J1[0]=='ReadMatchesTotal':ReadMatchesTotal=J1[1]
        if J1[0]=='ReadMatchesSense':ReadMatchesSense=J1[1]
        if J1[0]=='EndAntiSense':EndAntiSense=J1[1]
        if J1[0]=='MaxRead':MaxRead=J1[1]
        if J1[0]=='BaseMatchesTotal':BaseMatchesTotal=J1[1]
        if J1[0]=='ReferenceFileName':ReferenceFileName=J1[1]
        if J1[0]=='TargetStart':TargetStart=J1[1]
        if J1[0]=='ReadMatchesAntiSense':ReadMatchesAntiSense=J1[1]
        if J1[0]=='CompositionMatrix':CompositionMatrix=J1[1]
        if J1[0]=='DataUpper':DataUpper=J1[1]
        if J1[0]=='TransposeFile':TransposeFile=J1[1]
        if J1[0]=='SizeHistogramSense':SizeHistogramSense=J1[1]
        if J1[0]=='StartLenCollapse':StartLenCollapse=J1[1]
        if J1[0]=='DataPreCollapsed':DataPreCollapsed=J1[1]
        if J1[0]=='StartSense':StartSense=J1[1]
        if J1[0]=='Bin_Name':Bin_Name=J1[1]
        if J1[0]=='ExptDesc':ExptDesc=J1[1]
        if J1[0]=='BatchTable':BatchTable=J1[1]
        if J1[0]=='SizeHistogramTotal':SizeHistogramTotal=J1[1]
        if J1[0]=='StartTotal':StartTotal=J1[1]
        if J1[0]=='TabularOutput':TabularOutput=J1[1]
        if J1[0]=='Seed':Seed=J1[1]
        if J1[0]=='BaseMatchesSense':BaseMatchesSense=J1[1]
        if J1[0]=='MinRead':MinRead=J1[1]
        if J1[0]=='StartAntiSense':StartAntiSense=J1[1]
        if J1[0]=='StartCollapse':StartCollapse=J1[1]
        if J1[0]=='BaseBin':BaseBin=J1[1]
        if J1[0]=='ReadFileName':ReadFileName=J1[1]
        if J1[0]=='StoreIndex':StoreIndex=J1[1]
        if J1[0]=='IndexFilePath':IndexFilePath=J1[1]
        if J1[0]=='MinFreeGB':MinFreeGB=J1[1]


    ## Directories for user-specified input and output-- not really needed for most purposes

    ## Here are some general comments on the program
    ## Some Speed adjustments: Instead of AntiSense, use a "U" array that is a simple complement of T.
    ## QList and a number of other variables in the main finding loop are created during each loop.
    ## This could be avoided by having mutable lists (or variables) already existing before the loop.
    ## This includes the command where we check for a perfect sequence match by subsequencing T or U and comparing to S
    ## Also avoid all "Upper" statements and just insist that the input be upper case as soon as it is read from the disk

    Qprint(time0,LogFileName,'*'*30)
    if len(arguments1())>=1 and isfile(arguments1()[0]):
        fastqDumpPaths1+=arguments1()[0]
        Qprint(time0,LogFileName,'Script Name='+arguments1()[0])
        Qprint(time0,LogFileName,'Script Last-Modify time was '+asctime(localtime(stat(arguments1()[0])[8])))
    else:
        Qprint(time0,LogFileName,'No script name parameter was passed')
        Qprint(time0,LogFileName,'No script time parameter was passed')



    dashMM='-'*MaxMisMatch
    if ReadMatchesTotal:
        ReadMatchesSense=True
        ReadMatchesAntiSense=True



        
    ##  A set of routines that take input arguments in the form of a list, e.g.
    ##  ReadFile="MyFile.txt" MinMatch=30' and sets the relevant variables to this value.  Note that only variables defined
    ##  before this point in the program will be affected

    ## look for instructions, first in the arguments from the command line, then in a file passed in the command line,
    ## then in a file named "HQAlignBat.txt".  We will also save the entire run specification in a log
    ## file HQAlignLog.txt, along with the first several lines of each file.



    ## FileNames for program info 

    ##Use defaults, specifications, or dialog boxes to find the relevant files

    if not(ReferenceFileDir): ReferenceFileDir=getcwd()
    if ReferenceFileName.lower()=='ask':
        ReferenceFileName=UserAskFileNameWrite1("File with reference sequences?")
    ReferenceFile=open(ReferenceFileName,mode="rt")
    Qprint(time0,LogFileName,'ReferenceFileName='+ReferenceFile.name)
    ReferenceFileDir=dirname(abspath((ReferenceFile.name)))
    ReferenceFileName=basename(ReferenceFile.name)

    if not(MultipleDataSets):
        if ReadFileName.lower()=='ask':
            ReadFileName=UserAskFileNameRead1("File with high throughput sequence reads?")
        if ReadFileName.endswith('.sra'):
            ReadFileName,fastqDumpPath1=sraTranslate(ReadFileName,fastqDumpPath1)
        if ReadFileName.endswith('.gz'):
            if not gzipimported1:
                import gzip
            ReadFile=gzip.open(ReadFileName,"rt")
        else:
            ReadFile=open(ReadFileName,"rt")
        Qprint(time0,LogFileName,'ReadFileName='+ReadFile.name)
        ReadFileDir=dirname(abspath(ReadFile.name))
        ReadFileName=basename(ReadFile.name)
          
    if TabularOutput:
        if not(TabularOutFileDir): TabularOutFileDir=ReferenceFileDir
        if TabularOutputFileName.lower()=='ask':
            TabularOutputFileName=UserAskFileNameWrite1("File with high throughput sequence reads?")
        if TabularOutputFileName.lower()=='auto':
            TabularOutputFileName=join(TabularOutFileDir,ReferenceFileName +'_'+ExptDesc+now1+ 'HQ_Out.txt')
        outfile=open(TabularOutputFileName,mode="w")
        Qprint(time0,LogFileName,'TabularOutputFileName='+outfile.name)

    if PileUp:
        if not(PileUpFileDir): PileUpFileDir=ReferenceFileDir
        if PileUpFileName.lower()=='ask':
            PileUpFileName=UserAskFileNameWrite1("File for pile up output?")
        if PileUpFileName.lower()=='auto':
            PileUpFileName=join(PileUpFileDir,ReferenceFileName +'_'+ExptDesc+now1+'PileUp.m')
        PileUpFile=open(PileUpFileName,mode="w")
        Qprint(time0,LogFileName,'PileUpFileName='+PileUpFile.name)

    KeepReadNames=False  ##An option if one really wanted to keep the original sequencer read names, deprecated in favor of using the sequence
    if SimpleOutput:
        SimpleOutputFormat=SimpleOutputFormat.lower()
        bedformat=False
        pslformat=False
        hqaformat=True
        SimpleOutputSum=False
        if SimpleOutputFormat=='hqsum':
            SimpleOutputSum=True
        pslAddSeq=False
        if 'psl' in SimpleOutputFormat:
            pslformat=True
            hqaformat=False
            KeepReadNames=True
            if 'psls' in SimpleOutputFormat:
                pslAddSeq=True
        if 'bed' in SimpleOutputFormat:
            bedformat=True
            hqaformat=False
        if not(SimpleOutputFileDir): SimpleOutputFileDir=ReferenceFileDir
        if SimpleOutputFileName.lower()=='ask':
            SimpleOutputFileName=UserAskFileNameWrite1('.'+SimpleOutputFormat+" file for list of matches?")
        if SimpleOutputFileName.lower()=='auto':
            SimpleOutputFileName=join(SimpleOutputFileDir,ReferenceFileName +'_'+ExptDesc+ now1+ 'SimpleOutput.'+SimpleOutputFormat)
        SFile=open(SimpleOutputFileName,mode="w")
        Qprint(time0,LogFileName,'SimpleOutputFileName='+SFile.name)

        if SimpleOutputFormat=='hqsum':
            SFile.write('Number\tSource\tStart\tEnd\tOrientation\tMatch\tLength\tTag\tTarget\tGeneNumber\tExperiment_Description\r')
        if SimpleOutputFormat=='hqall':
            SFile.write('Source\tStart\tEnd\tOrientation\tMatch\tLength\tTag\tTarget\tGeneNumber\tExperiment_Description\r')
        if SimpleOutputFormat=='bed':
            SFile.write('track name=%s description=%s visibility=dense itemRGB=On priority=20\r'%(ExptDesc,ExptDesc)) ## a very basic set of bed/ucsc settings, change this for alternative display 
            #some colors for bed-file display on ucsc
            b0='0,0,255'
            b1='127,127,255'
            b2='220,220,255'
            r0='255,0,0'
            r1='255,127,127'
            r2='255,220,220'
            
    if BatchFilePath:
        try:
            copy2(BatchFilePath,join(BatchFileDir,BatchFileName.rsplit('.',1)[0]+'_'+ExptDesc+ now1+'.'+BatchFileType))
        except:
            Qprint(time0,LogFileName,'Unable to copy BatchFilePath')




    ## How to deal with multiple hits from the same query
    FirstHits_Only=False    ## only record the first hit by any given probe (so repetitive probes are all assigned to the first match in the genome)
    UniqueHits_Only=True   ## Require a unique hit, so toss all hits if there is more than one possible position in the reference file
    BestHits_Only=True     ## Best Hit(s) only (will take several hits if there are several, depending on the settings for FirstHits_Only,UniqueHits_Only,AllHits)

    '''Some Definitions:
        A "Match" is any target sequence that meets the cutoff criteria for a match (i.e. maximum of MinMatch baseswith a maximum of MM mismatches to start the read)
        A "Champion" is any target sequence that not only meets the cutoff criteria for a match but also does so more convincingly (with fewer mismathes) than any other candidate
        Note for the purposes of the following definitions, Champions are not assumed to be unique unless otherwise specified.

        Here are definitions of MultipleMatchMode:
        FirstMatch: Return 1st match in reference file that meets match criteria
        FirstChampion: Return 1st match in reference file that (i) meets match criteria, and (ii) has the best match score
        AllMatches: Return all matches in reference file that meet match criteria
        AllChampions: Return all matches in reference file that (i) meets match criteria, and (ii) tie for best match score
        UniqueMatchesOnly: Return only matches in the reference file that (i) meet the match criteria and (ii) where no other sites match criteria
        UniqueChampionsOnly: Return only matches in the reference file that (i) have best match scores and (iii) where no other target sites have that score
        SuperHits_n; Return only matches in the reference file that (i) meet match criteria and (ii) where any other match to the query has at least n additional mismatches'''


    SuperHits_Only=False    ## Only take hits that are at least #SuperHit better than the next best hit.  This rules out slight mismatches to very abundant sequences
    SuperHit=1              ## This requires a "Super" hit with mismatches<=MM-SuperHit , so that an imperfect hit could have been detected
    MMDetect=MaxMisMatch+0

    ## MaxMisMatch is the maximum number of mismatches allowed in any sequence that will be reported
    ## MMDetect is the maximum number that will be browsed.  This is particularly relevant for SuperHit analysis, where it is necessary to look at sequences with additional
    ## mismatches  (>MaxMisMatch) to be sure a superhit is not incorrectly assigned
    MultipleMatchMode=MultipleMatchMode.lower()
    if MultipleMatchMode.startswith('firstmatch'):
        FirstHits_Only=True; UniqueHits_Only=False; BestHits_Only=False
    if MultipleMatchMode.startswith('firstchamp'):
        FirstHits_Only=True; UniqueHits_Only=False; BestHits_Only=True
    if MultipleMatchMode.startswith('allmatch'):
        FirstHits_Only=False; UniqueHits_Only=False; BestHits_Only=False
    if MultipleMatchMode.startswith('allchamp'):
        FirstHits_Only=False; UniqueHits_Only=False; BestHits_Only=True
    if MultipleMatchMode.startswith('uniquematch'):
        FirstHits_Only=False; UniqueHits_Only=True; BestHits_Only=False
    if MultipleMatchMode.startswith('uniquechamp'):
        FirstHits_Only=False; UniqueHits_Only=True; BestHits_Only=True
    if MultipleMatchMode.startswith('superhit'):
        mmm0=MultipleMatchMode
        SuperHits_Only=True
        SuperHit='n'
        while SuperHit[0].isdigit():
            SuperHit=mmo[-1]+SuperHit
            mmo=mmo[:-1]
        SuperHit=SuperHit[1:-1]
        if not SuperHit:
            SuperHit='2'
        SuperHit=int(SuperHit)
        MMDetect=MaxMisMatch+SuperHit
        
    IndexRead=False        
    ttx11=strftime('_RefModDate-%m-%d-%y_RefModTime-%H-%M-%S_',localtime(stat(ReferenceFile.name)[8]))
    if IndexFilePath=='' or IndexFilePath.lower()=='auto':
        IndexFileName='HQAIndexFromRef_%s_%sSeed-%i_MM-%i_Gr-%i_%s%i.hqa2'%(ReferenceFileName,ttx11,Seed,MaxMisMatch,Granularity,str(GeneBin)[:1]+str(BinBin)[:1],StopAfterReference)
        IndexFilePath=join(ReferenceFileDir,IndexFileName)
    if isfile(IndexFilePath):
        Qprint(time0,LogFileName,'Trying Read From Precompiled Index.  IndexFileName='+IndexFileName)
        try:
            pf1=open(IndexFilePath,mode='rb')
            pp1=cPickle.Unpickler(pf1)
            Bins1=pp1.load()
            Bases1=pp1.load()
            Genes1=pp1.load()
            LT=pp1.load()
            LY=pp1.load()
            GeneNameA=pp1.load()
            RefAllUpper=pp1.load()
            FirstTry=array('l',[])
            FirstTry.fromfile(pf1,LY+1)
            NextTry=array('l',[])   
            NextTry.fromfile(pf1,LY+1)
            Gene_to_Start=array('L',[])
            Gene_to_Start.fromfile(pf1,Genes1)
            Gene_to_Len=array('L',[])
            Gene_to_Len.fromfile(pf1,Genes1)
            if GeneBin:
                if Genes1<256:
                    Address_to_GeneNum=array('B',[])
                elif Genes1<65536:
                    Address_to_GeneNum=array('H',[])
                else:
                    Address_to_GeneNum=array('L',[])
                Address_to_GeneNum.fromfile(pf1,LT+1)
            if BinBin:
                if Bins1<256:
                    Address_to_BinNum=array('B',[])
                elif Bins1<65536:
                    Address_to_BinNum=array('H',[])
                else:
                    Address_to_BinNum=array('L',[])
                Bin_to_GStart=array('L',[]) 
                Bin_to_TStart=array('L',[]) 
                Bin_to_Len=array('L',[])  
                Bin_to_Gene=array('L',[])  
                Bin_Precedence=array('L',[])
                Address_to_BinNum.fromfile(pf1,LT+1)  
                Bin_to_GStart.fromfile(pf1,Bins1) 
                Bin_to_TStart.fromfile(pf1,Bins1) 
                Bin_to_Len.fromfile(pf1,Bins1)  
                Bin_to_Gene.fromfile(pf1,Bins1)  
                Bin_Precedence.fromfile(pf1,Bins1)
            if float(sys.version[:3])<3.0:
                T=pf1.read(LT)
                U=pf1.read(LT)
            else:
                T=str(pf1.read(LT))[2:-1]
                U=str(pf1.read(LT))[2:-1]                
            pf1.close()
            if RefAllUpper:
                Tup=T
                Uup=U
            else:
                Tup=T.upper()
                Uup=U.upper()
            IndexRead=True
            Qprint(time0,LogFileName,'Read From Precompiled Index '+IndexFileName+' successful.')
        except:
            IndexRead=False
            Qprint(time0,LogFileName,'Read From Precompiled Index '+IndexFileName+' unsuccessful, commencing standard index build.')
    if not(IndexRead):

        ## Quickly open the source (reference) file to count bases and bins to set array sizes etc
        ## For the full output, the counts of hits are separated into bins, which can be the whole sequence, individual genes (defined by FASTA format), or fixed bins of size "Granularity"
        ## Bins1 is the total number of bins, and thes include
        ## TotalBin: Record a bin (at the beginning of the list) that covers everything
        ## GeneBin: Record bins of each gene (all hits)
        ## BinBin: Record bins of each segment with length 'granularity'
        ## BaseBin: Record a bin for every base

        ## Now assemble a more complete index, also T (an array of all references sequences with each "gene" separated by a *
        ## and U, which is the antisense of T
        Bins1=0
        StartSequence=True
        T=[]
        LX=0 ##  total number of bases parsed in this gene
        LXT=0 ## total number of bases parsed in all genes
        CurGene=0
        StartSequence=True
        GeneNameA=[]

        for D in ReferenceFile:
            if D.strip()=='':continue
            if (D[0]=='>') or StartSequence:
                if D[0]=='>':
                    En=D[1:].strip()
                else:
                    En=ReferenceFileName
                GeneNameA.append(En[:])
                CurGene+=1
                if not(StartSequence):
                    T.append(dashMM)
                T.append('*') ## so T[0] is an asterisk,
                T.append(dashMM)
                StartSequence=False
                if BinBin: Bins1+= 1+((LX-1) // Granularity)
                LX=0
                if D[0]=='>': continue
            E=D.translate(filterplus).replace(' ','')
            Badd=len(E)
            LX+=Badd
            LXT+=Badd
            T.append(E)
            if StopAfterReference>0 and LXT>=StopAfterReference:        
                if LXT>StopAfterReference:
                    T=T[:-LXT+StopAfterReference]
                LXT=StopAfterReference+0
                break
        T.append(dashMM)
        T.append('*')
        if BinBin: Bins1+= 1+((LX-1) // Granularity)

        Bases1=LXT+0
        Genes1=CurGene+0

        ReferenceFile.close()

        T=''.join(T)  ## joining a list of strings appears to be faster than appending line by line, with thanks to Oliver Crow
        GeneNameA=GeneNameA[:] ## to clean up memory fragmentation
        LT=len(T)

        ##Note that the indices in the all arrays below are zero based (so zero is the index for the first gene or first base).
        ##The values for Address to GeneNum zero based (so zero is the first gene)
        if GeneBin:
            if Genes1<256:
                Address_to_GeneNum=array('B',[0]*(LT+1))  ## which gene is a given address in, this yields the ordinal number (zero based) of a gene
            elif Genes1<65536:
                Address_to_GeneNum=array('H',[0]*(LT+1))  ## which gene is a given address in, this yields the ordinal number (zero based) of a gene
            else:
                Address_to_GeneNum=array('L',[0]*(LT+1))  ## which gene is a given address in, this yields the ordinal number (zero based) of a gene
        ##Note that the values in all of the arrays below are zero based), so 0 is the first character in T (an asterisk) and 1 is the first base
        Gene_to_Start=array('L',[0]*Genes1) ## where each gene starts in T
        Gene_to_Len=array('L',[0]*Genes1) ## length of each gene in T
        if BinBin:
            if Bins1<256:
                Address_to_BinNum=array('B',[0]*(LT+1))  ## which bin is a given address in
            elif Bins1<65536:
                Address_to_BinNum=array('H',[0]*(LT+1))  ## which bin is a given address in
            else:
                Address_to_BinNum=array('L',[0]*(LT+1))  ## which bin is a given address in
            Bin_to_GStart=array('L',[0]*Bins1) ## where the bin starts in the gene ; zero based
            Bin_to_TStart=array('L',[0]*Bins1) ## where the bin starts in the T array ; zero based
            Bin_to_Len=array('L',[0]*Bins1)  ## Bin_Length
            Bin_to_Gene=array('L',[0]*Bins1)  ## Ordinal number of the gene associated with each bin
            Bin_Precedence=array('L',[0]*Bins1)  ## Ordinal number of the bin among bins within a gene ; zero based

        ii1=0                ## position in T
        gg1=0                ## current gene
        hh1=-MaxMisMatch     ## current position in gene
        bb1=0                ## current bin
        cc1=-MaxMisMatch     ## current position in bin
        firstrun=True
        for C in T:
            if C=='*':
                if not(firstrun):
                    Gene_to_Len[gg1]=hh1-MaxMisMatch
                    gg1+=1
                    if BinBin and cc1>0:
                        Bin_to_Len[bb1]=cc1-MaxMisMatch
                        bb1+=1
                        cc1=-MaxMisMatch               
                else:
                    firstrun=False
                if ii1==LT-1: break
                Gene_to_Start[gg1]=ii1+1+MaxMisMatch
                hh1=-MaxMisMatch
                if GeneBin:
                    Address_to_GeneNum[ii1]=gg1
                if BinBin:
                    Bin_to_Gene[bb1]=gg1
                    Bin_to_GStart[bb1]=0
                    Bin_Precedence[bb1]=0
                    Bin_to_TStart[bb1]=ii1+1+MaxMisMatch
                    Address_to_BinNum[ii1]=bb1
            else:
                if BinBin:
                    if cc1==Granularity:
                        Bin_to_Len[bb1]=cc1
                        cc1=0
                        bb1+=1
                        Bin_Precedence[bb1]=1+Bin_Precedence[bb1-1]
                        Bin_to_Gene[bb1]=gg1
                        Bin_to_GStart[bb1]=hh1
                        Bin_to_TStart[bb1]=ii1
                    Address_to_BinNum[ii1]=bb1
                    cc1+=1
                if GeneBin:
                    Address_to_GeneNum[ii1]=gg1
                hh1+=1
            ii1+=1
                
        LY0=4**Seed
        Seed2m=int(LY0/4)
        FirstTry=array('l',[0]*(LY0+2*LT))  
        NextTry=array('l',[0]*(LY0+2*LT))   
        LX=0
        LY=LY0+0
        LY1=LY-1
        LS=0

        Qprint(time0,LogFileName,"Indexed referenced bases to: " + str(LXT))
        U=AntiSense(T,filterminus)
        LX=0
        Tryseed_a=0
        Tryseed_s=0
        ASeed2m=0*Seed2m
        TSeed2m=1*Seed2m
        CSeed2m=2*Seed2m
        GSeed2m=3*Seed2m
        Tup=T.upper()
        # save some memory if all upper case to start with
        if Tup==T:
            Tup=T
            Uup=U
            RefAllUpper=True
        else:
            Uup=U.upper()
            RefAllUpper=False

        for C in Uup:
            LS+=1
            if C=='*' or C=='-':
                LX=0
            else:
                Tryseed_s>>=2
                Tryseed_a=(Tryseed_a<<2) & LY1
                if C=='A':
                    Tryseed_s+=TSeed2m
                elif C=='T':
                    Tryseed_a+=1
                elif C=='C':
                    Tryseed_s+=GSeed2m
                    Tryseed_a+=2
                else:
                    Tryseed_s+=CSeed2m
                    Tryseed_a+=3
                LX+=1
                if LX>=Seed:
                        if FirstTry[Tryseed_a]==0:
                            FirstTry[Tryseed_a]=-LS+Seed-1
                        else:
                            FirstTry[LY]=FirstTry[Tryseed_a]
                            NextTry[LY]=NextTry[Tryseed_a]
                            FirstTry[Tryseed_a]=-LS+Seed-1
                            NextTry[Tryseed_a]=LY
                            LY+=1
                        if FirstTry[Tryseed_s]==0:
                            FirstTry[Tryseed_s]=LT-LS+1
                        else:
                            FirstTry[LY]=FirstTry[Tryseed_s]
                            NextTry[LY]=NextTry[Tryseed_s]
                            FirstTry[Tryseed_s]=LT-LS+1
                            NextTry[Tryseed_s]=LY
                            LY+=1

        FirstTry=FirstTry[:LY+1]
        NextTry=NextTry[:LY+1]
        ## Important- phasing of the FirstTry and NextTry arrays
        ## entries in FirstTry are indexed by base sequence for the first 4**Seed sequences, then operate above that as a linked list
        ## FirstTry==0 meand this seed sequence has never been encountered, FirstTry=-X means that there is an match to the seed starting
        ## at base X+1 (zero based) in Uup.  FirstTry=+X means there is a match to the seed starting at base X+1 (zero based) in Tup
        ## For entries numbered >=4**Seed, these are arbitrary addresses essentially in the form of a linked list.
        ## Next Try arrays tell the system where to look following the first try for any given seed sequence.


    Qprint(time0,LogFileName,"Finished Indexing")

    if TabularOutput:
        if Bin_Name:
            outfile.write('FeatureName')
            if TotalBin:
                outfile.write('\t'+ReferenceFileName)
            if GeneBin:
                for X in GeneNameA:
                    outfile.write('\t'+X)
            if BinBin:
                for (Y,Z) in zip(Bin_to_Gene,Bin_Precedence):
                    outfile.write('\t'+GeneNameA[Y]+'_bin_'+str(Z))
            if BaseBin:
                for b11 in xrange(LT):
                    if T[b11]!='*':
                        if T[b11]!='-':
                            if GeneBin:
                                outfile.write('\t'+GeneNameA[Address_to_GeneNum[b11]]+'_base_'+str(b11-Gene_to_Start[Address_to_GeneNum[b11]]+1))
                            else:
                                outfile.write('\t'+ReferenceFileName+'_base_'+str(b11-MaxMisMatch))
                        else:
                            if GeneBin:
                                outfile.write('\t'+GeneNameA[Address_to_GeneNum[b11]]+'_flank_'+str(b11-Gene_to_Start[Address_to_GeneNum[b11]]+1))
                            else:
                                outfile.write('\t'+ReferenceFileName+'_flank_'+str(b11-MaxMisMatch))
            outfile.write('\r')
        if Bin_Size: 
            outfile.write('FeatureSize')
            if TotalBin:
                outfile.write('\t'+str(Bases1))
            if GeneBin:
                for X in Gene_to_Len:
                    outfile.write('\t'+str(X))
            if BinBin:
                for X in Bin_to_Len:
                    outfile.write('\t'+str(X))
            if BaseBin:
                for b11 in xrange(LT):
                    if T[b11]!='*':
                        outfile.write('\t1')
            outfile.write('\r')
        if Bin_Start: 
            outfile.write('FeatureStart')
            if TotalBin:
                outfile.write('\t1')
            if GeneBin:
                for X in Gene_to_Start:
                    outfile.write('\t1')
            if BinBin:
                for X in Bin_to_GStart:
                    outfile.write('\t'+str(X+1))
            if BaseBin:
                for b11 in xrange(LT):
                    if T[b11]!='*':
                        if GeneBin:
                            outfile.write('\t'+str(b11-Gene_to_Start[Address_to_GeneNum[b11]]+1))
                        else:
                            outfile.write('\t'+str(b11-MaxMisMatch))
            outfile.write('\r')
        if Bin_End: 
            outfile.write('FeatureEnd')
            if TotalBin:
                outfile.write('\t'+str(Bases1))
            if GeneBin:
                for Y in Gene_to_Len:
                    outfile.write('\t'+str(Y))
            if BinBin:
                for (X,Y) in zip(Bin_to_GStart,Bin_to_Len):
                    outfile.write('\t'+str(X+Y))
            if BaseBin:
                for b11 in xrange(LT):
                    if T[b11]!='*':
                        if GeneBin:
                            outfile.write('\t'+str(b11-Gene_to_Start[Address_to_GeneNum[b11]]+1))
                        else:
                            outfile.write('\t'+str(b11-MaxMisMatch))
            outfile.write('\r')
        if Bin_Sequence: 
            outfile.write('FeatureSeq')
            if TotalBin:
                outfile.write('\t'+T[1+MaxMisMatch:-1-MaxMisMatch])
            if GeneBin:
                for (X,Y) in zip(Gene_to_Start,Gene_to_Len):
                    outfile.write('\t'+T[X:X+Y])
            if BinBin:
                for (X,Y) in zip(Bin_to_TStart,Bin_to_Len):
                    outfile.write('\t'+T[X:X+Y])
            if BaseBin:
                for X in xrange(LT):
                    if T[X]!='*':
                        outfile.write('\t'+T[X])
            outfile.write('\r')
        if Bin_Composition:
            for J in 'ACGT':
                outfile.write(J+'_in_Feature')
                if TotalBin:
                    outfile.write('\t'+str(Tup.count(J)))
                if GeneBin:
                    for (X,Y) in zip(Gene_to_Start,Gene_to_Len):
                        outfile.write('\t'+str(Tup[X+1:X+Y+1].count(J)))
                if BinBin:
                    for (X,Y) in zip(Bin_to_TStart,Bin_to_Len):
                        outfile.write('\t'+str(Tup[X+1:X+Y+1].count(J)))
                if BaseBin:
                    for X in xrange(LT):
                        if T[X]!='*':
                            outfile.write('\t'+str(int(Tup[X]==J)))
                outfile.write('\r')

    B1='ACGT'
    if not(MultipleDataSets):
        BatchTable='MultipleDataSets=False'
    if type(BatchTable)==str:
        NewBatchTable=[]
        BatchTable=BatchTable.split('\n')
        for CommandLine in BatchTable:
            NewBatchTable.append([])
            CommandLine=CommandLine.replace('\t',',')
            Commands=CommandLine.split(',')
            for e in Commands:
                if "=" in e:
                    e=e.split('=')
                    e[0]=e[0].strip().strip('"'+"'")
                    e[1]=e[1].strip().strip('"'+"'")
                    NewBatchTable[-1].append([e[0],e[1]])
        BatchTable=NewBatchTable

    for CommandLine in BatchTable:
        for J1 in CommandLine:
            try:
                CommandVarType1=type(HardLocals[J1[0]])
                if CommandVarType1==bool and J1[1].lower()=='false':
                    J1[1]=''
                try:
                    res1=CommandVarType1(J1[1])
                    Qprint(time0,LogFileName,"Setting "+str(J1[0])+'='+str(CommandVarType1(J1[1])))
                except:
                    try:
                        res1=CommandVarType1(float(J1[1]))
                        Qprint(time0,LogFileName,"Setting "+str(J1[0])+'='+str(J1[1]))
                    except:
                        Qprint(time0,LogFileName,"Error in setting "+str(J1d[0])+'='+str(J1[1])+'.  Will continue as best I can.  Check your batch command file.')
                        try:
                            res1=J1[1]
                        except:
                            pass
            except:
                Qprint(time0,LogFileName,"Error in setting "+str(J1[0])+'='+str(J1[1])+'.  Will continue as best I can.  Check your batch command file.')

            J1[1]=res1

            if J1[0]=='RequireEndLinker':RequireEndLinker=J1[1]
            if J1[0]=='ExtensionOK':ExtensionOK=J1[1]
            if J1[0]=='ReGran':ReGran=J1[1]
            if J1[0]=='SampleDesc':SampleDesc=J1[1]
            if J1[0]=='ReportInterval':ReportInterval=J1[1]
            if J1[0]=='StopAfterQuery':StopAfterQuery=J1[1]
            if J1[0]=='SimpleOutput':SimpleOutput=J1[1]
            if J1[0]=='PileUp':PileUp=J1[1]
            if J1[0]=='MinMatch':MinMatch=J1[1]
            if J1[0]=='CompositionIndexed':CompositionIndexed=J1[1]
            if J1[0]=='Composition_End':Composition_End=J1[1]
            if J1[0]=='MaxMisMatch':MaxMisMatch=J1[1]
            if J1[0]=='CompositionTupleLen':CompositionTupleLen=J1[1]
            if J1[0]=='KeepFullQuery':KeepFullQuery=J1[1]
            if J1[0]=='Base1SizeHistogram':Base1SizeHistogram=J1[1]
            if J1[0]=='TargetLength':TargetLength=J1[1]
            if J1[0]=='StartBarcode':
                StartBarcode=J1[1]
                FivePrimeTrim=0
                while StartBarcode and StartBarcode[0].upper()=='N':
                    FivePrimeTrim+=1
                    StartBarcode=StartBarcode[1:]
                if StartBarcode:
                    sbc=StartBarcode[0]                
            if J1[0]=='EndLinker':EndLinker=J1[1]
            if J1[0]=='VirtualSegmentLen':VirtualSegmentLen=J1[1]
            if J1[0]=='MaxRead':MaxRead=J1[1]
            if J1[0]=='DataUpper':DataUpper=J1[1]
            if J1[0]=='StartLenCollapse':StartLenCollapse=J1[1]
            if J1[0]=='DataPreCollapsed':DataPreCollapsed=J1[1]
            if J1[0]=='ExptDesc':ExptDesc=J1[1]
            if J1[0]=='MinRead':MinRead=J1[1]
            if J1[0]=='StartCollapse':StartCollapse=J1[1]
            if J1[0]=='ReadFileName':ReadFileName=J1[1]
            if J1[0]=='MultipleMatchMode':
                MultipleMatchMode=J1[1]
                if MultipleMatchMode.startswith('firstmatch'):
                    FirstHits_Only=True; UniqueHits_Only=False; BestHits_Only=False
                if MultipleMatchMode.startswith('firstchamp'):
                    FirstHits_Only=True; UniqueHits_Only=False; BestHits_Only=True
                if MultipleMatchMode.startswith('allmatch'):
                    FirstHits_Only=False; UniqueHits_Only=False; BestHits_Only=False
                if MultipleMatchMode.startswith('allchamp'):
                    FirstHits_Only=False; UniqueHits_Only=False; BestHits_Only=True
                if MultipleMatchMode.startswith('uniquematch'):
                    FirstHits_Only=False; UniqueHits_Only=True; BestHits_Only=False
                if MultipleMatchMode.startswith('uniquechamp'):
                    FirstHits_Only=False; UniqueHits_Only=True; BestHits_Only=True
                if MultipleMatchMode.startswith('superhit'):
                    mmm0=MultipleMatchMode
                    SuperHits_Only=True
                    SuperHit='n'
                    while SuperHit[0].isdigit():
                        SuperHit=mmo[-1]+SuperHit
                        mmo=mmo[:-1]
                    SuperHit=SuperHit[1:-1]
                    if not SuperHit:
                        SuperHit='2'
                    SuperHit=int(SuperHit)
                    MMDetect=MaxMisMatch+SuperHit

        if not(ReadFileDir): ReadFileDir=ReferenceFileDir
        if ReadFileName.lower()=='ask':
            FileQuery="File with high throughput sequencing reads"
            if StartBarcode: FileQuery+='  Barcode='+StartBarcode
            ReadFileName=UserGetFileNameRead1(FileQuery)
        if ReadFileName.endswith('.sra'):
            ReadFileName,fastqDumpPath1=sraTranslate(ReadFileName,fastqDumpPath1)
        if ReadFileName.endswith('.gz'):
            if not gzipimported1:
                import gzip
            ReadFile=gzip.open(ReadFileName,"rt")
        else:
            ReadFile=open(ReadFileName,"rt")
        Qprint(time0,LogFileName,'ReadFileName='+ReadFile.name)
        
        TargetStart=TargetStart.upper()
        if TargetStart=='NONE':
            TargetStart=''


        ## Set up arrays to store info from datamatching
        if TotalBin:
            if StartSense: Start_SenseT=0
            if StartAntiSense: Start_AntiSenseT=0
            if StartTotal: Start_TotalT=0
            if DyadSense: Dyad_SenseT=0
            if DyadAntiSense: Dyad_AntiSenseT=0
            if DyadTotal: Dyad_TotalT=0
            if EndSense: End_SenseT=0
            if EndAntiSense: End_AntiSenseT=0
            if EndTotal: End_TotalT=0
            if CoverageSense: Coverage_SenseT=0
            if CoverageAntiSense: Coverage_AntiSenseT=0
            if CoverageTotal: Coverage_TotalT=0
            if SizeHistogramSense: SizeHistogram_SenseT={}        ## not very efficient for large numbers of bins
            if SizeHistogramAntiSense: SizeHistogram_AntiSenseT={}       
            if SizeHistogramTotal: SizeHistogram_TotalT={}       
            if BaseMatchesSense: BaseMatches_SenseT={}        ## not very efficient for large numbers of bins
            if BaseMatchesAntiSense: BaseMatches_AntiSenseT={}      
            if BaseMatchesTotal: BaseMatches_TotalT={}        
            if ReadMatchesSense or PileUp:
                ReadMatches_SenseT={}        ## not very efficient for large numbers of bins
                MMSenseOffsetT=0
            if ReadMatchesAntiSense or PileUp:
                ReadMatches_AntiSenseT={}        
                MMAntiSenseOffsetT=0
            if CompositionMatrix: CompositionMatrixDT={} ## not very efficient for large numbers of bins
        if GeneBin:
            if StartSense: Start_SenseG=array('L',[0]*Genes1)
            if StartAntiSense: Start_AntiSenseG=array('L',[0]*Genes1)
            if StartTotal: Start_TotalG=array('L',[0]*Genes1)
            if DyadSense: Dyad_SenseG=array('L',[0]*Genes1)
            if DyadAntiSense: Dyad_AntiSenseG=array('L',[0]*Genes1)
            if DyadTotal: Dyad_TotalG=array('L',[0]*Genes1)
            if EndSense: End_SenseG=array('L',[0]*Genes1)
            if EndAntiSense: End_AntiSenseG=array('L',[0]*Genes1)
            if EndTotal: End_TotalG=array('L',[0]*Genes1)
            if CoverageSense: Coverage_SenseG=array('L',[0]*Genes1)
            if CoverageAntiSense: Coverage_AntiSenseG=array('L',[0]*Genes1)
            if CoverageTotal: Coverage_TotalG=array('L',[0]*Genes1)
            if SizeHistogramSense: SizeHistogram_SenseG=[{} for inde in range(Genes1)]        ## not very efficient for large numbers of bins
            if SizeHistogramAntiSense: SizeHistogram_AntiSenseG=[{} for inde in range(Genes1)]        
            if SizeHistogramTotal: SizeHistogram_TotalG=[{} for inde in range(Genes1)]        
            if BaseMatchesSense: BaseMatches_SenseG=[{} for inde in range(Genes1)]        ## not very efficient for large numbers of bins
            if BaseMatchesAntiSense: BaseMatches_AntiSenseG=[{} for inde in range(Genes1)]        
            if BaseMatchesTotal: BaseMatches_TotalG=[{} for inde in range(Genes1)]        
            if ReadMatchesSense or PileUp:
                ReadMatches_SenseG=[{} for inde in range(Genes1)]        
                MMSenseOffsetG=array('L',[0]*Genes1)
            if ReadMatchesAntiSense or PileUp:
                ReadMatches_AntiSenseG=[{} for inde in range(Genes1)]        
                MMAntiSenseOffsetG=array('L',[0]*Genes1)
            if CompositionMatrix: CompositionMatrixDG=[{} for inde in range(Genes1)] ## not very efficient for large numbers of bins
        if BinBin:
            if StartSense: Start_SenseB=array('L',[0]*Bins1)
            if StartAntiSense: Start_AntiSenseB=array('L',[0]*Bins1)
            if StartTotal: Start_TotalB=array('L',[0]*Bins1)
            if DyadSense: Dyad_SenseB=array('L',[0]*Bins1)
            if DyadAntiSense: Dyad_AntiSenseB=array('L',[0]*Bins1)
            if DyadTotal: Dyad_TotalB=array('L',[0]*Bins1)
            if EndSense: End_SenseB=array('L',[0]*Bins1)
            if EndAntiSense: End_AntiSenseB=array('L',[0]*Bins1)
            if EndTotal: End_TotalB=array('L',[0]*Bins1)
            if CoverageSense: Coverage_SenseB=array('L',[0]*Bins1)
            if CoverageAntiSense: Coverage_AntiSenseB=array('L',[0]*Bins1)
            if CoverageTotal: Coverage_TotalB=array('L',[0]*Bins1)
            if SizeHistogramSense: SizeHistogram_SenseB=[{} for inde in range(Bins1)]        ## not very efficient for large numbers of bins
            if SizeHistogramAntiSense: SizeHistogram_AntiSenseB=[{} for inde in range(Bins1)]        
            if SizeHistogramTotal: SizeHistogram_TotalB=[{} for inde in range(Bins1)]        
            if BaseMatchesSense: BaseMatches_SenseB=[{} for inde in range(Bins1)]        ## not very efficient for large numbers of bins
            if BaseMatchesAntiSense: BaseMatches_AntiSenseB=[{} for inde in range(Bins1)]        
            if BaseMatchesTotal: BaseMatches_TotalB=[{} for inde in range(Bins1)]        
            if ReadMatchesSense or PileUp:
                ReadMatches_SenseB=[{} for inde in range(Bins1)]        ## not very efficient for large numbers of bins
                MMSenseOffsetB=array('L',[0]*Bins1)
            if ReadMatchesAntiSense or PileUp:
                ReadMatches_AntiSenseB=[{} for inde in range(Bins1)]        
                MMAntiSenseOffsetB=array('L',[0]*Bins1)
            if CompositionMatrix: CompositionMatrixDB=[{} for inde in range(Bins1)] ## not very efficient for large numbers of bins
        if BaseBin:
            if StartSense: Start_SenseBa=array('L',[0]*LT)
            if StartAntiSense: Start_AntiSenseBa=array('L',[0]*LT)
            if StartTotal: Start_TotalBa=array('L',[0]*LT)
            if DyadSense: Dyad_SenseBa=array('L',[0]*LT)
            if DyadAntiSense: Dyad_AntiSenseBa=array('L',[0]*LT)
            if DyadTotal: Dyad_TotalBa=array('L',[0]*LT)
            if EndSense: End_SenseBa=array('L',[0]*LT)
            if EndAntiSense: End_AntiSenseBa=array('L',[0]*LT)
            if EndTotal: End_TotalBa=array('L',[0]*LT)
            if CoverageSense: Coverage_SenseBa=array('L',[0]*LT)
            if CoverageAntiSense: Coverage_AntiSenseBa=array('L',[0]*LT)
            if CoverageTotal: Coverage_TotalBa=array('L',[0]*LT)
            if SizeHistogramSense: SizeHistogram_SenseBa=[{} for inde in range(LT)]        ## not very efficient for large numbers of bins
            if SizeHistogramAntiSense: SizeHistogram_AntiSenseBa=[{} for inde in range(LT)]        
            if SizeHistogramTotal: SizeHistogram_TotalBa=[{} for inde in range(LT)]        
            if BaseMatchesSense: BaseMatches_SenseBa=[{} for inde in range(LT)]        ## not very efficient for large numbers of bins
            if BaseMatchesAntiSense: BaseMatches_AntiSenseBa=[{} for inde in range(LT)]        
            if BaseMatchesTotal: BaseMatches_TotalBa=[{} for inde in range(LT)]        
            if ReadMatchesSense or PileUp: ReadMatches_SenseBa=[{} for inde in range(LT)]        ## not very efficient for large numbers of bins
            if ReadMatchesAntiSense or PileUp: ReadMatches_AntiSenseBa=[{} for inde in range(LT)]        
            if CompositionMatrix: CompositionMatrixDBa=[{} for inde in range(LT)] ## not very efficient for large numbers of bins

        if StartLenCollapse: StartLenCollapseD=[{} for inde in range(LT)]
        if StartCollapse: StartCollapseA=array('B',[0]*LT)  ## Value will be 1 for sense hit, 2 for antisense hit, 3 for both, 0 for neither
        if SimpleOutput and SimpleOutputSum: SimpleOutputD={}
        ## initialize some variables for the main "finding" loop 
        P1=0
        P2=0
        M=0
        N=0

        LBC=len(StartBarcode)
        Experiment=SampleDesc+' ReadFile='+basename(ReadFile.name)
        Experiment+=' RefFile='+basename(ReferenceFile.name)
        if StartBarcode: Experiment+=' Barcode='+StartBarcode
        if EndLinker: Experiment+=' EndLinker='+EndLinker

        ## Warn the user if the seed is too long for the requested number of mismatches (some matches would be lost)
        if Seed*(MMDetect+1)>MinMatch:
            print("Warning, Seed too long for number of mismatches")
            print("Seed="+str(Seed))
            print("MinRead="+str(MinMatch))
            print("MisMatch="+str(MMDetect))
            print("Seed should be <1/(MMDetect+1) of MinRead")
            print("Will continue as is, but some matches my be missed")
            print("Suggest: For definitive pattern matching, Run again with Seed equal to "+str(MinRead//(MMDetect+1)))

        rangelist=range(0,(MMDetect+1)*Seed,Seed)
        rangearray=[]
        for jnde in range((MMDetect+1)*Seed):
            rangearray.append(range(jnde,jnde+Seed))
        FullAlign=(MMDetect>0) or ExtensionOK
        ##ReadName1=ExptDesc+''
        Mukltiplicity1=1
        for S in ReadFile:
            S0=S[0]
            if S0=='@' or S0=='>':
                if DataPreCollapsed and not(StartCollapse) and not(StartLenCollapse):
                    try:
                        Multiplicity1=int(S0.strip().split('-')[-1])
                    except:
                        print('Warning: Not able to parse multiplicity information for read named ' +S+'.  Setting multiplicity to 1 as default')
                        Multiplicity1=1
##                if KeepReadNames:
##                    ReadName1=S[1:]
                continue
            if S0=='+':
                next(ReadFile,'')
                continue
            if not(DataUpper):
                S=S.filter(filterplus).filter(filterplus,' ')
                S0=S[0]
            if FivePrimeTrim>0:
                S=S[FivePrimeTrim:]
            if StartBarcode:
                if not(S[0]==sbc):continue
                if not(S.startswith(StartBarcode)): continue
                S=S[LBC:]
            if EndLinker:
                L=S.rfind(EndLinker)
                if L==-1:
                    if RequireEndLinker:
                        continue
                    else:
                        S=S.rstrip()
                        L=len(S)
                else:
                    S=S[:L]
            else:
                S=S.rstrip()
                L=len(S)
            if L>MaxRead:
                S=S[:MaxRead]
                L=MaxRead+0
            if TargetLength and TargetLength!=L: continue
            if TargetStart and not(S.startswith(TargetStart)): continue            
            if L>=MinRead:
                LeastMM=MMDetect+1
                MM1=0
                QList=[]
                SampleDone=False
                for jnde in rangelist:
                    if jnde+Seed>L:break
                    Tryseed=0
                    for inde in rangearray[jnde]:
                        Tryseed<<=2
                        C=S[inde]
                        if C=='T':
                            Tryseed+=1
                        elif C=='C':
                            Tryseed+=2
                        elif C=='G':
                            Tryseed+=3
                    Q=FirstTry[Tryseed]
                    while Q!=0 and not(SampleDone):
                        Match=False
                        if (Q>0 and Tup[Q-1-jnde:Q+L-1-jnde]==S) or (Q<0 and Uup[-Q-1-jnde:-Q+L-1-jnde]==S):
                            Match=True
                            MM1=0
                            inde=L+0
                        if FullAlign and not(Match):
                            MM1=0
                            if Q>0:
                                for inde in range(L):
                                    if Q-1+inde-jnde<1:
                                         MM1+=1
                                    elif Tup[Q-1+inde-jnde]!=S[inde]:
                                        if Tup[Q-1+inde-jnde]=='*':
                                            if (inde<MMDetect) and (Q-1+inde-jnde)!=LT:
                                                MM1=inde+0
                                            else:
                                                break
                                        MM1+=1 
                                    elif inde>=MinMatch+MM1-1 and ExtensionOK:
                                        Match=True
                                        break
                                    if MM1>MMDetect:
                                        break
                                else:      
                                    Match=True   ## this statement only runs of we get to the end of the loop (inde=L-1) and a match hasn't been excluded by number of mismatches
                            else:
                                for inde in range(L):
                                    if -Q-1+inde-jnde<1:
                                         MM1+=1
                                    elif Uup[-Q-1+inde-jnde]!=S[inde]:
                                        if Uup[-Q-1+inde-jnde]=='*':
                                            if (inde<MMDetect) and (Q-1+inde-jnde)!=LT:
                                                MM1=inde+0
                                            else:
                                                break
                                        MM1+=1 
                                    elif inde>=MinMatch+MM1-1 and ExtensionOK:
                                        Match=True
                                        break
                                    if MM1>MMDetect:
                                        break
                                else:      
                                    Match=True   ## this statement only runs of we get to the end of the loop (inde=L-1) and a match hasn't been excluded by number of mismatches

                        if Match:
                            if Q>0:
                                QT1=Q-jnde
                            else:
                                QT1=Q+jnde
                            if not(QT1 in QList):
                                if FirstHits_Only and not(BestHits_Only):
                                    QList=[QT1]
                                    SampleDone=True
                                    break
                                if not(UniqueHits_Only) and not(BestHits_Only):
                                    QList.append(QT1)
                                if BestHits_Only and MM1==LeastMM:
                                    QList.append(QT1)
                                if BestHits_Only and MM1<LeastMM:
                                    QList=[QT1]
                                if UniqueHits_Only and not(BestHits_Only):
                                    if len(QList)>0:
                                        QList=[]
                                        SampleDone=True
                                        break
                                    else:
                                        QList=[QT1]
                                if SuperHits_Only:
                                    if MM1<=MaxMisMatch and MM1<=LeastMM-SuperHit:
                                        QList=[QT1]
                                    elif MM1<LeastMM+SuperHit:
                                        QList=[]
                                LeastMM=min(MM1,LeastMM)
                        Tryseed=NextTry[Tryseed]
                        if Tryseed==0:
                            Q=0
                            break
                        Q=FirstTry[Tryseed]
                    if SampleDone: break
                if UniqueHits_Only and len(QList)>1: QList=[]
                if FirstHits_Only and len(QList)>1: QList=QList[:1]
                for Q in QList:
                    LCov=L+0
                    if VirtualSegmentLen>0: LCov=VirtualSegmentLen+0
                    if Q>0:
                        if StartCollapse:
                            if StartCollapseA[Q-1] & 1==1: continue
                            StartCollapseA[Q-1]=StartCollapseA[Q-1] | 1
                        if StartLenCollapse:
                            if L in StartLenCollapseD[Q-1]: continue
                            StartLenCollapseD[Q-1][L]=1
                        M+=Multiplicity1
                        QS = Q-1  ## zero based position in T array of the start of aligment
                        QA=QS+MMDetect ## make sure we are not on a boundary
                        QE = Q+L-2 ## zero based position in T array of T the last base of aligment (so QE-QS=L-1)
                        if StartSense:
                            if TotalBin: Start_SenseT+=Multiplicity1
                            if GeneBin: Start_SenseG[Address_to_GeneNum[QA]]+=Multiplicity1
                            if not(GeneBin) or Address_to_GeneNum[QA]==Address_to_GeneNum[QS]:  
                                if BinBin: Start_SenseB[Address_to_BinNum[QS]]+=Multiplicity1
                                if BaseBin: Start_SenseBa[QS]+=Multiplicity1
                        if StartTotal: 
                            if TotalBin: Start_TotalT+=Multiplicity1
                            if GeneBin: Start_TotalG[Address_to_GeneNum[QA]]+=Multiplicity1
                            if not(GeneBin) or Gene_to_Start[Address_to_GeneNum[QA]]<=QS:  
                                if BinBin: Start_TotalB[Address_to_BinNum[QS]]+=Multiplicity1
                                if BaseBin: Start_TotalBa[QS]+=Multiplicity1
                        if DyadSense:
                            QD=QS+DyadOffset
                            if TotalBin: Dyad_SenseT+=Multiplicity1
                            if QD<LT and (not(GeneBin) or Address_to_GeneNum[QA]==Address_to_GeneNum[QD]):
                                if GeneBin: Dyad_SenseG[Address_to_GeneNum[QD]]+=Multiplicity1
                                if BinBin: Dyad_SenseB[Address_to_BinNum[QD]]+=Multiplicity1
                                if BaseBin: Dyad_SenseBa[QD]+=Multiplicity1
                        if DyadTotal:
                            QD=QS+DyadOffset
                            if TotalBin: Dyad_TotalT+=Multiplicity1
                            if QD<LT and (not(GeneBin) or Address_to_GeneNum[QA]==Address_to_GeneNum[QD]):
                                if GeneBin: Dyad_TotalG[Address_to_GeneNum[QD]]+=Multiplicity1
                                if BinBin: Dyad_TotalB[Address_to_BinNum[QD]]+=Multiplicity1
                                if BaseBin: Dyad_TotalBa[QD]+=Multiplicity1
                        if EndSense: 
                            if TotalBin: End_SenseT+=Multiplicity1
                            if GeneBin: End_SenseG[Address_to_GeneNum[QA]]+=Multiplicity1
                            if QE<LT and (not(GeneBin) or Address_to_GeneNum[QA]==Address_to_GeneNum[QE]):                            
                                if BinBin: End_SenseB[Address_to_BinNum[QE]]+=Multiplicity1
                                if BaseBin: End_SenseBa[QE]+=Multiplicity1
                        if EndTotal: 
                            if TotalBin: End_TotalT+=Multiplicity1
                            if GeneBin: End_TotalG[Address_to_GeneNum[QA]]+=Multiplicity1
                            if QE<LT and (not(GeneBin) or Address_to_GeneNum[QA]==Address_to_GeneNum[QE]):
                                if BinBin: End_TotalB[Address_to_BinNum[QE]]+=Multiplicity1
                                if BaseBin: End_TotalBa[QE]+=Multiplicity1
                        if CoverageSense: 
                            for I in range(LCov):
                                if QS+I>=LT-1:break
                                if TotalBin: Coverage_SenseT+=Multiplicity1
                                if not(GeneBin) or Address_to_GeneNum[QS+I]==Address_to_GeneNum[QA]:
                                    if GeneBin: Coverage_SenseG[Address_to_GeneNum[QS+I]]+=Multiplicity1
                                    if BinBin: Coverage_SenseB[Address_to_BinNum[QS+I]]+=Multiplicity1
                                    if BaseBin: Coverage_SenseBa[QS+I]+=Multiplicity1
                        if CoverageTotal: 
                            for I in range(LCov):
                                if QS+I>=LT-1:break
                                if TotalBin: Coverage_TotalT+=Multiplicity1
                                if not(GeneBin) or Address_to_GeneNum[QS+I]==Address_to_GeneNum[QA]:
                                    if GeneBin: Coverage_TotalG[Address_to_GeneNum[QS+I]]+=Multiplicity1
                                    if BinBin: Coverage_TotalB[Address_to_BinNum[QS+I]]+=Multiplicity1
                                    if BaseBin: Coverage_TotalBa[QS+I]+=Multiplicity1
                        if SizeHistogramSense:
                            La1=L+0
                            if Base1SizeHistogram: La1=(L+0,S[0])
                            if TotalBin:
                                if not (La1 in SizeHistogram_SenseT): SizeHistogram_SenseT[La1]=0
                                SizeHistogram_SenseT[La1]+=Multiplicity1
                            if GeneBin:
                                if not (La1 in SizeHistogram_SenseG[Address_to_GeneNum[QS]]): SizeHistogram_SenseG[Address_to_GeneNum[QS]][La1]=0
                                SizeHistogram_SenseG[Address_to_GeneNum[QA]][La1]+=Multiplicity1
                            if not(GeneBin) or Address_to_GeneNum[QS]==Address_to_GeneNum[QA]:
                                if BinBin:
                                    if not (La1 in SizeHistogram_SenseB[Address_to_BinNum[QS]]): SizeHistogram_SenseB[Address_to_BinNum[QS]][La1]=0
                                    SizeHistogram_SenseB[Address_to_BinNum[QS]][La1]+=Multiplicity1
                                if BaseBin:
                                    if not (La1 in SizeHistogram_SenseBa[QS]): SizeHistogram_SenseBa[QS][La1]=0
                                    SizeHistogram_SenseBa[QS][La1]+=Multiplicity1
                        if SizeHistogramTotal:
                            La1=L+0
                            if Base1SizeHistogram: La1=(L+0,S[0])
                            if TotalBin:
                                if not (La1 in SizeHistogram_TotalT): SizeHistogram_TotalT[La1]=0
                                SizeHistogram_TotalT[La1]+=Multiplicity1
                            if GeneBin:
                                if not (La1 in SizeHistogram_TotalG[Address_to_GeneNum[QS]]): SizeHistogram_TotalG[Address_to_GeneNum[QS]][La1]=0
                                SizeHistogram_TotalG[Address_to_GeneNum[QS]][La1]+=Multiplicity1
                            if not(GeneBin) or Address_to_GeneNum[QS]==Address_to_GeneNum[QA]:
                                if BinBin:
                                    if not (La1 in SizeHistogram_TotalB[Address_to_BinNum[QS]]): SizeHistogram_TotalB[Address_to_BinNum[QS]][La1]=0
                                    SizeHistogram_TotalB[Address_to_BinNum[QS]][La1]+=Multiplicity1
                                if BaseBin:
                                    if not (La1 in SizeHistogram_TotalBa[QS]): SizeHistogram_TotalBa[QS][La1]=0
                                    SizeHistogram_TotalBa[QS][La1]+=Multiplicity1
                        if BaseMatchesSense:
                            ## Note that this routine could be augmented to provide the user with information about the diversity of
                            ## different read starts in identifying a mismatch, whether both strands support a mismatch, and whether quality scores in FastQ support this
                            MMa1=0
                            for I in range(L):
                                if T[QS+I]=='*':break
                                if S[I].upper()!=Tup[I+QS]:
                                    MMa1+=1
                                    if MMa1>MaxMisMatch: break
                                BB='+'+Tup[I+QS]+'>'+S[I].upper()
                                if TotalBin:
                                    if not (BB in BaseMatches_SenseT): BaseMatches_SenseT[BB]=0
                                    BaseMatches_SenseT[BB]+=Multiplicity1
                                if GeneBin:
                                    if not (BB in BaseMatches_SenseG[Address_to_GeneNum[QS+I]]): BaseMatches_SenseG[Address_to_GeneNum[QS+I]][BB]=0
                                    BaseMatches_SenseG[Address_to_GeneNum[QS+I]][BB]+=Multiplicity1
                                if BinBin:
                                    if not (BB in BaseMatches_SenseB[Address_to_BinNum[QS+I]]): BaseMatches_SenseB[Address_to_BinNum[QS+I]][BB]=0
                                    BaseMatches_SenseB[Address_to_BinNum[QS+I]][BB]+=Multiplicity1
                                if BaseBin:
                                    if not (BB in BaseMatches_SenseBa[QS+I]): BaseMatches_SenseBa[QS+I][BB]=0
                                    BaseMatches_SenseBa[QS+I][BB]+=Multiplicity1
                        if BaseMatchesTotal: 
                            MMa1=0
                            for I in range(L):
                                if T[QS+I]=='*':break
                                if S[I].upper()!=Tup[I+QS]:
                                    MMa1+=1
                                    if MMa1>MaxMisMatch: break
                                BB='+'+Tup[I+QS]+'>'+S[I].upper()
                                if TotalBin:
                                    if not (BB in BaseMatches_TotalT): BaseMatches_TotalT[BB]=0
                                    BaseMatches_TotalT[BB]+=Multiplicity1
                                if GeneBin:
                                    if not (BB in BaseMatches_TotalG[Address_to_GeneNum[QS+I]]): BaseMatches_TotalG[Address_to_GeneNum[QS+I]][BB]=0
                                    BaseMatches_TotalG[Address_to_GeneNum[QS+I]][BB]+=Multiplicity1
                                if BinBin:
                                    if not (BB in BaseMatches_TotalB[Address_to_BinNum[QS+I]]): BaseMatches_TotalB[Address_to_BinNum[QS+I]][BB]=0
                                    BaseMatches_TotalB[Address_to_BinNum[QS+I]][BB]+=Multiplicity1
                                if BaseBin:
                                    if not (BB in BaseMatches_TotalBa[QS+I]): BaseMatches_TotalBa[QS+I][BB]=0
                                    BaseMatches_TotalBa[QS+I][BB]+=Multiplicity1
                        if ReadMatchesSense or PileUp:
                            MMX1=0
                            S1=''
                            for I in range(L):
                                if T[I+QS]=='*':
                                    S1+=S[I:].upper()
                                    break
                                if S[I].upper()!=Tup[I+QS]:
                                    MMX1+=1
                                    if MMX1>MaxMisMatch and not(KeepFullQuery): break
                                    S1+=S[I].upper()
                                else:
                                    S1+=S[I].lower()
                            if TotalBin:
                                S2=(S1,QS-1-MaxMisMatch)
                                if not (S2 in ReadMatches_SenseT):
                                    ReadMatches_SenseT[S2]=0
                                    MMSenseOffsetT=max(MMSenseOffsetT,QS+len(S1)-LT+2)
                                ReadMatches_SenseT[S2]+=Multiplicity1
                            if GeneBin:
                                GN1=Address_to_GeneNum[QS]
                                GP1=QS-Gene_to_Start[GN1]
                                S2=(S1,GP1)
                                if not (S2 in ReadMatches_SenseG[GN1]):
                                    ReadMatches_SenseG[GN1][S2]=0
                                    MMSenseOffsetG[GN1]=max(MMSenseOffsetG[GN1],GP1+len(S1)-Gene_to_Len[GN1])
                                ReadMatches_SenseG[GN1][S2]+=Multiplicity1
                            if BinBin:
                                BN1=Address_to_BinNum[QS]
                                BP1=QS-Bin_to_TStart[BN1]
                                S2=(S1,BP1)
                                if not (S2 in ReadMatches_SenseB[BN1]):
                                    ReadMatches_SenseB[BN1][S2]=0
                                    MMSenseOffsetB[BN1]=max(MMSenseOffsetB[BN1],BP1+len(S1)-Bin_to_Len[BN1])
                                ReadMatches_SenseB[BN1][S2]+=Multiplicity1
                            if BaseBin:
                                S2=(S1,QS)
                                if not (S2 in ReadMatches_SenseBa[QS]): ReadMatches_SenseBa[QS][S2]=0
                                ReadMatches_SenseBa[QS][S2]+=Multiplicity1
     
     
                        if CompositionMatrix:
                            COrigin=QS+0
                            if CompositionCenterEnd: COrigin=QS+L                            
                            for I in range(Composition_Start,Composition_End):
                                Qs1=COrigin+I
                                Qe1=COrigin+I+CompositionTupleLen
                                if Qs1<0 or Qs1>=LT or Qe1<0 or Qe1>=LT: continue
                                QK=Tup[Qs1:Qe1]                                
                                if CompositionIndexed:
                                    if TotalBin:
                                        if not (QK,I) in CompositionMatrixDT: CompositionMatrixDT[(QK,I)]=0
                                        CompositionMatrixDT[(QK,I)]+=Multiplicity1
                                    if GeneBin:
                                        if not (QK,I) in CompositionMatrixDG[Address_to_GeneNum[QS]]: CompositionMatrixDG[Address_to_GeneNum[QS]][(QK,I)]=0
                                        CompositionMatrixDG[Address_to_GeneNum[QS]][(QK,I)]+=Multiplicity1
                                    if BinBin:
                                        if not (QK,I) in CompositionMatrixDB[Address_to_BinNum[QS]]: CompositionMatrixDB[Address_to_BinNum[QS]][(QK,I)]=0
                                        CompositionMatrixDB[Address_to_BinNum[QS]][(QK,I)]+=Multiplicity1
                                    if BaseBin: ## can't think of why a base-level composition matrix would be useful, but will leave in this code
                                        if not (QK,I) in CompositionMatrixDBa[QS]: CompositionMatrixDBa[QS][(QK,I)]=0
                                        CompositionMatrixDBa[QS][(QK,I)]+=Multiplicity1
                                else:
                                    if TotalBin:
                                        if not QK in CompositionMatrixDT: CompositionMatrixDT[QK]=0
                                        CompositionMatrixDT[QK]+=Multiplicity1
                                    if GeneBin:
                                        if not QK in CompositionMatrixDG[Address_to_GeneNum[QS]]: CompositionMatrixDG[Address_to_GeneNum[QS]][QK]=0
                                        CompositionMatrixD[Address_to_GeneNum[QS]][QK]+=Multiplicity1
                                    if BinBin:
                                        if not QK in CompositionMatrixDB[Address_to_BinNum[QS]]: CompositionMatrixDB[Address_to_BinNum[QS]][QK]=0
                                        CompositionMatrixD[Address_to_BinNum[QS]][QK]+=Multiplicity1
                                    if BaseBin: ## can't think of why a base-level composition matrix would be useful, but will leave in this code
                                        if not (QK,I) in CompositionMatrixDBa[QS]: CompositionMatrixDBa[QS][QK]=0
                                        CompositionMatrixDBa[QS][QK]+=Multiplicity1
                                    
                        if SimpleOutput:
                            if GeneBin:
                                QR=QS-Gene_to_Start[Address_to_GeneNum[QS]]
                                SO1=GeneNameA[Address_to_GeneNum[QS]]+'\t'
                            else:
                                QR=QS-MaxMisMatch-1
                                SO1=ReferenceFileName+'\t'
                            ## aardvark, output details for bed, psl, pslm
                            if pslformat:
                                ##  A best approximation of a psl format
                                ##    matches - Number of bases that match that aren't repeats
                                ##    misMatches - Number of bases that don't match
                                ##    repMatches - Number of bases that match but are part of repeats
                                ##    nCount - Number of 'N' bases
                                ##    qNumInsert - Number of inserts in query
                                ##    qBaseInsert - Number of bases inserted in query
                                ##    tNumInsert - Number of inserts in target
                                ##    tBaseInsert - Number of bases inserted in target
                                ##    strand - '+' or '-' for query strand. For translated alignments, second '+'or '-' is for genomic strand
                                ##    qName - Query sequence name
                                ##    qSize - Query sequence size
                                ##    qStart - Alignment start position in query
                                ##    qEnd - Alignment end position in query
                                ##    tName - Target sequence name
                                ##    tSize - Target sequence size
                                ##    tStart - Alignment start position in target
                                ##    tEnd - Alignment end position in target
                                ##    blockCount - Number of blocks in the alignment (a block contains no gaps)
                                ##    blockSizes - Comma-separated list of sizes of each block
                                ##    qStarts - Comma-separated list of starting positions of each block in query
                                ##    tStarts - Comma-separated list of starting positions of each block in target
                                tsize1=Gene_to_Len[Address_to_GeneNum[QS]]
                                if MaxMisMatch==0 and not(KeepFullQuery) and not(ExtensionOK):
                                    for mult1 in xrange(Multiplicity1):
                                        SFile.write('%i\t%i\t0\t0\t0\t0\t0\t0\t+\t%s\t%i\t%i\t%i\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i'%(L,0,S+'_+'+S[0]+str(L),L,0,L,SO1,tsize1,QR,QR+L,1,L,0,QR))
                                        if pslAddSeq:
                                            SFile.write('\t%s'%(S))
                                        SFile.write('\r')
                                else:
                                    MMX1=0
                                    S1=''
                                    for I in range(L):
                                        if I+QS<LT and S[I].upper()!=Tup[I+QS]:
                                            MMX1+=1
                                            if MMX1>MaxMisMatch and not(KeepFullQuery): break
                                            S1+=S[I].upper()
                                        else:
                                            S1+=S[I].lower()
                                    L1=len(S1)
                                    for mult1 in xrange(Multiplicity1):
                                        SFile.write('%i\t%i\t0\t0\t0\t0\t0\t0\t+\t%s\t%i\t%i\t%i\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i'%(L1-MMX1,MMX1,S1+'_+'+S1[0]+str(L1),L1,0,L1,SO1,tsize1,QR,QR+L1,1,L1,0,QR))
                                        if pslAddSeq:
                                            SFile.write('\t%s'%(S))
                                        SFile.write('\r')
                                    
                                    
                            if bedformat:
                                ##  Here are the bed fields, the first three being required and the rest optional
                                ##    chrom -
                                ##    chromStart - zero-basd
                                ##    chromEnd -  not included
                                ##    name - Defines the name of the BED line. 
                                ##    score - A score between 0 and 1000. If the track line useScore attribute is set to 1 will determine the level of gray (higher numbers = darker gray). 
                                ##    strand - '+' or '-'.
                                ##    thickStart -  starting position at which  feature is drawn thickly 
                                ##    thickEnd -  ending position at which the feature is drawn thickly 
                                ##    itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On"
                                ##    blockCount - The number of blocks (exons) in the BED line.
                                ##    blockSizes - A comma-separated list of the block sizes. 
                                ##    blockStarts - A comma-separated list of block starts.
                                if VirtualSegmentLen>0:
                                    SFile.write(Multiplicity1*(SO1+'%i\t%i\t%s\t0\t+\t%i\t%i\t%s\r'%(QR,QR+VirtualSegmentLen,'+'+str(QR),QR,QR+L,b0)))
                                elif MaxMisMatch==0 and not(KeepFullQuery) and not(ExtensionOK):
                                    SFile.write(Multiplicity1*(SO1+'%i\t%i\t%s\t0\t+\t%i\t%i\t%s\r'%(QR,QR+L,'+'+S[0]+str(L)+'_'+S,QR,QR+L,b0)))
                                else:
                                    MMX1=0
                                    S1=''
                                    B1=[]
                                    B2=[]
                                    NumLower=0
                                    NumUpper=1
                                    for I in range(L):
                                        if I+QS<LT and S[I].upper()!=Tup[I+QS]:
                                            MMX1+=1
                                            if MMX1>MaxMisMatch and not(KeepFullQuery): break
                                            S1+=S[I].upper()
                                            if NumLower>0:
                                                B2.append(NumLower)
                                                NumLower=0
                                            NumUpper+=1
                                        else:
                                            S1+=S[I].lower()
                                            NumLower+=1
                                            if NumUpper>0:
                                                B1.append(QR+I)
                                                NumUpper=0
                                    if NumLower>0:
                                        B2.append(NumLower)
                                    L1=len(S1)
                                    if MMX1==0:
                                        SFile.write(Multiplicity1*(SO1+'%i\t%i\t%s\t0\t+\t%i\t%i\t%s\t%i\t%s\t%s\r'%(QR,QR+L1,'+'+str(QR)+'_'+S1,B1[0],B1[-1]+B2[-1],b0,len(B1),str(B2)[1:-1],str(B1)[1:-1])))
                                    elif MMX1==1:
                                        SFile.write(Multiplicity1*(SO1+'%i\t%i\t%s\t0\t+\t%i\t%i\t%s\t%i\t%s\t%s\r'%(QR,QR+L1,'+'+str(QR)+'_'+S1,B1[0],B1[-1]+B2[-1],b1,len(B1),str(B2)[1:-1],str(B1)[1:-1])))
                                    else:
                                        SFile.write(Multiplicity1*(SO1+'%i\t%i\t%s\t0\t+\t%i\t%i\t%s\t%i\t%s\t%s\r'%(QR,QR+L1,'+'+str(QR)+'_'+S1,B1[0],B1[-1]+B2[-1],b2,len(B1),str(B2)[1:-1],str(B1)[1:-1])))

                            if hqaformat:
                                MMX1=0
                                S1=''
                                for I in range(L):
                                    if I+QS<LT and S[I].upper()!=Tup[I+QS]:
                                        MMX1+=1
                                        if MMX1>MaxMisMatch and not(KeepFullQuery): break
                                        S1+=S[I].upper()
                                    else:
                                        S1+=S[I].lower()
                                SO1+=str(QR+1)+'\t'
                                SO1+=str(QR+len(S1))+'\t'
                                SO1+="'+'"+'\t'
                                if KeepFullQuery:
                                    SO1+=str(len(S1)-min(MaxMisMatch,MMX1))+'\t'
                                else:
                                    SO1+=str(len(S1)-MMX1)+'\t'
                                SO1+=str(len(S1))+'\t'
                                SO1+=S1+'\t'
                                SO1+=T[QS:QS+len(S1)]+'\t'
                                if GeneBin:
                                    SO1+=str(Address_to_GeneNum[QS]+1)+'\t'
                                else:
                                    SO1+=str(1)+'\t'                            
                                SO1+=Experiment+'\r'
                                if SimpleOutputSum:
                                    if not SO1 in SimpleOutputD:
                                        SimpleOutputD[SO1]=0
                                    SimpleOutputD[SO1]+=1
                                else:
                                    SFile.write(Multiplicity1*SO1)


                    if Q<0:
                        if StartCollapse:
                            if StartCollapseA[-Q-1] & 2==2: continue
                            StartCollapseA[-Q-1]=StartCollapseA[-Q-1] | 2
                        if StartLenCollapse:
                            if -L in StartLenCollapseD[-Q-1]: continue
                            StartLenCollapseD[-Q-1][-L]=1

                        N+=Multiplicity1
                        LCov=L+0
                        if VirtualSegmentLen>0: LCov=VirtualSegmentLen+0
                        QS = LT+Q   ## zero based position in T array of start of alignment (3' relative to T array)
                        QE = LT+Q-L+1  ## zero based position in T array of end of alignment (5' relative to T array).  So QE-QS=L-1

                        if StartAntiSense:
                            if TotalBin: Start_AntiSenseT+=Multiplicity1
                            if GeneBin: Start_AntiSenseG[Address_to_GeneNum[QS]]+=Multiplicity1
                            if BinBin: Start_AntiSenseB[Address_to_BinNum[QS]]+=Multiplicity1
                            if BaseBin: Start_AntiSenseBa[QS]+=Multiplicity1
                        if DyadAntiSense:
                            if TotalBin: Dyad_AntiSenseT+=Multiplicity1
                            QD=QS-DyadOffset
                            if QD>=0 and (not(GeneBin) or Address_to_GeneNum[QS]==Address_to_GeneNum[QD]):  ## don't count dyads outside of the current gene
                                if GeneBin: Dyad_AntiSenseG[Address_to_GeneNum[QD]]+=Multiplicity1
                                if BinBin: Dyad_AntiSenseB[Address_to_BinNum[QD]]+=Multiplicity1
                                if BaseBin: Dyad_AntiSenseBa[QD]+=Multiplicity1
                        if EndAntiSense: 
                            if TotalBin: End_AntiSenseT+=Multiplicity1
                            if QE>=0:
                                if GeneBin: End_AntiSenseG[Address_to_GeneNum[QE]]+=Multiplicity1
                                if BinBin: End_AntiSenseB[Address_to_BinNum[QE]]+=Multiplicity1
                                if BaseBin: End_AntiSenseBa[QE]+=Multiplicity1
                        if CoverageAntiSense: 
                            for I in range(LCov):
                                if QS-I<0:break
                                if TotalBin: Coverage_AntiSenseT+=Multiplicity1
                                if not(GeneBin) or (Address_to_GeneNum[QS]==Address_to_GeneNum[QS-I]):
                                    if GeneBin: Coverage_AntiSenseG[Address_to_GeneNum[QS-I]]+=Multiplicity1
                                    if BinBin: Coverage_AntiSenseB[Address_to_BinNum[QS-I]]+=Multiplicity1
                                    if BaseBin: Coverage_AntiSenseBa[QS-I]+=Multiplicity1
                        if SizeHistogramAntiSense:
                            La1=L+0
                            if Base1SizeHistogram: La1=(L+0,S[0])
                            if TotalBin:
                                if not (La1 in SizeHistogram_AntiSenseT): SizeHistogram_AntiSenseT[La1]=0
                                SizeHistogram_AntiSenseT[La1]+=Multiplicity1
                            if GeneBin:
                                if not (La1 in SizeHistogram_AntiSenseG[Address_to_GeneNum[QS]]): SizeHistogram_AntiSenseG[Address_to_GeneNum[QS]][La1]=0
                                SizeHistogram_AntiSenseG[Address_to_GeneNum[QS]][La1]+=Multiplicity1
                            if BinBin:
                                if not (La1 in SizeHistogram_AntiSenseB[Address_to_BinNum[QS]]): SizeHistogram_AntiSenseB[Address_to_BinNum[QS]][La1]=0
                                SizeHistogram_AntiSenseB[Address_to_BinNum[QS]][La1]+=Multiplicity1
                            if BaseBin:
                                if not (La1 in SizeHistogram_AntiSenseBa[QS]): SizeHistogram_AntiSenseBa[QS][La1]=0
                                SizeHistogram_AntiSenseBa[QS][La1]+=Multiplicity1
                        if BaseMatchesAntiSense: 
                            MMa1=0
                            for I in range(L):
                                if T[QS-I]=='*':break
                                if S[I].upper()!=Uup[I+LT-QS-1]:
                                    MMa1+=1
                                    if MMa1>MaxMisMatch: break
                                BB='-'+Uup[I+LT-QS-1]+'>'+S[I].upper()
                                if TotalBin:
                                    if not (BB in BaseMatches_AntiSenseT): BaseMatches_AntiSenseT[BB]=0
                                    BaseMatches_AntiSenseT[BB]+=Multiplicity1
                                if GeneBin:
                                    if not (BB in BaseMatches_AntiSenseG[Address_to_GeneNum[QS-I]]): BaseMatches_AntiSenseG[Address_to_GeneNum[QS-I]][BB]=0
                                    BaseMatches_AntiSenseG[Address_to_GeneNum[QS-I]][BB]+=Multiplicity1
                                if BinBin:
                                    if not (BB in BaseMatches_AntiSenseB[Address_to_BinNum[QS-I]]): BaseMatches_AntiSenseB[Address_to_BinNum[QS-I]][BB]=0
                                    BaseMatches_AntiSenseB[Address_to_BinNum[QS-I]][BB]+=Multiplicity1
                                if BaseBin:
                                    if not (BB in BaseMatches_AntiSenseBa[QS-I]): BaseMatches_AntiSenseBa[QS-I][BB]=0
                                    BaseMatches_AntiSenseBa[QS-I][BB]+=Multiplicity1

                        if StartTotal:
                            if TotalBin: Start_TotalT+=Multiplicity1
                            if GeneBin: Start_TotalG[Address_to_GeneNum[QS]]+=Multiplicity1
                            if BinBin: Start_TotalB[Address_to_BinNum[QS]]+=Multiplicity1
                            if BaseBin: Start_TotalBa[QS]+=Multiplicity1
                        if DyadTotal:
                            if TotalBin: Dyad_TotalT+=Multiplicity1
                            QD=QS-DyadOffset
                            if QD>=0 and (not(GeneBin) or Address_to_GeneNum[QS]==Address_to_GeneNum[QD]):  ## don't count dyads outside of the current gene
                                if GeneBin: Dyad_TotalG[Address_to_GeneNum[QD]]+=Multiplicity1
                                if BinBin: Dyad_TotalB[Address_to_BinNum[QD]]+=Multiplicity1
                                if BaseBin: Dyad_TotalBa[QD]+=Multiplicity1
                        if EndTotal: 
                            if TotalBin: End_TotalT+=Multiplicity1
                            if QE>=0:
                                if GeneBin: End_TotalG[Address_to_GeneNum[QE]]+=Multiplicity1
                                if BinBin: End_TotalB[Address_to_BinNum[QE]]+=Multiplicity1
                                if BaseBin: End_TotalBa[QE]+=Multiplicity1
                        if CoverageTotal: 
                            for I in range(LCov):
                                if QS-I<0:break
                                if TotalBin: Coverage_TotalT+=Multiplicity1
                                if (not(GeneBin) or Address_to_GeneNum[QS]==Address_to_GeneNum[QS-I]):
                                    if GeneBin: Coverage_TotalG[Address_to_GeneNum[QS-I]]+=Multiplicity1
                                    if BinBin: Coverage_TotalB[Address_to_BinNum[QS-I]]+=Multiplicity1
                                    if BaseBin: Coverage_TotalBa[QS-I]+=Multiplicity1
                        if SizeHistogramTotal:
                            La1=L+0
                            if Base1SizeHistogram: La1=(L+0,S[0])
                            if TotalBin:
                                if not (La1 in SizeHistogram_TotalT): SizeHistogram_TotalT[La1]=0
                                SizeHistogram_TotalT[La1]+=Multiplicity1
                            if GeneBin:
                                if not (La1 in SizeHistogram_TotalG[Address_to_GeneNum[QS]]): SizeHistogram_TotalG[Address_to_GeneNum[QS]][La1]=0
                                SizeHistogram_TotalG[Address_to_GeneNum[QS]][La1]+=Multiplicity1
                            if BinBin:
                                if not (La1 in SizeHistogram_TotalB[Address_to_BinNum[QS]]): SizeHistogram_TotalB[Address_to_BinNum[QS]][La1]=0
                                SizeHistogram_TotalB[Address_to_BinNum[QS]][La1]+=Multiplicity1
                            if BaseBin:
                                if not (La1 in SizeHistogram_TotalBa[QS]): SizeHistogram_TotalBa[QS][La1]=0
                                SizeHistogram_TotalBa[QS][La1]+=Multiplicity1

                        if BaseMatchesTotal: 
                            MMa1=0
                            for I in range(L):
                                if T[QS-I]=='*':break
                                if S[I].upper()!=Uup[I+LT-QS-1]:
                                    MMa1+=1
                                    if MMa1>MaxMisMatch: break
                                BB='-'+Uup[I+LT-QS-1]+'>'+S[I].upper()
                                if TotalBin:
                                    if not (BB in BaseMatches_TotalT): BaseMatches_TotalT[BB]=0
                                    BaseMatches_TotalT[BB]+=Multiplicity1
                                if GeneBin:
                                    if not (BB in BaseMatches_TotalG[Address_to_GeneNum[QS-I]]): BaseMatches_TotalG[Address_to_GeneNum[QS-I]][BB]=0
                                    BaseMatches_TotalG[Address_to_GeneNum[QS-I]][BB]+=Multiplicity1
                                if BinBin:
                                    if not (BB in BaseMatches_TotalB[Address_to_BinNum[QS-I]]): BaseMatches_TotalB[Address_to_BinNum[QS-I]][BB]=0
                                    BaseMatches_TotalB[Address_to_BinNum[QS-I]][BB]+=Multiplicity1
                                if BaseBin:
                                    if not (BB in BaseMatches_TotalBa[QS-I]): BaseMatches_TotalBa[QS-I][BB]=0
                                    BaseMatches_TotalBa[QS-I][BB]+=Multiplicity1
                        if ReadMatchesAntiSense or PileUp:
                            MMX1=0
                            S1=''
                            for I in range(L):
                                if QS-I==0 or Tup[QS-I]=='*':
                                    S1+=S[I:]  ## this may be longer than can be displayed
                                    break
                                if S[I].upper()!=Uup[LT-1-QS+I]:
                                    MMX1+=1
                                    if MMX1>MaxMisMatch and not(KeepFullQuery): break
                                    S1+=S[I].upper()
                                else:
                                    S1=S1+S[I].lower()
                                LS1=len(S1)
                            if TotalBin:
                                S2=(S1,QS-1-MaxMisMatch)
                                if not (S2 in ReadMatches_AntiSenseT):
                                    ReadMatches_AntiSenseT[S2]=0
                                    MMAntiSenseOffsetT=max(MMAntiSenseOffsetT,len(S1)-QS)
                                ReadMatches_AntiSenseT[S2]+=Multiplicity1
                            if GeneBin:
                                GN1=Address_to_GeneNum[QS]
                                GP1=QS-Gene_to_Start[GN1]
                                S2=(S1,GP1)
                                if not (S2 in ReadMatches_AntiSenseG[GN1]):
                                    ReadMatches_AntiSenseG[GN1][S2]=0
                                    MMAntiSenseOffsetG[GN1]=max(MMAntiSenseOffsetG[GN1],len(S1)-MaxMisMatch-GP1-1)
                                ReadMatches_AntiSenseG[GN1][S2]+=Multiplicity1
                            if BinBin:
                                BN1=Address_to_BinNum[QS]
                                BP1=QS-Bin_to_TStart[BN1]
                                S2=(S1,BP1)
                                if not (S2 in ReadMatches_AntiSenseB[BN1]):
                                    ReadMatches_AntiSenseB[BN1][S2]=0
                                    MMAntiSenseOffsetB[BN1]=max(MMAntiSenseOffsetB[BN1],len(S1)-BP1-1)
                                ReadMatches_AntiSenseB[BN1][S2]+=Multiplicity1
                            if BaseBin:
                                S2=(S1,0)
                                if not (S2 in ReadMatches_AntiSenseBa[QS]): ReadMatches_AntiSenseBa[QS][S2]=0
                                ReadMatches_AntiSenseBa[QS][S2]+=Multiplicity1

                        if CompositionMatrix:
                            COrigin=QS
                            if CompositionCenterEnd: COrigin=QE                            
                            for I in range(Composition_Start,Composition_End):
                                Qs1=COrigin-I
                                Qe1=COrigin-I-CompositionTupleLen
                                if Qs1<0 or Qs1>=LT or Qe1<0 or Qe1>=LT: continue
                                QK=Uup[LT-1-Qs1:LT-1-Qe1]                                
                                if CompositionIndexed:
                                    if TotalBin:
                                        if not (QK,I) in CompositionMatrixDT: CompositionMatrixDT[(QK,I)]=0
                                        CompositionMatrixDT[(QK,I)]+=Multiplicity1
                                    if GeneBin:
                                        if not (QK,I) in CompositionMatrixDG[Address_to_GeneNum[QS]]: CompositionMatrixDG[Address_to_GeneNum[QS]][(QK,I)]=0
                                        CompositionMatrixDG[Address_to_GeneNum[QS]][(QK,I)]+=Multiplicity1
                                    if BinBin:
                                        if not (QK,I) in CompositionMatrixDB[Address_to_BinNum[QS]]: CompositionMatrixDB[Address_to_BinNum[QS]][(QK,I)]=0
                                        CompositionMatrixDB[Address_to_BinNum[QS]][(QK,I)]+=Multiplicity1
                                    if BaseBin:
                                        if not (QK,I) in CompositionMatrixDBa[QS]: CompositionMatrixDBa[QS][(QK,I)]=0
                                        CompositionMatrixDBa[QS][(QK,I)]+=Multiplicity1
                                else:
                                    if TotalBin:
                                        if not QK in CompositionMatrixDT: CompositionMatrixDT[QK]=0
                                        CompositionMatrixDT[QK]+=Multiplicity1
                                    if GeneBin:
                                        if not QK in CompositionMatrixDG[Address_to_GeneNum[QS]]: CompositionMatrixDG[Address_to_GeneNum[QS]][QK]=0
                                        CompositionMatrixDG[Address_to_GeneNum[QS]][QK]+=Multiplicity1
                                    if BinBin:
                                        if not QK in CompositionMatrixDG[Address_to_BinNum[QS]]: CompositionMatrixDG[Address_to_BinNum[QS]][QK]=0
                                        CompositionMatrixDG[Address_to_BinNum[QS]][QK]+=Multiplicity1
                                    if BaseBin:
                                        if not (QK,I) in CompositionMatrixDBa[QS]: CompositionMatrixDBa[QS][QK]=0
                                        CompositionMatrixDBa[QS][QK]+=Multiplicity1
                        if SimpleOutput:
                            ## aardvark, output details for bed, psl, bedgraph
                            if GeneBin:
                                QR=QS-Gene_to_Start[Address_to_GeneNum[QS]]
                                SO1=GeneNameA[Address_to_GeneNum[QS]]+'\t'
                            else:
                                QR=QS-MaxMisMatch-1
                                SO1=ReferenceFileName+'\t'
                            if pslformat:
                                ##  A best approximation of a psl format
                                ##    matches - Number of bases that match that aren't repeats
                                ##    misMatches - Number of bases that don't match
                                ##    repMatches - Number of bases that match but are part of repeats
                                ##    nCount - Number of 'N' bases
                                ##    qNumInsert - Number of inserts in query
                                ##    qBaseInsert - Number of bases inserted in query
                                ##    tNumInsert - Number of inserts in target
                                ##    tBaseInsert - Number of bases inserted in target
                                ##    strand - '+' or '-' for query strand. For translated alignments, second '+'or '-' is for genomic strand
                                ##    qName - Query sequence name
                                ##    qSize - Query sequence size
                                ##    qStart - Alignment start position in query
                                ##    qEnd - Alignment end position in query
                                ##    tName - Target sequence name
                                ##    tSize - Target sequence size
                                ##    tStart - Alignment start position in target
                                ##    tEnd - Alignment end position in target
                                ##    blockCount - Number of blocks in the alignment (a block contains no gaps)
                                ##    blockSizes - Comma-separated list of sizes of each block
                                ##    qStarts - Comma-separated list of starting positions of each block in query
                                ##    tStarts - Comma-separated list of starting positions of each block in target
                                tsize1=Gene_to_Len[Address_to_GeneNum[QS]]
                                if MaxMisMatch==0 and not(KeepFullQuery) and not(ExtensionOK):
                                    for mult1 in xrange(Multiplicity1):
                                        SFile.write('%i\t%i\t0\t0\t0\t0\t0\t0\t-\t%s\t%i\t%i\t%i\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i'%(L,0,S+'_-'+S[0]+str(L),L,0,L,SO1,tsize1,QR-L+1,QR+1,1,L,0,QR-L+1))
                                        if pslAddSeq:
                                            SFile.write('\t%s'%(S))
                                        SFile.write('\r')
                                else:
                                    MMX1=0
                                    S1=''
                                    for I in range(L):
                                        if I+QS<LT and S[I].upper()!=Tup[I+QS]:
                                            MMX1+=1
                                            if MMX1>MaxMisMatch and not(KeepFullQuery): break
                                            S1+=S[I].upper()
                                        else:
                                            S1+=S[I].lower()
                                    L1=len(S1)
                                    for mult1 in xrange(Multiplicity1):
                                        SFile.write('%i\t%i\t0\t0\t0\t0\t0\t0\t-\t%s\t%i\t%i\t%i\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i'%(L1-MMX1,MMX1,S1+'_-'+S1[0]+str(L1),L1,0,L1,SO1,tsize1,QR-L1+1,QR+1,1,L1,0,QR-L1+1))
                                        if pslAddSeq:
                                            SFile.write('\t%s'%(S))
                                        SFile.write('\r')                                    
                                    
                            if bedformat:
                                ##  Here are the bed fields, the first three being required and the rest optional
                                ##    chrom -
                                ##    chromStart - zero-basd
                                ##    chromEnd -  not included
                                ##    name - Defines the name of the BED line. 
                                ##    score - A score between 0 and 1000. If the track line useScore attribute is set to 1 will determine the level of gray (higher numbers = darker gray). 
                                ##    strand - '+' or '-'.
                                ##    thickStart -  starting position at which  feature is drawn thickly 
                                ##    thickEnd -  ending position at which the feature is drawn thickly 
                                ##    itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On"
                                ##    blockCount - The number of blocks (exons) in the BED line.
                                ##    blockSizes - A comma-separated list of the block sizes. 
                                ##    blockStarts - A comma-separated list of block starts.
                                if VirtualSegmentLen>0:
                                    SFile.write(Multiplicity1*(SO1+'%i\t%i\t%s\t0\t-\t%i\t%i\t%s\r'%(QR-VirtualSegmentLen+1,QR+1,'-'+str(QR),QR-L+1,QR+1,r0)))
                                elif MaxMisMatch==0 and not(KeepFullQuery) and not(ExtensionOK):
                                    SFile.write(Multiplicity1*(SO1+'%i\t%i\t%s\t0\t-\t%i\t%i\t%s\r'%(QR-L+1,QR+1,'-'+S[0]+str(L)+'_'+S,QR-L+1,QR+1,r0)))
                                else:
                                    MMX1=0
                                    S1=''
                                    B1=[]
                                    B2=[]
                                    NumLower=0
                                    NumUpper=1
                                    for I in range(L):
                                        if QS-I>0 and S[I].upper()!=Uup[LT-1-QS+I]:
                                            MMX1+=1
                                            if MMX1>MaxMisMatch and not(KeepFullQuery): break
                                            S1+=S[I].upper()
                                            if NumLower>0:
                                                B2.append(NumLower)
                                                B1.append(QR-I+1)
                                                NumLower=0
                                            NumUpper+=1
                                        else:
                                            S1+=S[I].lower()
                                            NumLower+=1
                                            if NumUpper>0:
                                                NumUpper=0
                                    if NumLower>0:
                                        B2.append(NumLower)
                                        B1.append(QR-I)
                                    L1=len(S1)
                                    B2=B2[::-1]
                                    B1=B1[::-1]
                                    if MMX1==0:
                                        SFile.write(Multiplicity1*(SO1+'%i\t%i\t%s\t0\t-\t%i\t%i\t%s\t%i\t%s\t%s\r'%(QR-L1+1,QR+1,'-'+str(QR)+'_'+S1,B1[0],B1[-1]+B2[-1],r0,len(B1),str(B2)[1:-1],str(B1)[1:-1])))
                                    elif MMX1==1:
                                        SFile.write(Multiplicity1*(SO1+'%i\t%i\t%s\t0\t-\t%i\t%i\t%s\t%i\t%s\t%s\r'%(QR-L1+1,QR+1,'-'+str(QR)+'_'+S1,B1[0],B1[-1]+B2[-1],r1,len(B1),str(B2)[1:-1],str(B1)[1:-1])))
                                    else:
                                        SFile.write(Multiplicity1*(SO1+'%i\t%i\t%s\t0\t-\t%i\t%i\t%s\t%i\t%s\t%s\r'%(QR-L1+1,QR+1,'-'+str(QR)+'_'+S1,B1[0],B1[-1]+B2[-1],r2,len(B1),str(B2)[1:-1],str(B1)[1:-1])))

                            if hqaformat:
                                MMX1=0
                                S1=''
                                for I in range(L):
                                    if QS-I>0 and S[I].upper()!=Uup[LT-1-QS+I]:
                                        MMX1+=1
                                        if MMX1>MaxMisMatch and not(KeepFullQuery): break
                                        S1+=S[I].upper()
                                    else:
                                        S1+=S[I].lower()

                                SO1+=str(QR+1)+'\t'
                                SO1+=str(QR-len(S1)+2)+'\t'
                                SO1+="'-'"+'\t'
                                if KeepFullQuery:
                                    SO1+=str(len(S1)-min(MaxMisMatch,MMX1))+'\t'
                                else:
                                    SO1+=str(len(S1)-MMX1)+'\t'
                                SO1+=str(len(S1))+'\t'
                                SO1+=S1+'\t'
                                SO1+=U[LT-1-QS:LT-1-QS+len(S1)]+'\t'
                                if GeneBin:
                                    SO1+=str(Address_to_GeneNum[QS]+1)+'\t'
                                else:
                                    SO1+=str(1)+'\t'                            
                                SO1+=Experiment+'\r'
                                if SimpleOutputSum:
                                    if not SO1 in SimpleOutputD:
                                        SimpleOutputD[SO1]=0
                                    SimpleOutputD[SO1]+=1
                                else:
                                    SFile.write(Multiplicity1*SO1)

                                    
                                    


                P1+=Multiplicity1
                P2+=Multiplicity1
                if P2>=ReportInterval:
                    Report0= "sense match="+str(M)+"; antisense match="+str(N)+"; total >="+str(MinRead)+" with barcode, linker="+str(P1)
                    Qprint(time0,LogFileName,Report0)
                    P2=0
                    if StopAfterQuery>0 and P1>=StopAfterQuery:break   

        
        ReadFile.close()
        if P2 and P2>0:
            Report0= "sense match="+str(M)+"; antisense match="+str(N)+"; total >="+str(MinRead)+" with barcode, linker="+str(P1)
            Qprint(time0,LogFileName,Report0)

        if SimpleOutput and SimpleOutputSum:
            SOL1=list(zip(SimpleOutputD.keys(),SimpleOutputD.values()))
            SOL1=sorted(SOL1,key=lambda x:x[1])
            SOL1=sorted(SOL1,key=lambda x:int(x[0].split('\t')[4]))
            SOL1=sorted(SOL1,key=lambda x:int(x[0].split('\t')[2]))
            SOL1=sorted(SOL1,key=lambda x:int(x[0].split('\t')[1]))
            SOL1=sorted(SOL1,key=lambda x:x[0].split('\t')[3])
            SOL1=sorted(SOL1,key=lambda x:x[0].split('\t')[8])
            SOL1=sorted(SOL1,key=lambda x:x[0].split('\t')[9])
            for SOL11 in SOL1:
                SFile.write(str(SOL11[1])+'\t'+SOL11[0])
            
        if TabularOutput:

            if StartSense:
                outfile.write('Start_Sense'+' '+Experiment)
                if TotalBin:
                    outfile.write('\t'+str(Start_SenseT))
                if GeneBin:
                    for X in Start_SenseG: outfile.write('\t'+str(X))                    
                if BinBin:
                    for X in Start_SenseB: outfile.write('\t'+str(X))                    
                if BaseBin:
                    for (N,X) in enumerate(Start_SenseBa):
                        if T[N]!='*':
                            outfile.write('\t'+str(X))                    
                outfile.write('\r')

            if StartAntiSense:
                outfile.write('Start_AntiSense'+' '+Experiment)
                if TotalBin:
                    outfile.write('\t'+str(Start_AntiSenseT))
                if GeneBin:
                    for X in Start_AntiSenseG: outfile.write('\t'+str(X))                    
                if BinBin:
                    for X in Start_AntiSenseB: outfile.write('\t'+str(X))                    
                if BaseBin:
                    for (N,X) in enumerate(Start_AntiSenseBa):
                        if T[N]!='*':
                            outfile.write('\t'+str(X))                    
                outfile.write('\r')

            if StartTotal:
                outfile.write('Start_Total'+' '+Experiment)
                if TotalBin:
                    outfile.write('\t'+str(Start_TotalT))
                if GeneBin:
                    for X in Start_TotalG: outfile.write('\t'+str(X))                    
                if BinBin:
                    for X in Start_TotalB: outfile.write('\t'+str(X))                    
                if BaseBin:
                    for (N,X) in enumerate(Start_TotalBa):
                        if T[N]!='*':
                            outfile.write('\t'+str(X))                    
                outfile.write('\r')

            if EndSense:
                outfile.write('End_Sense'+' '+Experiment)
                if TotalBin:
                    outfile.write('\t'+str(End_SenseT))
                if GeneBin:
                    for X in End_SenseG: outfile.write('\t'+str(X))                    
                if BinBin:
                    for X in End_SenseB: outfile.write('\t'+str(X))                    
                if BaseBin:
                    for (N,X) in enumerate(End_SenseBa):
                        if T[N]!='*':
                            outfile.write('\t'+str(X))                    
                outfile.write('\r')

            if EndAntiSense:
                outfile.write('End_AntiSense'+' '+Experiment)
                if TotalBin:
                    outfile.write('\t'+str(End_AntiSenseT))
                if GeneBin:
                    for X in End_AntiSenseG: outfile.write('\t'+str(X))                    
                if BinBin:
                    for X in End_AntiSenseB: outfile.write('\t'+str(X))                    
                if BaseBin:
                    for (N,X) in enumerate(End_AntiSenseBa):
                        if T[N]!='*':
                            outfile.write('\t'+str(X))                    
                outfile.write('\r')

            if EndTotal:
                outfile.write('End_Total'+' '+Experiment)
                if TotalBin:
                    outfile.write('\t'+str(End_TotalT))
                if GeneBin:
                    for X in End_TotalG: outfile.write('\t'+str(X))                    
                if BinBin:
                    for X in End_TotalB: outfile.write('\t'+str(X))                    
                if BaseBin:
                    for (N,X) in enumerate(End_TotalBa):
                        if T[N]!='*':
                            outfile.write('\t'+str(X))                    
                outfile.write('\r')

            if DyadSense:
                outfile.write('Dyad_Sense'+' '+Experiment)
                if TotalBin:
                    outfile.write('\t'+str(Dyad_SenseT))
                if GeneBin:
                    for X in Dyad_SenseG: outfile.write('\t'+str(X))                    
                if BinBin:
                    for X in Dyad_SenseB: outfile.write('\t'+str(X))                    
                if BaseBin:
                    for (N,X) in enumerate(Dyad_SenseBa):
                        if T[N]!='*':
                            outfile.write('\t'+str(X))                    
                outfile.write('\r')

            if DyadAntiSense:
                outfile.write('Dyad_AntiSense'+' '+Experiment)
                if TotalBin:
                    outfile.write('\t'+str(Dyad_AntiSenseT))
                if GeneBin:
                    for X in Dyad_AntiSenseG: outfile.write('\t'+str(X))                    
                if BinBin:
                    for X in Dyad_AntiSenseB: outfile.write('\t'+str(X))                    
                if BaseBin:
                    for (N,X) in enumerate(Dyad_AntiSenseBa):
                        if T[N]!='*':
                            outfile.write('\t'+str(X))                    
                outfile.write('\r')

            if DyadTotal:
                outfile.write('Dyad_Total'+' '+Experiment)
                if TotalBin:
                    outfile.write('\t'+str(Dyad_TotalT))
                if GeneBin:
                    for X in Dyad_TotalG: outfile.write('\t'+str(X))                    
                if BinBin:
                    for X in Dyad_TotalB: outfile.write('\t'+str(X))                    
                if BaseBin:
                    for (N,X) in enumerate(Dyad_TotalBa):
                        if T[N]!='*':
                            outfile.write('\t'+str(X))                    
                outfile.write('\r')

            if CoverageSense:
                outfile.write('Coverage_Sense'+' '+Experiment)
                if TotalBin:
                    outfile.write('\t'+str(Coverage_SenseT))
                if GeneBin:
                    for X in Coverage_SenseG: outfile.write('\t'+str(X))                    
                if BinBin:
                    for X in Coverage_SenseB: outfile.write('\t'+str(X))                    
                if BaseBin:
                    for (N,X) in enumerate(Coverage_SenseBa):
                        if T[N]!='*':
                            outfile.write('\t'+str(X))                    
                outfile.write('\r')

            if CoverageAntiSense:
                outfile.write('Coverage_AntiSense'+' '+Experiment)
                if TotalBin:
                    outfile.write('\t'+str(Coverage_AntiSenseT))
                if GeneBin:
                    for X in Coverage_AntiSenseG: outfile.write('\t'+str(X))                    
                if BinBin:
                    for X in Coverage_AntiSenseB: outfile.write('\t'+str(X))                    
                if BaseBin:
                    for (N,X) in enumerate(Coverage_AntiSenseBa):
                        if T[N]!='*':
                            outfile.write('\t'+str(X))                    
                outfile.write('\r')

            if CoverageTotal:
                outfile.write('Coverage_Total'+' '+Experiment)
                if TotalBin:
                    outfile.write('\t'+str(Coverage_TotalT))
                if GeneBin:
                    for X in Coverage_TotalG: outfile.write('\t'+str(X))                    
                if BinBin:
                    for X in Coverage_TotalB: outfile.write('\t'+str(X))                    
                if BaseBin:
                    for (N,X) in enumerate(Coverage_TotalBa):
                        if T[N]!='*':
                            outfile.write('\t'+str(X))                    
                outfile.write('\r')

            if SizeHistogramSense:
                outfile.write('SizeHistogram_Sense'+' '+Experiment)
                if TotalBin:
                    outfile.write('\t'+sdi1(SizeHistogram_SenseT))
                if GeneBin:
                    for X in SizeHistogram_SenseG: outfile.write('\t'+sdi1(X))                    
                if BinBin:
                    for X in SizeHistogram_SenseB: outfile.write('\t'+sdi1(X))                    
                if BaseBin:
                    for (N,X) in enumerate(SizeHistogram_SenseBa):
                        if T[N]!='*':
                            outfile.write('\t'+sdi1(X))                    
                outfile.write('\r')

            if SizeHistogramAntiSense:
                outfile.write('SizeHistogram_AntiSense'+' '+Experiment)
                if TotalBin:
                    outfile.write('\t'+sdi1(SizeHistogram_AntiSenseT))
                if GeneBin:
                    for X in SizeHistogram_AntiSenseG: outfile.write('\t'+sdi1(X))                    
                if BinBin:
                    for X in SizeHistogram_AntiSenseB: outfile.write('\t'+sdi1(X))                    
                if BaseBin:
                    for (N,X) in enumerate(SizeHistogram_AntiSenseBa):
                            if T[N]!='*':
                                outfile.write('\t'+sdi1(X))
                outfile.write('\r')

            if SizeHistogramTotal:
                outfile.write('SizeHistogram_Total'+' '+Experiment)
                if TotalBin:
                    outfile.write('\t'+sdi1(SizeHistogram_TotalT))
                if GeneBin:
                    for X in SizeHistogram_TotalG: outfile.write('\t'+sdi1(X))                    
                if BinBin:
                    for X in SizeHistogram_TotalB: outfile.write('\t'+sdi1(X))                    
                if BaseBin:
                    for (N,X) in enumerate(SizeHistogram_TotalBa):
                        if T[N]!='*':
                            outfile.write('\t'+sdi1(X))                    
                outfile.write('\r')

            if BaseMatchesSense:
                outfile.write('BaseMatches_Sense'+' '+Experiment)
                if TotalBin:
                    outfile.write('\t'+sdi1(BaseMatches_SenseT))
                if GeneBin:
                    for X in BaseMatches_SenseG: outfile.write('\t'+sdi1(X))                    
                if BinBin:
                    for X in BaseMatches_SenseB: outfile.write('\t'+sdi1(X))                    
                if BaseBin:
                    for (N,X) in enumerate(BaseMatches_SenseBa):
                        if T[N]!='*':
                            outfile.write('\t'+sdi1(X))                    
                outfile.write('\r')

            if BaseMatchesAntiSense:
                outfile.write('BaseMatches_AntiSense'+' '+Experiment)
                if TotalBin:
                    outfile.write('\t'+sdi1(BaseMatches_AntiSenseT))
                if GeneBin:
                    for X in BaseMatches_AntiSenseG: outfile.write('\t'+sdi1(X))                    
                if BinBin:
                    for X in BaseMatches_AntiSenseB: outfile.write('\t'+sdi1(X))                    
                if BaseBin:
                    for (N,X) in enumerate(BaseMatches_AntiSenseBa):
                        if T[N]!='*':
                            outfile.write('\t'+sdi1(X))                    
                outfile.write('\r')

            if BaseMatchesTotal:
                outfile.write('BaseMatches_Total'+' '+Experiment)
                if TotalBin:
                    outfile.write('\t'+sdi1(BaseMatches_TotalT))
                if GeneBin:
                    for X in BaseMatches_TotalG: outfile.write('\t'+sdi1(X))                    
                if BinBin:
                    for X in BaseMatches_TotalB: outfile.write('\t'+sdi1(X))                    
                if BaseBin:
                    for (N,X) in enumerate(BaseMatches_TotalBa):
                        if T[N]!='*':
                            outfile.write('\t'+sdi1(X))                    
                outfile.write('\r')

            if ReadMatchesSense:
                outfile.write('ReadMatches_Sense'+' '+Experiment)
                if TotalBin:
                    outfile.write('\t'+sdi1(ReadMatches_SenseT))
                if GeneBin:
                    for X in ReadMatches_SenseG: outfile.write('\t'+sdi1(X))                    
                if BinBin:
                    for X in ReadMatches_SenseB: outfile.write('\t'+sdi1(X))                    
                if BaseBin:
                    for (N,X) in enumerate(ReadMatches_SenseBa):
                        if T[N]!='*':
                            outfile.write('\t'+sdi1(X))                    
                outfile.write('\r')

            if ReadMatchesAntiSense:
                outfile.write('ReadMatches_AntiSense'+' '+Experiment)
                if TotalBin:
                    outfile.write('\t'+sdi1(ReadMatches_AntiSenseT))
                if GeneBin:
                    for X in ReadMatches_AntiSenseG: outfile.write('\t'+sdi1(X))                    
                if BinBin:
                    for X in ReadMatches_AntiSenseB: outfile.write('\t'+sdi1(X))                    
                if BaseBin:
                    for (N,X) in enumerate(ReadMatches_AntiSenseBa):
                        if T[N]!='*':
                            outfile.write('\t'+sdi1(X))                    
                outfile.write('\r')



            if CompositionMatrix:
                CompositionMatrixD=[]
                if TotalBin:
                    CompositionMatrixD=[CompositionMatrixDT]    
                if GeneBin:
                    CompositionMatrixD+=CompositionMatrixDG   
                if BinBin:
                    CompositionMatrixD+=CompositionMatrixDB    
                if BaseBin:
                    for (N,X) in enumerate(CompositionMatrixDBa):
                        if T[N]!='*':
                            CompositionMatrixD.append(X)    
                if CompositionIndexed:
                    for I in range(Composition_Start,Composition_End):
                        for inde in range(4**CompositionTupleLen):
                            Tupl=''
                            ij1=inde
                            for jnde in range(CompositionTupleLen):
                                Tupl=Tupl+B1[ij1 % 4]
                                ij1=ij1 // 4
                            outfile.write(Tupl+'_'+str(I))
                            for Dtup1 in CompositionMatrixD:
                                if (Tupl,I) in Dtup1:
                                    outfile.write('\t'+str(Dtup1[(Tupl,I)]))
                                else:
                                    outfile.write('\t0')
                            outfile.write('\r')
                else:
                    for inde in range(4**CompositionTupleLen):
                        Tupl=''
                        ij1=inde
                        for jnde in range(CompositionTupleLen):
                            Tupl=Tupl+B1[ij1 % 4]
                            ij1=ij1 // 4
                        outfile.write(Tupl+'_'+str(I))
                        for Dtup1 in CompositionMatrixD:
                            if Tupl in Dtup1:
                                outfile.write('\t'+str(Dtup1[Tupl]))
                            else:
                                outfile.write('\t0')
                        outfile.write('\r')

        if PileUp:
            ReadMatches_Sense=[]
            ReadMatches_AntiSense=[]
            Bin_SequenceA=[]
            IniOffsetA=[]
            FinOffsetA=[]
            IniOffsetB=[]
            FinOffsetB=[]
            GeneNamePA=[]
            FBA=[]
            if TotalBin:
                ReadMatches_Sense=[ReadMatches_SenseT]
                ReadMatches_AntiSense=[ReadMatches_AntiSenseT]
                Bin_SequenceA.append(T[1:-1])
                IniOffsetA=[MMAntiSenseOffsetT]
                FinOffsetA=[MMSenseOffsetT]
                IniOffsetB=[MaxMisMatch]
                FinOffsetB=[MaxMisMatch]
                if Genes1>1:
                    GeneNamePA=[ReferenceFileName]
                else:
                    GeneNamePA=[GeneNameA[0]]
                FBA=[1]
            if GeneBin:
                ReadMatches_Sense+=ReadMatches_SenseG
                ReadMatches_AntiSense+=ReadMatches_AntiSenseG
                for (X,Y) in zip(Gene_to_Start,Gene_to_Len):
                    Bin_SequenceA.append(T[X-MaxMisMatch:X+Y+MaxMisMatch])
                    FBA.append(1)
                IniOffsetA+=MMAntiSenseOffsetG
                FinOffsetA+=MMSenseOffsetG
                IniOffsetB+=[MaxMisMatch]*Genes1
                FinOffsetB+=[MaxMisMatch]*Genes1
                GeneNamePA+=GeneNameA
            if BinBin:
                ReadMatches_Sense+=ReadMatches_SenseB
                ReadMatches_AntiSense+=ReadMatches_AntiSenseB
                for (X,Y,Z) in zip(Bin_to_TStart,Bin_to_Len,Bin_to_GStart):
                    Bin_SequenceA.append(T[X:X+Y])
                    FBA.append(Z+1)
                IniOffsetA+=MMAntiSenseOffsetB
                FinOffsetA+=MMSenseOffsetB
                IniOffsetB+=[0]*Bins1
                FinOffsetB+=[0]*Bins1
                for inde11 in range(Bins1):
                    if GeneBin:
                        GeneNamePA.append(GeneNameA[Bin_to_Gene[inde11]])
                    else:
                        GeneNamePA.append(ReferenceFileName)
            for X1 in range(len(ReadMatches_Sense)):
                S11=Bin_SequenceA[X1]
                N11=GeneNamePA[X1][:12]
                AntiS11=AntiSense(S11,filterminus)[::-1]
                LS11=len(S11)
                LN11=len(N11)
                Ruler2=' '*(IniOffsetA[X1]+IniOffsetB[X1])+Ruler(FBA[X1],FBA[X1]+LS11-1-S11.count('-'),N11)
                J=0
                if ReadMatches_Sense:
                    PileUpFile.write('\r#Pile Up-- Sense.  Sequence='+N11+'.'+Experiment+'\r\r')
                    RulerDrawn=False
                    X2=list(ReadMatches_Sense[X1].keys())
                    if X2: X2.sort(key=lambda xx:(xx[1],xx[0]))
                    X4=[]
                    X5=[]
                    for X3 in X2:
                        if X5==[]:
                            LLX=X3[1]+9999
                        else:
                            LLX=min(X5)
                        if PileSep==-1 or X3[1]+IniOffsetB[X1]+IniOffsetA[X1]<=LLX+PileSep:
                            X4.append(' '*(X3[1]+IniOffsetB[X1]+IniOffsetA[X1])+X3[0]+' +'+str(X3[1]+FBA[X1])+' '+str(len(X3[0]))+X3[0][0].upper()+' #'+str(ReadMatches_Sense[X1][X3]))
                            X5.append(len(X4[-1]))
                        else:
                            LLXPos=X5.index(LLX)+0
                            X4[LLXPos]+=' '*(X3[1]+IniOffsetB[X1]+IniOffsetA[X1]-LLX)+X3[0]+' +'+str(X3[1]+FBA[X1])+' '+str(len(X3[0]))+X3[0][0].upper()+' #'+str(ReadMatches_Sense[X1][X3])
                            X5[LLXPos]=len(X4[LLXPos])
                    for Z1 in X4:
                        J+=1
                        if J==40:
                            J=0
                            PileUpFile.write(' '*IniOffsetA[X1]+S11+' '*FinOffsetA[X1]+'\r')
                            PileUpFile.write(Ruler2+'\r')
                            PileUpFile.write(' '*IniOffsetA[X1]+AntiS11+'\r')
                            RulerDrawn=True
                        PileUpFile.write(Z1+'\r')
                    if not(RulerDrawn):
                        PileUpFile.write(' '*IniOffsetA[X1]+S11+' '*FinOffsetA[X1]+'\r')
                        PileUpFile.write(Ruler2+'\r')
                        PileUpFile.write(' '*IniOffsetA[X1]+AntiS11+'\r')
                       
                if ReadMatches_AntiSense:
                    PileUpFile.write('\r#Pile Up-- AntiSense '+Experiment+'  Sequence identity >' + N11+'\r')
                    PileUpFile.write(' '*IniOffsetA[X1]+S11+' '*FinOffsetA[X1]+'\r')
                    PileUpFile.write(Ruler2+'\r')
                    PileUpFile.write(' '*IniOffsetA[X1]+AntiS11+'\r\r')
                    X2=list(ReadMatches_AntiSense[X1].keys())
                    if X2: X2.sort(key=lambda xx:(xx[1],xx[0]))
                    X4=[]
                    X5=[]
                    for X3 in X2:
                        X31=max(X3[1]-len(X3[0])+1+IniOffsetB[X1]+IniOffsetA[X1],0)
                        if X5==[]:
                            LLX=X31+9999
                        else:
                            LLX=min(X5)
                        if PileSep==-1 or X31<=LLX+PileSep:
                            X4.append(' '*(X31)+X3[0][::-1]+' -'+str(X3[1]+FBA[X1])+' '+str(len(X3[0]))+X3[0][0].upper()+' #'+str(ReadMatches_AntiSense[X1][X3]))
                            X5.append(len(X4[-1]))
                        else:
                            LLXPos=X5.index(LLX)
                            X4[LLXPos]+=' '*(X31-LLX)+X3[0][::-1]+' -'+str(X3[1]+FBA[X1])+' '+str(len(X3[0]))+X3[0][0].upper()+' #'+str(ReadMatches_AntiSense[X1][X3])
                            X5[LLXPos]=len(X4[LLXPos])
                    for Z1 in X4:
                        PileUpFile.write(Z1+'\r')
                        J+=1
                        if J==40:
                            J=0
                            PileUpFile.write(' '*IniOffsetA[X1]+S11+' '*FinOffsetA[X1]+'\r')
                            PileUpFile.write(Ruler2+'\r')
                            PileUpFile.write(' '*IniOffsetA[X1]+AntiS11+'\r')
                if J!=0:
                    PileUpFile.write(' '*IniOffsetA[X1]+S11+' '*FinOffsetA[X1]+'\r')
                    PileUpFile.write(Ruler2+'\r')
                    PileUpFile.write(' '*IniOffsetA[X1]+AntiS11+'\r')
                PileUpFile.write('\r'*2)

    if TabularOutput:
        outfile.close()
    if PileUp:
        PileUpFile.close()
    if SimpleOutput:
        SFile.close()

    PickleBytesTotal=30+sum(map(len,GeneNameA))+4*(len(FirstTry)+len(NextTry)+2*Genes1+5*Bins1)+2*LT
    if not(RefAllUpper): PickleBytesTotal+=2*LT
    if GeneBin:
        if Genes1<256:
            PickleBytesTotal+=LT+1
        elif Genes1<65536:
            PickleBytesTotal+=2*(LT+1)
        else:
            PickleBytesTotal+=4*(LT+1)
    if BinBin:
        if Bins1<256:
            PickleBytesTotal+=LT+1
        elif Bins1<65536:
            PickleBytesTotal+=2*(LT+1)
        else:
            PickleBytesTotal+=4*(LT+1)
    
    IndexWritten=False    
    if StoreIndex and not(isfile(IndexFilePath)) and (LT>1000000):
        FreeSpace1=get_free_space(ReferenceFileDir)
        if FreeSpace1>MinFreeGB*(2**30)+PickleBytesTotal:
            Qprint(time0,LogFileName,'Attempting to write Precompiled Index '+IndexFileName)

            try:
                pf1=open(IndexFilePath,mode='wb')
                pp1=cPickle.Pickler(pf1)
                pp1.dump(Bins1)
                pp1.dump(Bases1)
                pp1.dump(Genes1)
                pp1.dump(LT)
                pp1.dump(LY)
                pp1.dump(GeneNameA)
                pp1.dump(RefAllUpper)
                FirstTry.tofile(pf1)
                NextTry.tofile(pf1)
                Gene_to_Start.tofile(pf1)
                Gene_to_Len.tofile(pf1)
                if GeneBin:
                    Address_to_GeneNum.tofile(pf1)
                if BinBin:
                    Address_to_BinNum.tofile(pf1)  
                    Bin_to_GStart.tofile(pf1) 
                    Bin_to_TStart.tofile(pf1) 
                    Bin_to_Len.tofile(pf1)  
                    Bin_to_Gene.tofile(pf1)  
                    Bin_Precedence.tofile(pf1)
                if float(sys.version[:3])<3.0:
                    pf1.write(T)
                    pf1.write(U)
                else:
                    pf1.write(bytes(T,'ascii'))
                    pf1.write(bytes(U,'ascii'))                    
                pf1.close()
                IndexWritten=True
                Qprint(time0,LogFileName,'Successful write of Precompiled Index '+IndexFileName)

            except:
                IndexWritten=False
                Qprint(time0,LogFileName,'Failed in writing Precompiled Index '+IndexFileName)
        else:
            Qprint(time0,LogFileName,'Failed in writing precompiled index due to insufficient disk space, Needed=%i,Available=%i,MinFreeGB=%i'%(PickleBytesTotal,FreeSpace1,MinFreeGB))
        
    if TransposeFile and TabularOutput:
        outfileTN=outfile.name
        outfileUN=outfileTN.strip('.txt')+'TR.txt'
        Qprint (time0,LogFileName,Transpose(outfileTN,outfileUN,100000000))

    ## aardvark: check freespace and see if we can store the index


        
## Debug option 2... to get a dump of local variables in the main loop, uncomment the following five lines

##    DF1=open('dump.txt',mode='w')
##    LocalsD1=locals().keys()
##    for V3 in LocalsD1:
##        DF1.write(str(V3)+'\t'+str(locals()[V3])+'\r')
##    DF1.close()

## Debug option 1: to run as a function, uncomment the next line;        
## MainHQ()

## Debug option 1: to run as a profiled function uncomment the following two lines     

##import cProfile
##cProfile.run('MainHQ()')
##

##11/16/12  Fixed a bug that caused blank cells in an excel text file to overflow number of expected characters
