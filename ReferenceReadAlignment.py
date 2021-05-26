
# coding: utf-8

# ## Burrows Wheeler transformation
# 

# In[1]:


def rotations(t):
    ''' Return list of rotations of input string t '''
    words = list(t)
    retList = []

    for i in range(len(words)):
        if(i % 10000 == 0):
            print("rot step", i)
            
        word = t[-1] + t[:-1]
        new = ''.join(word)
        t = new
        retList.append(new)
        i += 1
    return retList

def bwm(t):
    ''' Return lexicographically sorted list of t’s rotations '''
    x = sorted(rotations(t))
    return x

def bwtViaBwm(t):
    ''' Given T, returns BWT(T) by way of the BWM '''
    x = ''.join(map(lambda x: x[-1], bwm(t)))
    return x


# In[2]:


t = 'banana$'
b = bwtViaBwm(t)
b


# In[3]:


# confirm that the representation of the first column above is sensible
print('\n'.join(bwm(t)))


# ## FM Index

# ### First column

# In[4]:


def firstCol(tots):
    ''' Return map from character to the range of rows prefixed by the character. '''
    first = {}
    totc = 0
    for c, count in sorted(tots.items()):
        first[c] = (totc, totc + count)
        totc += count
    return first


# ### Occurrences

# In[5]:


def occBwt(bw):
    tots = dict()
    occ = dict()
    for c in bw:
        if c not in tots:
            tots[c] = 0
            occ[c] = []
            occ[c].append(0)
        tots[c] += 1
    
    i = 0
    for c in bw:
        # u occ[c] onda dodaj +1 u odnosu na prosli
        occ[c].append(occ[c][i] + 1)
        
        # u svaki ostali occ prepisi predhodnu vrijednost 
        for ch in occ:
            if(ch != c):
                occ[ch].append(occ[ch][i])
        
        i+=1
        
    return occ, tots


# In[6]:


occ, tots = occBwt(b)
occ


# In[7]:


c = firstCol(tots)
c


# ### Suffix array

# In[8]:


def suffixArray(s):
    """ Given T return suffix array SA(T).  We use Python's sorted
        function here for simplicity, but we can do better. """
    satups = sorted([(s[i:], i) for i in range(len(s))])
    # Extract and return just the offsets
    return map(lambda x: x[1], satups)


# In[9]:


sa = list(suffixArray(t))
sa


# #### Searching for rows in bw matrix which start with given substring/read [x,y)

# In[10]:


def rangeInMatrix(c, occ, substring):    
    start = c.get(substring[-1],(-100,-100))[0]
    end = c.get(substring[-1],(-100,-100))[1]
    
    if(start < 0 or end < 0 or start >= end ):
        return -1,-1
    
    for i in range(len(substring)-2,-1,-1):        
        if(start < 0 or end < 0 or start >= end ):
            return -1,-1
        
        start = c.get(substring[i],(-100,-100))[0] + occ.get(substring[i])[start]
        end = c.get(substring[i],(-100,-100))[0] + occ.get(substring[i])[end]

    if(start < 0 or end < 0 or start >= end ):
        return -1,-1
        
    return start, end


# In[11]:


start, end = rangeInMatrix(c, occ, 'ana') 
print(start, end)


# #### Searching for position of substring/read in initial string/reference genome

# In[12]:


posOfSubstring = []
while start < end:
    posOfSubstring.append(sa[start])
    start +=1

print(sorted(posOfSubstring))


# <br/>

# ## Global alignment

# In[13]:


def scoringMatrix(a, b, match, mismatch, gap):
    if a == b: return match
    if a == '_' or b == '_' : return gap
    return mismatch


# In[14]:


import numpy

def globalAlignment(x, y, s, match, mismatch, gap):
    D = numpy.zeros((len(x) + 1, len(y) + 1), dtype=int)
    
    for i in range(1, len(x) + 1):
        D[i,0] = D[i-1,0] + s(x[i-1], '_', match, mismatch, gap)  
    for j in range(1, len(y)+1):
        D[0,j] = D[0,j-1] + s('_', y[j-1], match, mismatch, gap)
    
    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            D[i,j] = max(D[i-1,j]   + s(x[i-1], '_', match, mismatch, gap),
                         D[i,j-1]   + s('_', y[j-1], match, mismatch, gap), 
                         D[i-1,j-1] + s(x[i-1], y[j-1], match, mismatch, gap))
            
    # function returns table and global alignment score
    #alignment score is in cell (n,m) of the matrix
    return D, D[len(x),len(y)]


# In[15]:


def traceback(x,y,V,s, match, mismatch, gap):
    
    # initializing starting position cell(n,m)
    i=len(x)
    j=len(y)
    
    # initializing strings we use to represent alignments in x, y, edit transcript and global alignment
    ax, ay, am, tr = '', '', '', ''
    
    # exit condition is when we reach cell (0,0)
    while i > 0 or j > 0:
        # calculating diagonal, horizontal and vertical scores for current cell
        d, v, h = -100, -100, -100
        
        if i > 0 and j > 0:
            delta = 1 if x[i-1] == y[j-1] else 0
            d = V[i-1,j-1] + s(x[i-1], y[j-1], match, mismatch, gap)  # diagonal movement   
        if i > 0: v = V[i-1,j] + s(x[i-1], '_', match, mismatch, gap)  # vertical movement
        if j > 0: h = V[i,j-1] + s('_', y[j-1], match, mismatch, gap)  # horizontal movement
          
        
        # backtracing to next (previous) cell
        if d >= v and d >= h:
            ax += x[i-1]
            ay += y[j-1]
            if delta == 1:
                tr += 'M'
                am += '|'
            else:
                tr += 'R'
                am += ' '
            i -= 1
            j -= 1
        elif v >= h:
            ax += x[i-1]
            ay += '_'
            tr += 'D'
            am += ' '
            i -= 1
        else:
            ay += y[j-1]
            ax += '_'
            tr += 'I'
            am += ' '
            j -= 1
            
    alignment='\n'.join([ax[::-1], am[::-1], ay[::-1]])
    return alignment, tr[::-1]


# In[16]:


x='TACGTCAGC'
y='TATGTCATGC'


# In[17]:


D, alignmentScore = globalAlignment(x, y, scoringMatrix,1,-2,-7)


# In[18]:


print(alignmentScore)


# In[19]:


print(D)


# In[20]:


alignment, transcript = traceback(x, y, D, scoringMatrix,1,-2,-7)


# In[21]:


print(alignment)
print(transcript)


# ## Seed & Extend

# In[22]:


from Bio import SeqIO
from Bio.Seq import Seq
import os


# In[23]:


# import argparse

# parser = argparse.ArgumentParser()
# parser.add_argument("--fasta", help="Path to the FASTA file")
# parser.add_argument("--fastq", help="Path to the FASTQ file")
# parser.add_argument("--margin", help="Length of the margin for ")
# parser.add_argument("--match", help="Array")

# parser.add_argument("--seed-length", help="Length of the seed (defaults to 10)")



# args = parser.parse_args()
# print(args.fasta)
# print("Hello world")


# In[24]:


fastaFilename = '/sbgenomics/project-files/reference.fasta'
fastqFilename = '/sbgenomics/project-files/reads.fastq'

lenOfSeed = 10
margin = 2

matchs = [0, 1, 2]
mismatchs = [-3, -2]
gaps = [-7, -5]


# In[25]:


referenceGenome = ''
fasta_sequences = SeqIO.parse(open(fastaFilename),'fasta')
for fasta in fasta_sequences:
    name, referenceGenome = fasta.id, str(fasta.seq)

names = []
reads = []
fastq_sequences = SeqIO.parse(open(fastqFilename),'fastq')#fastq_sequences = SeqIO.parse(open(fastqFilename),'fastq')
for fastq in fastq_sequences:
    name, sequence = fastq.id, str(fastq.seq)
    names.append(name)
    reads.append(sequence)


# In[26]:


def reverseComplement(read):
    seq = Seq(read)
    return seq.reverse_complement()


# In[27]:


def seedAndExtend(occ, c, sa, r, posOfSubstring, transcripts, alignmentScores, match, mismatch, gap):
    start, end = rangeInMatrix(c, occ, r[0:lenOfSeed]) 
    
    while start < end:
        posOfSubstring.append(sa[start])
        start +=1
    
    for startPos in posOfSubstring:
        startPos = startPos +lenOfSeed
        endPos = startPos + len(r) - lenOfSeed + margin
        
        if(endPos > len(referenceGenome)):
            endPos = len(referenceGenome)
        
        if(startPos > endPos):
            continue
        
        D, alignmentScore = globalAlignment(r[lenOfSeed:], referenceGenome[startPos:endPos], scoringMatrix, match, mismatch, gap)
        alignment, transcript = traceback(r[lenOfSeed:], referenceGenome[startPos:endPos], D, scoringMatrix, match, mismatch, gap)
        transcripts.append(transcript)
        alignmentScores.append(alignmentScore)
    
    return posOfSubstring, transcripts, alignmentScores


# In[28]:


import pickle

def getFileName(match, mismatch, gap):
    return "output_" + str(match) + "_" + str(mismatch) + "_" + str(gap) + ".txt"

def preCompute(referenceGenome):    
    bwTransformation = bwtViaBwm(referenceGenome+'$')
    occ, tots = occBwt(bwTransformation)
    c = firstCol(tots)
    sa = list(suffixArray(referenceGenome))
    
    pickle.dump(c, open( "c.dumpfile", "wb" ))
    pickle.dump(occ, open( "occ.dumpfile", "wb" ))
    pickle.dump(sa, open( "sa.dumpfile", "wb" ))
    
    return c, occ, sa


def getTriplets(c, occ, sa, match, mismatch, gap):
    i = 0
    fileName = getFileName(match, mismatch, gap)
    resultTuples = []
    
    with open(fileName, "w") as file1:
        for read in reads:          
            posOfSubstring = []
            transcripts = []
            alignmentScores = []
            posOfSubstring, transcripts, alignmentScores = seedAndExtend(occ,c,sa, read, posOfSubstring, transcripts, alignmentScores, match, mismatch, gap)
            posOfSubstring, transcripts, alignmentScores = seedAndExtend(occ,c,sa, reverseComplement(read), posOfSubstring, transcripts, alignmentScores, match, mismatch, gap)
            
            if len(posOfSubstring) > 0:
                alignmentScores, transcripts, posOfSubstring = (list(t) for t in zip(*sorted(zip(alignmentScores, transcripts, posOfSubstring), reverse = True)))
            else:
                #DEBUG
                print("NOT FOUND -->",names[i], read)

            resultTuple = [read, posOfSubstring, transcripts, alignmentScores]
            file1.write(str(resultTuple)+ '\n')
            
            resultTuples.append(resultTuple)
            i +=1
            
    return resultTuples


# #### read pickle files

# In[29]:


def readPickleFiles():
    c = pickle.load( open( "c.dumpfile", "rb" ))
    occ = pickle.load(open( "occ.dumpfile", "rb" ))
    sa = pickle.load(open( "sa.dumpfile", "rb" )) 
    return c, occ, sa


# In[30]:


# c, occ, sa = preCompute(referenceGenome)
c, occ, sa = readPickleFiles()


# In[31]:


print(c)


# In[32]:


def getMatrixRow(numrot,refGenome):
    row = refGenome
    i = 0 
    while i< numrot:
        row = row[1:] + row[0]
        i+=1
    return row


# In[33]:


#DEBUG
print(names[3608]) #error on this read
                    
#read:         TGCGGTGGCTTACGCCTGTAATCCCAGCCGAGGCAGGCAGAATGAGCTCATGAGTTCGAGACCAGCATTGCAAAACCCTGTCTCTACTAAAAAAACAAAA
#reverse read: TTTTGTTTTT TTAGTAGAGACAGGGTTTTGCAATGCTGGTCTCGAACTCATGAGCTCATTCTGCCTGCCTCGGCTGGGATTACAGGCGTAAGCCACCGCA
#                        _TTAGTAGAGACAGGGTTTTGCAATGCTGGTCTCGAACTCATGAGCTCAT_T_CTGCCTGCCTCGGCTGGGATTACAGGCGTAAGCCACCGCA


# In[34]:


for match in matchs:
    for mismatch in mismatchs:
        for gap in gaps:
            print("Continuing")
            #result = getTriplets(c, occ, sa, match, mismatch, gap)
            #outputs[getFileName(match, mismatch, gap)] = result


# ## Result visualization

# In[36]:


# Reading data from our .txt file
import re
import ast
from datetime import datetime


def readDataFromFile(filename):
    data = []
    
    i = 0
    with open(filename,"r") as file:
        for row in file.readlines():
            scores = ''
            
            row = row.split('[')

            scores = '['+ row[4][0:-2]
            scores = ast.literal_eval(scores)

            data.append(scores)
            
        
    return data


# In[37]:


#Reading scores from SAM file

bwaMemFileName = '/sbgenomics/project-files/bwa_mem_reads.sam'
bwaNameToScore = {}
        
import pysam
samfile = pysam.AlignmentFile(bwaMemFileName, "rb")

for read in samfile:
    match = re.search("'AS', (\d+)", str(read))
    if(match):
        bwaNameToScore[read.qname] = int(match.group(1))


# In[38]:


import matplotlib.pyplot as plt
import numpy as np


# In[39]:


# plot graph
outputs = {}

for match in matchs:
    for mismatch in mismatchs:
        for gap in gaps:

            fileName = getFileName(match, mismatch, gap)
        
            data = readDataFromFile(fileName)
            outputs[fileName] = data
                
            y = []
            yBwa = []
            x = 0
            curr = 0
            
            for readScores in data:
                #avoid unmatched reads
                if(len(readScores) > 0):
                    bestScore = max(readScores)
                    
                    y.append(bestScore)
                    yBwa.append(bwaNameToScore[names[curr][0:-2]])
                    x = x + 1
                curr += 1
            
            xNp = np.arange(start=0, stop=x, step=1)
            titleString = "Custom, " + str(match) + ", " + str(mismatch) + ", " + str(gap)
            plt.figure(figsize=(10,10))
            plt.plot(xNp[0:20], y[0:20], label = "Custom")
            plt.plot(xNp[0:20], yBwa[0:20], label = "BWA MEM")
            plt.xlabel('Reads')
            plt.ylabel('Alignment score')
            plt.title(titleString)
            plt.legend()
            plt.show()
            
            plt.savefig("Graph_" + fileName + ".png")


# In[40]:


# write table to csv

import csv

with open('results.csv', 'w', newline='') as csvfile:
    csvWriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
    row = []
    row.append("SeqId")
    row.append("Sequence")
    row.append("BWA MEM")
    for match in matchs:
        for mismatch in mismatchs:
            for gap in gaps:
                row.append("Custom " + str(match) + " " + str(mismatch) + " " + str(gap))
    
    #headers
    csvWriter.writerow(row)
    
    for i in range(len(names)):
        row = []
        
        row.append(names[i][0:-2])
        row.append(reads[i])
        
        #append BWA MEM output
        row.append(bwaNameToScore[names[i][0:-2]])
        
        #append Custom outputs (all 12)
        for match in matchs:
            for mismatch in mismatchs:
                for gap in gaps:
                    fileName = getFileName(match, mismatch, gap)
                
                    if(len(outputs[fileName][i]) > 0):
                        row.append(outputs[fileName][i][0])
                    else:
                        row.append(None)
        
        csvWriter.writerow(row)

