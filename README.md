# merMapper
---

# merMapper: A tool for mapping the mer frequency of CDR3 regions in antibodies

CDR3 regions play a crucial role in determining the specificity and affinity of antibodies. However, the high mutational capacity and diversity of these regions pose challenges for developing therapeutic agents. Here, we present MerMapper, a bioinformatics tool that maps the mer frequency of CDR3 regions across different contexts, enabling the identification of amino acid residues that may require modification to improve the humanness of therapeutics. MerMapper decomposes CDR3 sequences into mers, assigns a score to each mer based on its frequency in a given context, and facilitates comparison across different species and germline sequences. MerMapper provides maps of human, mouse and camelid antibodies. it also provides utilities to create context specific custom maps from cdr3 regions stored in plain text.  We demonstrate the utility of MerMapper by calculating the mean scores of a set of crd3 populations of different origin. When scored against maps of different genetic and species contexts. We found a ~2 fold decrease in average scores when scoring human cdrh3 against first a human map and secondly against a camelid map. Meanwhile we found a 15 fold increase in mers belonging to the 1 permille most rare mers. This indicates that, while merpopulations of different species are largely overlapping, each context has specific sparsely populated no go zones. 

-----------

# downloade merMapper

```
git clone https://github.com/abmap/merMapper
```

----------

# Use cases

lorem

## Residue prediction

Below is an example of how to load a map and search the merspace to identify feasibel mers in the neighborhood of an input sequence.

```
! pip install numpy
! git clone https://github.com/abmap/merMapper

import itertools, multiprocessing, copy
import numpy as np

mapFile = "mer6Lama.npz"  #7merScoresBR.npz   "mer6Lama.npz"
seq = 'AREGGAAPGARREWYLDL'
MaxMutations = 3  # be carefull the space to explore grows rapidly. in genenral dont exceed 3 mutations 
mask = [0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,0,0,0] 
mask = [1]*len(seq) 
# if provided as an input to genVariantes(), the mask 
# denotes what positions are allowed for mutations. 0 = mutations not allowed 

##############

numAmino = np.array(list("URHKDESTNQCWGPAVILMFY"))
aminoNum = {AA:e for e,AA in enumerate(numAmino)}
def asVekt(seq): return [aminoNum[AA] for AA in seq]
def asSeq(seq): return "".join(numAmino[seq])

def selfStack(vekt,dim):
        return np.vstack(np.array([vekt[t:len(vekt)+1-dim+t] for t in range(dim)]).T) 
    
def getMers(vekt,dim):
        return tuple([tuple(vekt[t:len(vekt)+1-dim+t]) for t in range(dim)])
    
def getScors(seq): 
    global merMap
    if len(seq) < dim: return [],[]
    mers =  selfStack(seq,dim)
    return  merMap[getMers(seq,dim)], mers

def getSeqScore(seq):
    global merMap
    if len(seq) < dim: return None
    return min(merMap[getMers(seq,dim)])

def strRound(val,decs):
    val = str(round(val,decs))
    if len(val) == 1: val+="."
    while len(val)<(decs+2): val+="0"
    return(val)

def mutsForLocs(loc):
    global mutTemplate, mutMask
    locComb = copy.deepcopy(mutTemplate)
    locComb[:,loc] = mutMask
    return locComb

def genVariants(seq,nMutations = 1,mask = None): 
    global mutTemplate, mutMask
    
    if mask == None: mask = [True]*len(seq)
    else: mask = np.array(mask) == True
    if len(mask)!=len(seq): raise Exception("Sequece and mask lenght mismatch!")
    if sum(mask)<nMutations: 
        nMutations = sum(mask)
        print("Mask allows only",sum(mask)," mutations")

    mutIdx = np.arange(len(mask))[mask]
    AAs = (np.arange(1,21).astype(np.uint8))
    
    mutMask = np.array([comb for comb in itertools.product(AAs, repeat = nMutations)])
    mutlocations = np.array([comb for comb in itertools.combinations(mutIdx,nMutations)])
    mutTemplate = np.vstack([np.full(len(mutMask),AA) for AA in seq]).T
    mutTemplate = np.uint8(mutTemplate)
    return np.stack(np.vstack([mutsForLocs(loc) for loc in mutlocations]))

###############

seq = asVekt(seq) # switch to integer representation use asSeq() to switch back
merMap = np.load(mapFile)['arr_0']             
dim = len(np.shape(merMap))

variants = genVariants(seq,MaxMutations, mask)
print("Scoring", len(variants),"variants...") 
minScores = np.array(multiprocessing.Pool().map(getSeqScore,variants,chunksize=1000)) #

# find best mutation 
bestVariant = variants[np.argmax(minScores)]

# print comparison between the best mutation and the original sequence
print("Org seq:     " , asSeq(seq))
print("best variant:", asSeq(bestVariant))

marks = np.array([" "]*len(seq))
marks[np.array(seq) != np.array(bestVariant)] = "^"
print("             ","".join(marks))

print("\nMers before and after mutations to best variant:") 
scoresMut, mersMut = getScors(bestVariant) 
scores, mers = getScors(seq) 
for mer,score,merMut,scoreMut in zip(mers,scores,mersMut,scoresMut): 
    toPrint = str(asSeq(mer)) +" " + strRound(score,5)
    if sum(mer != merMut):
        toPrint += " > " + asSeq(merMut) + " " + strRound(scoreMut,5)
    print(toPrint)

print("\nFound "+str(sum(minScores>0.05))+ 
      " ("+strRound(sum(minScores>0.05)/len(variants)*100,2)+
      "%) vanriants with no mers scoring less than 0.05")

```

This will return the variant within `nMutations` residue substitutions that will raise the lowest mer frequence score the moast.

```console
Scoring 6528000 variants...
Org seq:      AREGGAAPGARREWYLDL
best variant: ARGGGAAPGAGRSWYLDL
                ^       ^ ^     

Mers before and after mutations to best variant:
AREGGA 0.03219 > ARGGGA 0.25020
REGGAA 0.00000 > RGGGAA 0.14488
EGGAAP 0.00000 > GGGAAP 0.21155
GGAAPG 0.10096
GAAPGA 0.06751
AAPGAR 0.17022 > AAPGAG 0.20292
APGARR 0.05235 > APGAGR 0.10990
PGARRE 0.00000 > PGAGRS 0.13248
GARREW 0.00000 > GAGRSW 0.13248
ARREWY 0.10096 > AGRSWY 0.26747
RREWYL 0.00000 > GRSWYL 0.24504
REWYLD 0.00000 > RSWYLD 0.06751
EWYLDL 0.00000 > SWYLDL 0.08009

Found 18 (0.00%) vanriants with no mers scoring less than 0.05
```
-----

## Build your own maps

### Citation

```
tba
```  
