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

import numpy as np
import copy 

#select sequence map and search strategi
# human cdr3s: 
# seq = 'AKDRNHYGSGSYFDY' 
# seq = 'LYRAEPPRGGYYFDY'
# seq = 'ARDRALRGRSWFDP'
# seq = 'AKDGDYYDSSGALDY'
# seq = 'AREGGAAPGARREWYLDL'
# seq = 'ARDLVIRGIGASDY'
# seq = 'ASGYFDHFFDF' 
seq = 'TISIGSGGVIGGFDY'

 # lama cdr3
# seq = 'AAGPPRLCTLSVWTVYDY'  

mapFile = "mer6Lama.npz"
MaxMutations = 2  # be carefull the space to explore grows ~ x^n. in genenral dont exceed 3 for local search, and 2 for global 
localSearch = False # searce only the space around the mer with the lowest score.

# let the code run ...

# for t in [1,39,346,432,35,653,10000, 1003]:
#     print(asSeq(seqs[t]))

numAmino = np.array(list("URHKDESTNQCWGPAVILMFY"))
aminoNum = {AA:e for e,AA in enumerate(numAmino)}
def asVekt(seq): return [aminoNum[AA] for AA in seq]
def asSeq(seq): return "".join(numAmino[seq])


def selfStack(vekt,dim):
        return np.vstack(np.array([vekt[t:len(vekt)+1-dim+t] for t in range(dim)]).T) 
    
def getScors(seq,merMap): 
    dim = len(np.shape(merMap))
    if len(seq) < dim: return [],[]
    mers =  selfStack(seq,dim)
    return [merMap[tuple(mer)] for mer in mers], mers

def genMutations(seq,startPos,endPos,nMutations = 1): 

    mutations =  [copy.deepcopy(seq) for mut in range((endPos-startPos)*20)]
    mutations = []
    for pos in range(startPos,endPos):
      for AA in np.arange(20)+1: 
          mutant = copy.deepcopy(seq) 
          mutant[pos] = AA
          mutations.append(mutant)
          
    if nMutations == 1: return mutations
    else: 
        superMuts = []
        for mut in mutations: 
           superMuts += genMutations(mut,startPos,endPos,nMutations = nMutations-1)
        return superMuts   
    
def strRound(val,decs):
    val = str(round(val,decs))
    if len(val) == 1: val+="."
    while len(val)<(decs+2): val+="0"
    return(val)

seq = asVekt(seq)
merMap = np.load(mapFile)['arr_0']             

# calc scores for mers                     
scores, mers = getScors(seq,merMap) 

# define region to explore
if localSearch:
    startPos = np.argmin(scores)
    endPos = startPos+len(np.shape(merMap))
else: startPos = 0; endPos = len(seq)
   
mutations = genMutations(seq,startPos,endPos,MaxMutations) 

freqs = [np.min(getScors(mut,merMap)[0]) for mut in mutations]
bestMutation = mutations[np.argmax(freqs)]


print("Org seq:  " , asSeq(seq))
print("Mutations:", asSeq(bestMutation))


marks = np.array([" "]*len(seq))
marks[np.array(seq) != np.array(bestMutation)] = "↑"
print("          ","".join(marks))


print("")
print("Mers before and after mutations:") 
scoresMut, mersMut = getScors(bestMutation,merMap) 
for mer,score,merMut,scoreMut in zip(mers,scores,mersMut,scoresMut): 
    toPrint = str(asSeq(mer)) +" " + strRound(score,5)
    if sum(mer != merMut):
        toPrint += " → " + asSeq(merMut) + " " + strRound(scoreMut,5)
    print(toPrint)

        

        
```

This will return the `nMutations` residue substitutions that will raise the minimum mer frequence the moast.

```console
Org seq:   TISIGSGGVIGGFDY
Mutations: TISRGSGGVGGGFDY
              ↑     ↑     

Mers before and after mutations:
TISIGS 0.00000 → TISRGS 0.21676
ISIGSG 0.03219 → ISRGSG 0.19986
SIGSGG 0.05235 → SRGSGG 0.20589
IGSGGV 0.22170 → RGSGGV 0.24678
GSGGVI 0.03219 → GSGGVG 0.27576
SGGVIG 0.00000 → SGGVGG 0.79968
GGVIGG 0.26747 → GGVGGG 0.17857
GVIGGF 0.00000 → GVGGGF 0.12551
VIGGFD 0.00000 → VGGGFD 0.17450
IGGFDY 0.14488 → GGGFDY 0.19329
```
-----

## Build your own maps

### Citation

```
tba
```  
