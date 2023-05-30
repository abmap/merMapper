import numpy as np
import copy 

#select sequence map and search strategi
# human cdr3s: 
# seq = 'AKDRNHYGSGSYFDY' 
# seq = 'LYRAEPPRGGYYFDY'
# seq = 'ARDRALRGRSWFDP'
# seq = 'AKDGDYYDSSGALDY'
# seq = 'TISIGSGGVIGGFDY'
# seq = 'AREGGAAPGARREWYLDL'
# seq = 'ARDLVIRGIGASDY'
# seq = 'ASGYFDHFFDF'  

 # lama cdr3
# seq = 'AAGPPRLCTLSVWTVYDY'  

mapFile = "mer6Lama.npz"
MaxMutations = 2  # be carefull the space to explore grows ~ x^n. in genenral dont exceed 3 for local search, and 2 for global 
localSearch = False # searce only the space around ther with the lowest score

# let the code run ...git

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

        