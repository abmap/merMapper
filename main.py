import itertools, multiprocessing, copy
import numpy as np

mapFile = "mer6Lama.npz"  #7merScoresBR.npz   "mer6Lama.npz"
seq = 'AREGGAAPGARREWYLDL'
MaxMutations = 3  # be carefull the space to explore grows rapidly. in genenral dont exceed 3 mutations 
mask = [0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,0,0,0] 
mask = [1]*len(seq) 
# if provided as an input to genVariantes(), 
# the mask denotes what positions are allowed for mutations. 0 = mutations not allowed 

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
