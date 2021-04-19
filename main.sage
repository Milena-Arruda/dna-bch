import numpy as np

def read_dna_database(path,N = None):
    '''
    Loading the database from a given path. 
    N sequences will be choosen randomly from database.
    '''
    database = []
    with open(path, 'r') as f:
        seq = ''
        for line in f:
            if line!= '\n' and line[0]!= '>':
                seq += line[:-1]
            elif seq != '':
                if set(seq.upper()) == {'A','C','G','T'}:
                    database.append(seq.upper())
                seq = ''

    if N is None:
        N = 500
    else:
        pass

    ind = np.random.choice(len(database),N,replace=False) 
    dna_nuc = []
    for i in ind:
        dna_nuc.append(database[i]) 
    
    return dna_nuc   

def dna_to_gf(dna,label):
    '''
    This function returns the DNA elements mapped according to 
    the alphabet defined in label
    '''
    dnaDict = {'A': label[0], 'C': label[1],
               'G': label[2], 'T': label[3]}

    return [dnaDict[i] for i in dna]

def coprime(k,n):
    ''' 
    This function returns the elements of a given vector k 
    that are coprime of n 
    '''
    return[i for i in k if gcd(i, n) == 1]

def cyclotomic(k,m,q):
    ''' 
    This function returns the minimun expoent of cyclotomic 
    coset of alpha^k 
    '''
    return min([mod(k*q^e,n) for e in range(m)])

def class_cyclotomics(n,m,q): 
    ''' 
    This function returns the vector with the minimun expoent 
    of each one of cyclotomic coset of alpha^k for k in range 
    from 1 to n-1 
    '''
    return sorted(set([cyclotomic(k,m,q) for k in range(1,n)]))

def range_ell(n,m,q): 
    ''' 
    This function returns the range of ell satisfying the 
    conditions discussed in section 3.B 
    '''
    A = class_cyclotomics(n,m,q)
    primes = coprime(A,n)
    return sorted(set([min([cyclotomic(i,m,q),cyclotomic(-i,m,q)]) for i in primes])) 

def solve(a,b,n):
    '''
    This function returns the solution of the following 
    equation: a*x = b mod n 
    '''
    solution = []
    t = gcd(a,n)
    for i in range(n):
        if (a*i) % n == b % n:
            if i not in solution:
                solution.append(i) 
            else:
                break
    return solution

def dictionary_of_codes(n,d):
    '''
    This function returns a dictionary with all codes with 
    length n and design distance d. 
    The key words of dictionary are the generator polynomial. 
    The values of dictionary are a list with two elements:
    1) C (the BCH code according to SageMath structure)
    2) 0 (this numeric element will be used to count the 
    number of DNA sequences that code C identifies)
    '''
    dictCodes = {}
  
    for ell in range_ell(n,m,q):
        for beta in range(n):
            C = codes.BCHCode(Fq,n,d,offset=beta,jump_size=ell)
            g = C.generator_polynomial()
            if g not in dictCodes:
                dictCodes[g] = [C,0]

    return dictCodes

def find_bch_codes(rx):
    '''
    This function returns all the generator polynomials of 
    BCH codes that identifies a DNA sequence. 
    Furthermore, this function updates the numeric element of 
    dictionary of BCH codes.
    '''
    gBCH = []
    position = []

    for code in dictCodes:
        C = dictCodes[code][0]
        beta = int(C.offset())
        ell = int(C.jump_size())
    
        if rx(a^(order*beta)) == 0 and rx(a^(order*(beta+ell))) == 0:
            # codeword of C
            dictCodes[code][1] += 1
            gBCH.append(i)
            position.append('x')

        elif rx(a^(order*beta)) != 0 and rx(a^(order*(beta+ell))) != 0:
            Lambda = rx(a^(order*(beta+ell)))/rx(a^(order*beta))
            multiplicativegroup = [a^(order*i) for i in range(n)]
            if Lambda in multiplicativegroup:
                pos = solve(ell,multiplicativegroup.index(Lambda),n)
            if pos:
                for p in pos:
                    rxcurr = rx[p]
                    rxnew = rxcurr - rx(a^(order*beta))/a^(order*beta*p)
                    rrx = list(r)
                    rrx[p] = rxnew
                    if rxnew in label:
                        if FPoly(rrx)(a^(order*(beta+ell))) == 0 and FPoly(rrx)(a^(order*beta)) == 0:   
                            dictCodes[code][1] += 1
                            gBCH.append(code)
                            position.append((rx - FPoly(rrx)).degree())
                
    return gBCH, position

def m_find(n,q):
    '''
    This function returns the minimum m such that 
    q^m - 1 mod n = 0
    '''
    m = 2
    while (q^m-1) % n != 0:
        m += 1
    
    return m

#*************************************************************
#                          main                      
#*************************************************************
#path = '' #Indicate the pathway of database
#N = 500
#database = read_dna_database(path, N)

N = 1
database = ['ctgatccttcaagcg'] #Example of Section IV.B
## STEP 0: BCH CODE PARAMETERS AND ALGEBRAIC STRUCTURE
n, d = len(database[0]), 3
t = int(np.floor((d-1)/2))
q = 4
m = m_find(n,q)
order = int((q^m-1)/n)

Fq.<b> = GF(q,prefix='b')
label = [0,1,b,b+1]
Fqm.<a> = GF(q^m,prefix='a')
FPoly.<x> = Fqm[]

dictBCH = {}
print('Please wait while the program is loading...')
dictCodes = dictionary_of_codes(n,d)
position = []

numberCodes = []
flag = 0
P = 0.749
for i in range(N):
    seq = database[i].upper()

    print('analyzing sequence ' + str(i))
    ## STEP 1: DNA SEQUENCE
    r = dna_to_gf(seq,label) 
    rx = FPoly(r)

    ## STEP 2: FIND BCH CODES
    gBCH, pos = find_bch_codes(rx) #gBCH: Table 1 column 3 and pos: Table 1 column 4
    numberCodes.append(len(gBCH))
    if pos != []:
        position.append(pos)


print('Searching common codes between DNA sequences...')
## STEP 3: COMMON CODES
dictCodesValues = [i[1] for i in dictCodes.values()]
M = max(dictCodesValues)
g = list(dictCodes.keys())[list(dictCodesValues).index(M)] 

## STEP 4: STATISTICAL ANALYSIS
totalSeq = len(database)
porcent = float(M/(totalSeq)*100)
print("%.2f" % porcent + "% of sequences identified by")
print("[" + str(n) + "," + str(n-g.degree()) + ",3] C_BCH with")
print("g(x) = "+ str(g))

## STEP 5: BOXPLOT INFORMATION
#np.quantile(numberCodes, 0.5)              
#np.quantile(numberCodes, 0.25)             
#np.quantile(numberCodes, 0.75)             
#max(numberCodes)                           
#min(numberCodes) 