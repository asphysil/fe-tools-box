from itertools import count
import numpy as np 
from random import seed
from random import sample
from random import shuffle, randrange

def GtypeAFE_order(nn, nb, latt, atms_order):
    dist = np.zeros(nb)    
    atm = np.zeros((nb,3))
    dist_pbc = np.zeros(27)
    v2 = np.zeros(3)

    for i in range(0, nb):
        #atms_order[i][0] = 1.0
        atm[i][0] = atms_order[i][1]
        atm[i][1] = atms_order[i][2]
        atm[i][2] = atms_order[i][3]

    dmin = min_dist(nb, atm) + 1.0

    for i in range(0, nb):
        if int(atms_order[i][0]) == 1:
            v1 = atms_order[i,1:4]
            
            for j in range(0, nb):
                m = 0
                for lx in range(-1, 2):
                    for ly in range(-1,2):
                        for lz in range(-1,2):
                            v2[0] = atms_order[j][1] + lx*latt[0][0]
                            v2[1] = atms_order[j][2] + ly*latt[1][1]
                            v2[2] = atms_order[j][3] + lz*latt[2][2]
                            dist_pbc[m] = np.linalg.norm(v1-v2)
                            m = m + 1
                dist[j] = np.min(dist_pbc) 
            for k in range(0, nn+1):
                m = np.argmin(dist)
                d = dist[m]
                #print(d, dmin)
                if d <= dmin:
                    dist[m] = 10000.0
                    if d > 1.0:
                        if int(atms_order[m][0]) == 1:
                           #print('ok')
                           atms_order[m][0] = 2.0
        else:
            pass
    #print(atms_order)
    return atms_order

#-----------------------------G-type-AFM-order-version2--------------------------
def GtypeAFE_order_v2(nn, nb, latt, atms):

    dmin = min_dist(nb, atms) + 1.0 # minimum distance between two atoms
    dist = np.zeros(nb)    
    atms_order = np.zeros((nb,4)) #
    dist_pbc = np.zeros(27)
    v2 = np.zeros(3)
    vec_n = np.zeros(3)

    for i in range(0, nb):
        atms_order[i][0] = 1.0
        atms_order[i][1] = atms[i][0]
        atms_order[i][2] = atms[i][1]
        atms_order[i][3] = atms[i][2]

    

    for i in range(0, nb):
        if int(atms_order[i][0]) == 1:
            v1 = atms[i]
            
            for j in range(0, nb):
                m = 0
                for lx in range(-1, 2):
                    for ly in range(-1,2):
                        for lz in range(-1,2):
                            vec_n[0] = lx 
                            vec_n[1] = ly 
                            vec_n[2] = lz 
                            T = np.matmul(np.transpose(latt), vec_n)
                            #print(T) 
                            v2[0] = atms[j][0] + T[0]
                            v2[1] = atms[j][1] + T[1]
                            v2[2] = atms[j][2] + T[2]

                            dist_pbc[m] = np.linalg.norm(v1-v2)
                            m = m + 1
                dist[j] = np.min(dist_pbc) 
            for k in range(0, nn+1):
                m = np.argmin(dist)
                d = dist[m]
                #print(d, dmin)
                if d <= dmin:
                    dist[m] = 10000.0
                    if d > 1.0:
                        if int(atms_order[m][0]) == 1:
                           #print('ok')
                           atms_order[m][0] = 2.0
        else:
            pass
    #print(atms_order)
    return atms_order

def GtypeAFE_order_v3(nn, nb, latt, atms_order):
    dist = np.zeros(nb)    
    atm = np.zeros((nb,3))
    dist_pbc = np.zeros(27)
    v2 = np.zeros(3)

    for i in range(0, nb):
        if int(atms_order[i][0])==2:
           atms_order[i][0] = 1
        
        atm[i][0] = atms_order[i][1]
        atm[i][1] = atms_order[i][2]
        atm[i][2] = atms_order[i][3]

    dmin = min_dist(nb, atm) + 1.0

    for i in range(0, nb):
        if int(atms_order[i][0]) == 1:
            v1 = atms_order[i,1:4]
            
            for j in range(0, nb):
                m = 0
                for lx in range(-1, 2):
                    for ly in range(-1,2):
                        for lz in range(-1,2):
                            v2[0] = atms_order[j][1] + lx*latt[0][0]
                            v2[1] = atms_order[j][2] + ly*latt[1][1]
                            v2[2] = atms_order[j][3] + lz*latt[2][2]
                            dist_pbc[m] = np.linalg.norm(v1-v2)
                            m = m + 1
                dist[j] = np.min(dist_pbc) 
            for k in range(0, nn+1):
                m = np.argmin(dist)
                d = dist[m]
                #print(d, dmin)
                if d <= dmin:
                    dist[m] = 10000.0
                    if d > 1.0:
                        if int(atms_order[m][0]) == 1:
                           #print('ok')
                           atms_order[m][0] = 2.0
        else:
            pass
    #print(atms_order)
    return atms_order

def rock_salt_order_random(natms, ndopant, atms):
    m = int(natms/2)
    na = np.zeros(m, dtype=np.int)
    nb = np.zeros(m, dtype=np.int)
    m1 =0
    m2 =0 
    for i in range(0, natms):
        m = int(atms[i][0])
        if m ==1:
           na[m1] = i 
           m1 = m1 + 1
           atms[i][0] = 1
        else:
            nb[m2] = i
            m2 = m2 + 1
            atms[i][0] = 1

    #print('--', ndopant, m2)

    sequence_a = [i for i in range(m1)]
    np.random.seed(randrange(0, 101, 2))
    #seed(seed1)
    shuffle(sequence_a)

    np.random.seed(randrange(0, 1001, 2))
    #seed(seed2)
    sequence_b = [i for i in range(m2)]
    shuffle(sequence_b)

    for i in range(0, ndopant):
        m = nb[sequence_b[i]]
        atms[m][0] = 3


    return atms 

def rock_salt_order_fix(shift, latt, natms, ndopant, atms):
    m = int(natms/2)
    na = np.zeros(m, dtype=np.int)
    nb = np.zeros(m, dtype=np.int)
    m1 =0
    m2 =0 
    for i in range(0, natms):
        m = int(atms[i][0])
        if m ==1:
           na[m1] = i 
           m1 = m1 + 1
           atms[i][0] = 1
        else:
            nb[m2] = i
            m2 = m2 + 1
            atms[i][0] = 1

    #print('--', ndopant, m2)

    #sequence_a = [i for i in range(m1)]
    #np.random.seed(randrange(0, 101, 2))
    #seed(seed1)
    #shuffle(sequence_a)

    #np.random.seed(randrange(0, 1001, 2))
    #seed(seed2)
    sequence_b = [i for i in range(m2)]
    #shuffle(sequence_b)

    for i in range(0, ndopant):
        m = nb[sequence_b[i+shift]]
        atms[m][0] = 3


    return atms

def ti_ba_order(natms, nti, latt, atms, atmti):

    nhalf = int(natms/2)
    na1 = np.zeros(nhalf, dtype=np.int)
    na2 = np.zeros(nhalf, dtype=np.int)

    dist_pbc = np.zeros(27)

    dist = np.ones(nhalf)*10000

    nba = np.zeros(nti, dtype=np.int)

    w = np.zeros(3)
    m1 =0
    m2 =0 
    for i in range(0, natms):
        m = int(atms[i][0])
        if m ==1:
           na1[m1] = i 
           m1 = m1 + 1
           atms[i][0] = 1
        else:
            na2[m2] = i
            m2 = m2 + 1
            atms[i][0] = 1
    print(m1, m2) 

    for i in range(0, nti):
        v1 = atmti[i]
        for j in range(0, m1):
            if (na1[j]>=0):
                v2=atms[na1[j],1:4]
                #print(atms[na[j]], v2)
                nn =0 
                for lx in range(-1, 2):
                    for ly in range(-1,2):
                        for lz in range(-1,2):
                            w[0] = v2[0] + lx*latt[0][0]
                            w[1] = v2[1] + ly*latt[1][1]
                            w[2] = v2[2] + lz*latt[2][2]
                            dist_pbc[nn] = np.linalg.norm(v1-w)
                            nn = nn + 1
                dist[j] = np.min(dist_pbc)

        nba[i] = na1[np.argmin(dist)]
        na1[np.argmin(dist)] = -1
        #print('Min dist Ti-Ba is ', np.min(dist), dist, np.argmin(dist))
        dist = np.ones(nhalf)*10000
    #print(nba)
    for i in range(0, nti):
        atms[nba[i]][0] = 2
    return atms 



def rock_salt_order(nn, nb, nba, seed1, seed2, latt, atm):
    dmin = min_dist(nb, atm) + 1.0

    dist = np.zeros(nb)    
    atms_order = np.zeros((nb,4))
    dist_pbc = np.zeros(27)
    v2 = np.zeros(3)

    for i in range(0, nb):
        atms_order[i][0] = 1.0
        atms_order[i][1] = atm[i][0]
        atms_order[i][2] = atm[i][1]
        atms_order[i][3] = atm[i][2]

    sequence_a = [i for i in range(nb)]
    #print(sequence)
    # select a subset without replacement
    #subset = sample(sequence, nba)
    #print(subset)
    np.random.seed(randrange(0, 101, 2))
    #seed(seed1)
    shuffle(sequence_a)

    np.random.seed(randrange(100, 1001, 2))
    #seed(seed2)
    sequence_b = [i for i in range(nb)]
    #print(sequence_a, sequence_b)
    # select a subset without replacement
    #subset = sample(sequence, nba)
    #print(subset)
    shuffle(sequence_b)
    #ncount = 0 

    for i in range(0, nba):
        ma1 = sequence_a[i]
        ncount = 0
        if i > nb:
            print("please check concentration value", i)
            exit()
        if int(atms_order[sequence_a[i]][0]) == 1:
            v1 = atm[sequence_a[i]]
            for j in range(0, nb):
                mb1 = sequence_b[j]  
                m = 0
                for lx in range(-1, 2):
                    for ly in range(-1,2):
                        for lz in range(-1,2):
                            v2[0] = atm[mb1][0] + lx*latt[0][0]
                            v2[1] = atm[mb1][1] + ly*latt[1][1]
                            v2[2] = atm[mb1][2] + lz*latt[2][2]
                            dist_pbc[m] = np.linalg.norm(v1-v2)
                            m = m + 1
                d = np.min(dist_pbc) 
                if d > dmin:
                    if int(atms_order[ma1][0]) == 1:
                        #print('ok')
                        atms_order[ma1][0] = 3
                        ncount = 1
                       # print('ok, ma1')
                if ncount==1:
                    break
        else:
            nba = nba + 1

    #print(atms_order, len(atms_order))
    return atms_order


#----------------------------------------------------------
def ti_ba_neighbour_order(nbi, nti, latt, atmbi, atmti):
    dmin = min_ab_dist(nbi, nti, atmbi, atmti) + 1.0
  
    atms_order = np.zeros((nbi,4))
    for i in range(0, nbi):
        atms_order[i][0] = 1.0
        atms_order[i][1] = atmbi[i][0]
        atms_order[i][2] = atmbi[i][1]
        atms_order[i][3] = atmbi[i][2]

    ti_bi = np.zeros(nti*8)

    v2 = np.zeros(3)
    m=0
    for i in range(0, nti):
        v1=atmti[i]
        for j in range(0, nbi):
            v2= atmbi[j]
            dv = v1-v2
            dist = np.linalg.norm(dv)
            if dist <=dmin:
               #print(j) 
               ti_bi[m] = j
               m = m + 1
    ti_bi_uniq= np.zeros(m)
    for i in range(0,m):
        ti_bi_uniq[i] = ti_bi[i]
    ti_bi_uniq=np.unique(ti_bi_uniq)

    #print(ti_bi_uniq)
    shuffle(ti_bi_uniq)

    for i in range(0, nti):
        #print(ti_bi_uniq[i])
        atms_order[int(ti_bi_uniq[i])][0]=2.0
    return atms_order

def min_dist(n, atm):

    dist = np.zeros(n-1)
    
    for i in range(0,1):
        v1 = atm[i]
        m = 0
        for j in range(0, n):
            if j != i:
                #print(j)
                v2 = atm[j]
                dv = v1-v2
                dist[m] = np.linalg.norm(dv)
                m = m +1
                
            
    d = np.min(dist)
    return d

def min_ab_dist(nbi, nti, atmbi, atmti):
    
    dist=np.zeros(nbi)

    for i in range(0, 1):
        v1 = atmti[i]
        m = 0
        for j in range(0, nbi):
            v2 = atmbi[j]
            dv = v1-v2
            dist[m] = np.linalg.norm(dv)
            m = m + 1
    dmin=np.min(dist)
    print(dmin)
    return dmin 

