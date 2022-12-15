#using sparse matrices
import random
import numpy as np
from scipy import rand
import scipy.optimize
import time
import cvxopt


start = time.time()

def vickery_auction (bids):
    n = len(bids)
    if n == 0:
        return((-1,0))
    else:
        winner = 0
        first_price = bids[0]
        second_price = 0
        for i in range(1,n):
            if bids[i] >= first_price :
                second_price = first_price
                first_price = bids[i]
                winner = i
            elif bids[i] > second_price :
                second_price = bids[i]
        return((winner,second_price))


               
def vickery_with_reserves (bids,reserves):
    n = len(bids)
    n_res = len(reserves)
    if n!= n_res : 
        print("bids and reserves are not the same length")
        return("error")
    if n == 0:
        return((-1,0))
    else:
        new_bids = []
        index_new_bids = []
        for i in range (n):
            if bids[i] >= reserves[i]:
                new_bids.append(bids[i])
                index_new_bids.append(i)
        (temp_winner,temp_price_vickery) = vickery_auction(new_bids)
        if temp_winner != -1 :
            winner = index_new_bids[temp_winner]
            price = max(temp_price_vickery,reserves[winner])
        else : 
            winner = -1
            price = 0
        return(winner,price) 

def effect_of_reserves_on_bids(numerous_bids,reserves):
    n_rounds = len(numerous_bids)
    if n_rounds == 0 : return(0)  #revenue
    else :
        revenue = 0
        for i in range(n_rounds):
            (temp_winner,temp_price) = vickery_with_reserves(numerous_bids[i],reserves)
            revenue += temp_price
        return(revenue)


n_bidders = 3
n_rounds = 10



def generate_bids(n_bidders,n_rounds):
    list_of_bids = []
    for i in range (n_rounds):
        actual_bids = []
        for j in range (n_bidders):
            actual_bids.append(random.randint(0,100)) #can be replaced with a function in i and j
        list_of_bids.append(actual_bids)
    return list_of_bids

ex_set_of_bids = generate_bids(n_bidders,n_rounds)

#print(ex_set_of_bids)

#ex_set_of_bids = [[86, 30], [37, 3], [16, 51]]

for bids in ex_set_of_bids:
    bids.append(0)
    bids.append(0)

n_bidders +=2  

print(ex_set_of_bids)


def calculate_R (set_of_bids):
    R = [0,100_000_000]                  #we are supposed to add 0 and +infinity i will consider that 100_000_000 = +infinity
    for i in range (n_rounds):
        for j in range (n_bidders):
            if not (set_of_bids[i][j] in R):
                R.append(set_of_bids[i][j])
    R.sort()
    return R

ex_R = calculate_R (ex_set_of_bids)

print(ex_R)

ex_a = 0

def calculate_P_a (set_of_bids,a,R):
    P_a = []    #P_a is going to be a list of 4-tuple, those that are in P_a
    for b_1 in range (n_bidders):
        for b_2 in range (n_bidders):   # same notations for b_1 and b_2 as in Derakhshan's article (begining of page 7)
            if b_1 != b_2 :
                for r_1 in R:
                    for r_2 in R:
                        if set_of_bids[a][b_1] >= set_of_bids[a][b_2]:
                            if set_of_bids[a][b_1] >= r_1:
                                if set_of_bids[a][b_2] >= r_2:
                                    P_a.append((b_1,b_2,r_1,r_2))
    return P_a

ex_P_a =  calculate_P_a(ex_set_of_bids,ex_a,ex_R)

P_a_0 = calculate_P_a(ex_set_of_bids,0,ex_R)
P_a_1 = calculate_P_a(ex_set_of_bids,1,ex_R)
P_a_2 = calculate_P_a(ex_set_of_bids,2,ex_R)



def calculate_all_the_P_a (set_of_bids,R):
    all_the_P_a = []
    for a in range (n_rounds):
        all_the_P_a.append(calculate_P_a(set_of_bids,a,R))
    return (all_the_P_a)

ex_all_the_P_a = calculate_all_the_P_a (ex_set_of_bids,ex_R)


def calculate_number_of_profiles (all_the_P_a):
    number_of_profiles = 0
    for P_a in (all_the_P_a) :
        number_of_profiles += len(P_a)
    return number_of_profiles

def calculate_number_of_profiles_of_each_P_a (all_the_P_a):
    number_of_profiles_of_each_P_a = []
    for P_a in (all_the_P_a) :
        number_of_profiles_of_each_P_a.append(len(P_a))
    return number_of_profiles_of_each_P_a

ex_number_of_profiles_of_each_P_a = calculate_number_of_profiles_of_each_P_a (ex_all_the_P_a)

ex_number_of_profiles = calculate_number_of_profiles(ex_all_the_P_a)

def calculate_rev_a_p (set_of_bids,a,p):
    return (max (set_of_bids[a][p[1]],p[2]))

def calculate_rev_a_p_list (set_of_bids,all_the_P_a):     #it will output the revenue list, as a line
    rev_list = []
    for a in range (n_rounds):
        for p in all_the_P_a[a]:
            rev_list.append(calculate_rev_a_p(set_of_bids,a,p))
    return rev_list

ex_rev_a_p_list = calculate_rev_a_p_list(ex_set_of_bids,ex_all_the_P_a)

def calculate_rev_a_p_vector  (set_of_bids,all_the_P_a,R):
    rev_list = calculate_rev_a_p_list (set_of_bids, all_the_P_a)
    rev_vector_partial = np.array(rev_list)
    zeros_of_rev = np.zeros(n_bidders*len(R))
    rev_vector = np.concatenate((rev_vector_partial,zeros_of_rev))
    #rev_vector_trasposed = np.atleast_2d(rev_list)
    #rev_vector = np.transpose (rev_vector_trasposed)
    return rev_vector

def sp_rev (set_of_bids,all_the_P_a,R):
    return(cvxopt.matrix(calculate_rev_a_p_vector(set_of_bids,all_the_P_a,R)))

a = sp_rev(ex_set_of_bids,ex_all_the_P_a,ex_R)



ex_rev_vector = calculate_rev_a_p_vector(ex_set_of_bids,ex_all_the_P_a, ex_R)

#print(ex_rev_vector)


def calculate_I_0_0 (all_the_P_a):
    I_0_0_list = []
    number_of_profiles_of_each_P_a = calculate_number_of_profiles_of_each_P_a (all_the_P_a)
    for a in range(n_rounds):
        line_a = []
        for j in range (a):
            line_a += [0]*number_of_profiles_of_each_P_a[j]
        line_a += [1]*number_of_profiles_of_each_P_a[a]
        for j in range (a+1,n_rounds):
            line_a += [0]*number_of_profiles_of_each_P_a[j]
        I_0_0_list.append(line_a)
    return (np.array(I_0_0_list))

def sp_I_0_0 (all_the_P_a):
    return (cvxopt.sparse(cvxopt.matrix((calculate_I_0_0(all_the_P_a)))))

#ex_I_0_0 = sp_I_0_0(ex_all_the_P_a)

#print(ex_I_0_0)

def calculate_I_0_1 (R):
    return (np.zeros((n_rounds,n_bidders*len(R))))

def sp_I_0_1 (R):
    return(cvxopt.spmatrix([],[],[],(n_rounds,n_bidders*(len(R)))))

#ex_I_0_1 = calculate_I_0_1(ex_R)

#print(sp_I_0_1(ex_R))

def calculate_I_1_0 (all_the_P_a,R):
    I_1_0_list = []
    for a in range (n_rounds):
        for b in range (n_bidders):
            for r in R:
                line_b_r_a = []
                for j in range (a):
                    for p in all_the_P_a[j]: line_b_r_a.append(0)
                for p in all_the_P_a[a] :
                    if (p[0] == b and p[2] == r) or (p[1] == b and p[3] == r):
                        line_b_r_a.append(1)
                    else: line_b_r_a.append(0)
                for j in range (a+1,n_rounds):
                    for p in all_the_P_a[j]: line_b_r_a.append(0)
                I_1_0_list.append(line_b_r_a)
    return (np.array(I_1_0_list))

def sp_I_1_0 (all_the_P_a,R):
    values = []
    xindices = []
    yindices = []
    xindex = 0
    for a in range (n_rounds):
        for b in range (n_bidders):
            for r in R:
                yindex = a * n_bidders * len(R)
                for p in all_the_P_a[a] :
                    if (p[0] == b and p[2] == r) or (p[1] == b and p[3] == r):
                        values.append(1)
                        xindices.append(xindex)
                        yindices.append(yindex)
                    yindex += 1
                xindex += 1
    return(cvxopt.spmatrix(values,xindices,yindices,(n_bidders*n_rounds*len(R), calculate_number_of_profiles(all_the_P_a))))


#ex_I_1_0 = calculate_I_1_0(ex_all_the_P_a,ex_R)  # for this one I believe using sparse matrix would be better (this line takes some time,)

#print(sp_I_1_0(ex_all_the_P_a,ex_R))

#np.set_printoptions(threshold = sys.maxsize)
#print(ex_I_1_0)

def calculate_I_1_1 (R):
    I_1_1_list = []
    for a in range (n_rounds):
        for b in range (n_bidders):
            for r in R:
                line_b_r_a = []
                for bb in range (n_bidders):
                    for rr in R:
                        if (b == bb and r == rr):
                            line_b_r_a.append(-1)
                        else: line_b_r_a.append(0)
                I_1_1_list.append(line_b_r_a)
    return (np.array(I_1_1_list))

def sp_I_1_1 (R):
    partial_I_1_1 = []
    for a in range (n_rounds):
        partial_I_1_1.append(cvxopt.spdiag(cvxopt.matrix(np.ones(n_bidders*len(R)))))
    return(cvxopt.sparse([partial_I_1_1]))

#print(sp_I_1_1(ex_R))

        

#ex_I_1_1 = calculate_I_1_1(ex_R)

#print(ex_I_1_1)

def calculate_I (all_the_P_a,R):
    I_0_0 = calculate_I_0_0(all_the_P_a)
    I_0_1 = calculate_I_0_1(R)
    I_1_0 = calculate_I_1_0(all_the_P_a,R)
    I_1_1 = calculate_I_1_1(R)
    I_0 = np.concatenate((I_0_0,I_0_1),axis=1)
    I_1 = np.concatenate((I_1_0,I_1_1),axis=1)
    I = np.concatenate((I_0,I_1),axis=0)
    return I

# with the cvx solver that I am using I will need to add the non-negativity constraint on the sap because it is not automatically done by the solver
def sp_I_2_0 (all_the_P_a):
    return(cvxopt.spdiag(cvxopt.matrix((-1)*np.ones(calculate_number_of_profiles(all_the_P_a)))))

def sp_I_2_1 (all_the_P_a,R):
    return(cvxopt.spmatrix([],[],[],(calculate_number_of_profiles(all_the_P_a),len(R)*n_bidders)))

def sp_I (all_the_P_a,R):
    I_0_0 = sp_I_0_0(all_the_P_a)
    I_0_1 = sp_I_0_1(R)
    I_1_0 = sp_I_1_0(all_the_P_a,R)
    I_1_1 = sp_I_1_1(R)
    I_2_0 = sp_I_2_0(all_the_P_a)
    I_2_1 = sp_I_2_1(all_the_P_a,R)
    return (cvxopt.sparse([[I_0_0,I_1_0,I_2_0],[I_0_1,I_1_1,I_2_1]]))

#print(sp_I(ex_all_the_P_a,ex_R))

timebeforecalI = time.time()
#ex_I = calculate_I(ex_all_the_P_a,ex_R)
timeaftercalI = time.time()

def calculate_b (R):
    b_0 = np.ones((n_rounds))
    b_1 = np.zeros((n_rounds * n_bidders * len (R)))
    b = np.concatenate((b_0,b_1))
    return b

def sp_b (all_the_P_a,R):
    b_0 = np.ones(n_rounds)
    b_1 = np.ones((n_rounds * n_bidders * len (R)))
    b_2 = np.zeros(calculate_number_of_profiles(all_the_P_a))
    return (cvxopt.matrix(np.concatenate((b_0,b_1,b_2))))

#print(sp_b(ex_all_the_P_a,ex_R))

#ex_b = calculate_b(ex_R)

def calculate_E (R,number_of_profiles):
    E_0 = np.zeros((n_bidders,number_of_profiles))
    list_E_1 = []
    for i in range (n_bidders):
        line_E_1 = []
        for b in range (i):
            for r in range (len(R)):
                line_E_1.append(0)
        for r in range (len(R)):
            line_E_1.append(1)
        for b in range (i+1,n_bidders):
            for r in range(len(R)):
                line_E_1.append(0)
        list_E_1.append(line_E_1)
    E_1 = np.array(list_E_1)
    E = np.concatenate((E_0,E_1),axis = 1)
    return E

def sp_E (R,number_of_profiles):
    E_0 = cvxopt.spmatrix([],[],[],(n_bidders,number_of_profiles))
    values = []
    xindices = []
    yindices = []
    for i in range (n_bidders):
        for r in range (len(R)):
            values.append(1)
            xindices.append(i)
            yindices.append(r + i*len(R))
    E_1 = cvxopt.spmatrix(values,xindices,yindices,(n_bidders,n_bidders*len(R)))
    return (cvxopt.sparse([[E_0],[E_1]]))

#print (sp_E(ex_R,ex_number_of_profiles))
        

#ex_E = calculate_E(ex_R,ex_number_of_profiles)

#print(ex_E)

def calculate_c ():
    return np.ones(n_bidders)

def sp_c ():
    return(cvxopt.matrix(np.ones(n_bidders)))

#print(sp_c())

#ex_c = calculate_c()

beforelinprog = time.time()

#res = scipy.optimize.linprog((-1)*ex_rev_vector, A_ub=ex_I, b_ub=ex_b, A_eq=ex_E, b_eq=ex_c, bounds=(0,None), method='highs', callback=None, options=None, x0=None)

sp_res = cvxopt.solvers.lp(c=sp_rev(ex_set_of_bids,ex_all_the_P_a,ex_R),G=sp_I(ex_all_the_P_a,ex_R),h=sp_b(ex_all_the_P_a,ex_R),A=sp_E(ex_R,ex_number_of_profiles),b=sp_c())

afterlinprog = time.time()

#print((-1)*res.fun)

print(sp_res['primal objective'])

#print(res.x)

# with np.printoptions(threshold=np.inf): print(res.x)

def prolpr (res,R,set_of_bids,all_the_P_a):
    n_all_the_P_a = calculate_number_of_profiles(all_the_P_a)
    s_q = res.x
    q = s_q[n_all_the_P_a:]
    n_R = len(R)
    reserve=[]    # will corespond to r^r of the pro-LPr algo
    for i in range (n_bidders):
        u = np.random.rand()
        s = 0
        b = False
        for ind_r in range (n_R) :
            if not b:
                s+= q[n_R*i + ind_r]
                if s >= u :
                    reserve.append(R[ind_r])
                    b = True
        if not b : print("probleme")
    r_0 = [0]*n_bidders
    rev_r_0 = effect_of_reserves_on_bids(set_of_bids,r_0)
    rev_reserve = effect_of_reserves_on_bids(set_of_bids,reserve)
    if rev_r_0 >= rev_reserve :
        return (rev_r_0,r_0)
    else:
        return(rev_reserve,reserve)

#esp_star,reserves_porlpr =  prolpr(res,ex_R,ex_set_of_bids,ex_all_the_P_a)


"""
print(reserves_porlpr)
    


print(effect_of_reserves_on_bids(ex_set_of_bids,reserves=reserves_porlpr))

print(effect_of_reserves_on_bids(ex_set_of_bids,reserves=([0]*n_bidders)))

print(beforelinprog - start)
print(afterlinprog - beforelinprog)

print(timeaftercalI -timebeforecalI)
"""