"""NetDes-Duo - paired gene regulatory network inference from single-cell
trajectories: Hill-function ODE fitting plus edge selection over one network
(NetDes) or two conditions jointly (NetDes-Duo).
"""

from __future__ import print_function
from itertools import combinations
import numpy as np
import itertools
import random
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import least_squares
from scipy.optimize import minimize
from scipy.optimize import Bounds
from scipy.optimize import fmin_l_bfgs_b
import matplotlib
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from sklearn.model_selection import KFold
from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report
from sklearn.svm import SVC
from scipy import interpolate
np.set_printoptions(suppress=True)
np.set_printoptions(precision=3)
import math
import copy
import warnings
from matplotlib.backends.backend_pdf import PdfPages




def modelHill(pars, u):# pars=(lambda, 5*n+1, r0/max(R)) here u=inp_Signal/max(inp_signal), n is Uniform (0,1)take 0.001<x[0]=xbasal/max(inp_signal)<=0.999
    #return np.log2((pars[0] + (1-pars[0])/(1+(u/pars[2])**(5*(pars[1])+1)))/(1+pars[0]))
    return np.log2((pars[0] + (1-pars[0])/(1+(u/pars[2])**(5*(pars[1])+1))))


def regulate_term(x):
    num=round(len(x)/3)
    kkm=[]
    for i in range(num):
        kkm.append(abs(x[3*i]-1))
        results=np.sum(kkm)
    return(results)




def fmodelOneHilltest(x, *args):
    y = args[-1] # The observed data is always the last argument
    N = y.size

    # Iterate over the input arguments (excluding the observed data)
    # and use the Hill equation to calculate the model for each input
    mods = []
    for i, a in enumerate(args[:-1]):
        xm = x[i*3:(i+1)*3] # Extract the relevant values from x for this input
        mods.append(modelHill(xm, a))
    # Sum up the models for all inputs and calculate the MSE
    mod_pred = x[-1] + sum(mods)
    resid = (mod_pred - y)**2
    func_val = sum(resid)
    mse = 1/N * func_val

    return mse



def fit_results(x, *args):

    # Iterate over the input arguments (excluding the observed data)
    # and use the Hill equation to calculate the model for each input
    mods = []
    for i, a in enumerate(args):
        xm = x[i*3:(i+1)*3] # Extract the relevant values from x for this input
        mods.append(modelHill(xm, a))

    # Sum up the models for all inputs and calculate the MSE
    mod_pred = x[-1] + sum(mods)

    return mod_pred


def fmodelOneHilltestw(x, *args):
    y = args[-2]
    w = args[-1] # The observed data and w is always the last two argument
    N = y.size

    # Iterate over the input arguments (excluding the observed data)
    # and use the Hill equation to calculate the model for each input
    mods = []
    k=[]
    for i, a in enumerate(args[:-2]):
        xm = x[i*3:(i+1)*3] # Extract the relevant values from x for this input
        km = w*((x[i*3]-1)**2)
        mods.append(modelHill(xm, a))
        k.append(km)

    # Sum up the models for all inputs and calculate the MSE
    mod_pred = x[-1] + sum(mods)  ## log2(x[-1]) to x[-1]
    resid = (mod_pred - y)**2
    func_val = sum(resid)
    mse = (1/N * func_val)+sum(k)

    return mse


def get_args(gene,ob,target,w):# gene is a dataframe_tf12  ob is the name of the regulators
    argsss=[]
    for i in range(len(ob)):
        resu=np.array(gene.loc[ob[i]].iloc[ 0:],dtype='float64')
        argsss.append(resu.flatten())
    resu2=np.array(gene.loc[target].iloc[ 0:],dtype='float64')
    argsss.append(resu2.flatten())
    for i in range(len(argsss)-1):
        #min_val = np.min(argsss[i])
        #argsss[i] = argsss[i] - min_val
        argsss[i] = argsss[i]/max(argsss[i])
    argsss.append(w)
     # Convert argsss to a tuple and return it
    argsss_tuple = tuple(argsss)
    return argsss_tuple



def get_args2(gene,ob,target):# gene is a dataframe_tf12  ob is the name of the regulators
    argsss=[]
    for i in range(len(ob)):
        resu=np.array(gene.loc[ob[i]].iloc[ 0:],dtype='float64')
        argsss.append(resu.flatten())
    resu2=np.array(gene.loc[target].iloc[ 0:],dtype='float64')
    argsss.append(resu2.flatten())
    for i in range(len(argsss)-1):
        #min_val = np.min(argsss[i])
        #argsss[i] = argsss[i] - min_val
        argsss[i] = argsss[i]/max(argsss[i])
     # Convert argsss to a tuple and return it
    argsss_tuple = tuple(argsss)
    return argsss_tuple



def stPoints(n): # n should be the number of the regulators
    # generate the starting point of the minimizer by random.
    st=[]
    for i in range(n):
        a1 = random.uniform(0,2)
        a2 = random.uniform(0.3,0.7)
        a3 = random.uniform(0.01,0.99)
        st.append(a1)
        st.append(a2)
        st.append(a3)
    b=random.uniform(-2,2)
    st.append(b)
    return(st)


def generate_same_range(point):
    n = (len(point) - 1) // 3
    same_range_point = []
    for j in range(n):
        a1 = point[3*j]
        if a1 >= 2:
            a1 = random.uniform(2, 4)
        else:
            a1 = random.uniform(0.25, 0.5)
        a2 = random.uniform(0.3, 0.7)
        a3 = random.uniform(0.3,0.7)
        same_range_point.append(a1)
        same_range_point.append(a2)
        same_range_point.append(a3)
    b = random.uniform(-2, 2)
    same_range_point.append(b)
    return same_range_point



def getarray(x,y,w):## x should be the results of the kf.split
    arr=[]
    for i in range(len(x.T)):
        arr.append(x.T[i])
    arr.append(y)
    arr.append(w)
    arr_tuple = tuple(arr)

    return arr_tuple



def getarray2(x,y):## x should be the results of the kf.split
    arr=[]
    for i in range(len(x.T)):
        arr.append(x.T[i])
    arr.append(y)
    arr_tuple = tuple(arr)

    return arr_tuple



def getarray3(x):## x should be the results of the kf.split
    arr=[]
    for i in range(len(x.T)):
        arr.append(x.T[i])
    arr_tuple = tuple(arr)

    return arr_tuple




def checkstPoints(St):
    rs=[]
    for i in range(int(len(St)/3)):
        if St[3*i]>1:
            rs=rs+[1]
        else:
            rs=rs+[0]
    return rs



def lamlist(St):
    lam=[]
    for i in range(int(len(St)/3)):
        lam=lam+[10**St[3*i]]
    return lam


def get_combinations(my_list):
    combinations = []
    for i in range(1, len(my_list) + 1):
        for subset in itertools.combinations(my_list, i):
            combinations.append(list(subset))
    return(combinations)


def spresults_find(mselist,roclist):
    pc=0
    i=0
    counts=mselist.sort_values(by='mean').iloc[:,7].value_counts()
    rocint2= counts.idxmax()
    while pc==0 and i==0:
        if mselist.sort_values(by='mean').iloc[i,7]==rocint2:
            pc=i
            break
        else:
            i=i+1
    return(mselist.sort_values(by='mean').iloc[i,6])

def wfind(yfit):
    yfit1=yfit.copy()
    yfit1=np.append(yfit,0)
    yfit2=np.insert(yfit, 0, 0)
    yslop=yfit1-yfit2
    yslop1=np.append(yslop,0)
    yslop2=np.insert(yslop, 0, 0)
    ytest=yslop1-yslop2
    slop=(max(yfit)-min(yfit))/100
    cutslop=0.33*slop
    ## work for the 2
    wp=[]
    for i in range(len(yslop)-2):
        if yslop[i+1]>cutslop and ytest[i+2]>0:
            finding=i+1
            wp.append(finding)
    wp2=np.append(wp,0)
    wp1=np.insert(wp, 0, 0)
    choosetest=wp2-wp1>1
    positions = np.where(choosetest)[0]
    if len(positions)==0:
        return(0)
    else:
        i=wp[positions[np.argmin(yfit[wp2[positions]])]]
        if (yfit1[i+1]-min(yfit1))/(max(yfit1)-min(yfit))>0.25:
            return(0)
        else:
            return(i)


def influence(res):
    lambd=[]
    for i in range(round(len(res)/3)):
            lambd.append(res[3*i:(3*i+3)])
    return(lambd)


def remove_groups(lst):
    for i in  range(len(parameters)-4, -1, -3):
        if lst[i] == 0:
            del lst[i:i+3]
    return lst


def influence2(res):
    lambd=[]
    for i in range(len(res)):
        if i % 3 == 0:
            lambd.append(res[i])
    return(lambd)


def power_of_two(x):
    return 2 ** x


def transtolog2(x):
    return np.log2(x)


def search_list(list1, list2,x1,x2):
    result = []
    for i in range(len(list1)):
        if list1[i]>x1 and list1[i]<x2:
            result.append(list2[i])
    return result


def sum_lists(list1, list2,w):
    """
    Takes two lists of equal length as inputs and returns a new list that contains the sum of their corresponding elements.
    """
    if len(list1) != len(list2):
        print("Error: lists are not of equal length")
        return None
    else:
        return [x + w*y for x, y in zip(list1, list2)]


def stpoints_nums(MSES,n):
    stdatabases_val=[]
    ii=0
    if len(MSES)>100:
        for ii3 in range(50):
            test_st=np.array(MSES.sort_values(by='mean').iloc[[ii3],[6]])[0][0]
            for i in range(n):
                stdatabases_val.append(generate_same_range(test_st))
    else:
        for ii3 in range(int(len(MSES)/2)):
            test_st=np.array(MSES.sort_values(by='mean').iloc[[ii3],[6]])[0][0]
            for i in range(n):
                stdatabases_val.append(generate_same_range(test_st))
    return(stdatabases_val)



def boundary(t,act=[],rep=[]): # get the boundary
    unit=[(0.01, 100), (0, 1), (0.01, 0.99)]
    last=[(-6,6)]
    bou=[]
    for i in range(t):
        bou=bou+unit
    bou=bou+last
    for i2 in act:
        bou[i2*3]=(1.00001,100)
    for i2 in rep:
        bou[i2*3]=(0.01,0.99999)
    return(bou)


def stPoints2(n, act=[], rep=[]):
    points = []
    for i in range(2 ** n):
        # Check if any bits in 'act' should be set to 0
        for j in act:
            if i & (1 << j):
                break
        else:
            # Check if any bits in 'rep' should be set to 1
            for j in rep:
                if not (i & (1 << j)):
                    break
            else:
                st = []
                for j in range(n):
                    a1 = random.uniform(0.25, 0.5) if i & (1 << j) else random.uniform(2, 4)
                    a2 = random.uniform(0.3, 0.7)
                    a3 = random.uniform(0.01,0.99)
                    st.extend([a1, a2, a3])
                b = random.uniform(-2, 2)
                st.append(b)
                points.append(st)
    return points





def CV_fit(ob, t_array, stdatabase, boundary_cv,cc=5,w=0):
    X=t_array[0:len(ob)]
    X=np.stack((X),axis=1)
    y=t_array[len(ob)]
    kf = KFold(n_splits=cc)
    nd=0
    MSE_test_term=pd.DataFrame(index=range(len(stdatabase)), columns=range(5))
    for train_index, test_index in kf.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        args1=getarray(X_train,y_train,w)
        argss=getarray2(X_train,y_train)
        argst=getarray3(X_test)
        for th in range(len(stdatabase)):
            stPoint=stdatabase[th]
            resCl=minimize(fmodelOneHilltestw,stPoint, args=args1, method='L-BFGS-B', jac='2-point', bounds=boundary_cv)
            y_pred=fit_results(resCl.x,*argst)
            MSE_test_term.iloc[th,nd]=mean_squared_error(y_test,y_pred)
        nd=nd+1
    return(MSE_test_term)


def MSE_stPoint_one_gene(gene_list, network, expression, log_expression, ii, pick = "all"):
    ## running the gene one by one
    choose=gene_list[ii]
    ob = list(network[network['Target']==choose]['Source'])
    if choose in ob:
        ob.remove(choose)
    tfexpression_use=expression.copy()
    tfexpression_use.loc[choose]=log_expression.loc[choose]
    #### Force
    restri_A = list(network[(network['Target'] == choose) & (network['Interaction'] == 2)]['Source'])
    restri_R = list(network[(network['Target'] == choose) & (network['Interaction'] == 0)]['Source'])
    positions_A = [ob.index(item) for item in restri_A]
    positions_R = [ob.index(item) for item in restri_R]
    random.seed(1)
    nofeasB=boundary(len(ob),act=positions_A,rep=positions_R)
    test1=get_args2(tfexpression_use,ob=ob,target=gene_list[ii]) # get the initial array
    stdatabase=stPoints2(len(ob),act=positions_A,rep=positions_R)# sampling SP
    if pick != 'all' and 2**len(ob)>=500:
        a=int(pick*2**len(ob))
        random_sample = random.sample(range(2**len(ob)), a)
        stdatabase=[stdatabase[i] for i in random_sample]
    MSE_test_term= CV_fit(ob=ob,t_array=test1,stdatabase=stdatabase,boundary_cv=nofeasB)
    means=MSE_test_term.mean(axis=1)
    MSE_test_term['mean']=means

    return stdatabase, MSE_test_term

def MSE_stPoint(gene_list,network,expression,log_expression,pick='all'):
    """
    This function computes the Mean Squared Error (MSE) for various sampled starting points (SP) 
    for each gene in a given list. The goal is to evaluate how well the sampled starting points 
    fit the gene expression data based on a given regulatory network and expression matrix.

    Parameters:
    ----------
    gene_list : list
        List of target genes for which MSE will be calculated.

    network : pandas.DataFrame
        A dataframe representing the regulatory network. It contains interactions where each 
        row specifies a source gene, a target gene, and their interaction type. The 'Interaction' 
        column should have values like 2 (activation), 0 (repression) or 1 unknown.

    expression : pandas.DataFrame
        A dataframe of scaled gene expression data where rows represent genes, and columns 
        represent samples or conditions.

    log_expression : pandas.DataFrame
        A dataframe of scaled gene expression data in log scale, where rows represent genes, 
        and columns represent samples or conditions.

    pick : str or float, optional, default 'all'
        Specifies how many starting points to sample. If 'all', all possible starting points 
        will be used. If a float between 0 and 1, it determines the proportion of starting 
        points to sample, with random sampling applied.

    Returns:
    --------
    list :
        A list containing:
        - all_st: A list of starting points for each gene in `gene_list`.
        - MSE_st: A list of dataframes, where each dataframe contains MSE values calculated 
          for each sampled starting point for the corresponding gene.

    """
    all_st=[]
    MSE_st=[]
    for ii in range(len(gene_list)):
        ## running the gene one by one
        choose=gene_list[ii]
        print(choose)
        ob = list(network[network['Target']==choose]['Source'])
        if choose in ob:
            ob.remove(choose)
        tfexpression_use=expression.copy()
        tfexpression_use.loc[choose]=log_expression.loc[choose]
        #### Force
        restri_A = list(network[(network['Target'] == choose) & (network['Interaction'] == 2)]['Source'])
        restri_R = list(network[(network['Target'] == choose) & (network['Interaction'] == 0)]['Source'])
        positions_A = [ob.index(item) for item in restri_A]
        positions_R = [ob.index(item) for item in restri_R]
        random.seed(1)
        nofeasB=boundary(len(ob),act=positions_A,rep=positions_R)
        test1=get_args2(tfexpression_use,ob=ob,target=gene_list[ii]) # get the initial array
        stdatabase=stPoints2(len(ob),act=positions_A,rep=positions_R)# sampling SP
        if pick != 'all' and 2**len(ob)>=500:
            a=int(pick*2**len(ob))
            random_sample = random.sample(range(2**len(ob)), a)
            stdatabase=[stdatabase[i] for i in random_sample]
        all_st.append(stdatabase)# save stpoint
        MSE_test_term= CV_fit(ob=ob,t_array=test1,stdatabase=stdatabase,boundary_cv=nofeasB)
        means=MSE_test_term.mean(axis=1)
        MSE_test_term['mean']=means
        MSE_st.append(MSE_test_term)
    for k in range(len(gene_list)):
        MSE_st[k]['SP']=all_st[k]
    return [all_st, MSE_st]



def CV_tuning(ob,t_array,nums,stdatabase,boundary_cv,cc=5):
    X=t_array[0:len(ob)]
    X=np.stack((X),axis=1)
    y=t_array[len(ob)]
    kf = KFold(n_splits=cc)
    nd=0
    logwstep=np.linspace(start=-15,stop=2, num=nums-1,dtype='float64')
    wstep=[0]*nums
    wstep[1:nums]=np.exp(logwstep)
    MSE_min=[0]*len(wstep)
    for i in range(len(wstep)):
        w=wstep[i]
        MSE_test_term=pd.DataFrame(index=range(len(stdatabase)), columns=range(cc))
        nd=0
        for train_index, test_index in kf.split(X):
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]
            args1=getarray(X_train,y_train,w)
            argss=getarray2(X_train,y_train)
            argst=getarray3(X_test)
            for th in range(len(stdatabase)):
                stPoint=stdatabase[th]
                resCl=minimize(fmodelOneHilltestw,stPoint, args=args1, method='L-BFGS-B', jac='2-point', bounds=boundary_cv)
                y_pred=fit_results(resCl.x,*argst)
                MSE_test_term.iloc[th,nd]=mean_squared_error(y_test,y_pred)
            nd=nd+1
        means=MSE_test_term.mean(axis=1)
        MSE_min[i]=min(means)
    return MSE_min


def MSE_top_one_gene(gene_list,network,expression,log_expression,MSE_table,ii,cc=5,num_tuning=20):
    choose=gene_list[ii]
    ob = list(network[network['Target']==choose]['Source'])
    if choose in ob:
        ob.remove(choose)
    tfexpression_use=expression.copy()
    tfexpression_use.loc[choose]=log_expression.loc[choose]
    #### 1. choose the regulators gene   ### get a series of stPoints
    restri_A = list(network[(network['Target'] == choose) & (network['Interaction'] == 2)]['Source'])
    restri_R = list(network[(network['Target'] == choose) & (network['Interaction'] == 0)]['Source'])
    positions_A = [ob.index(item) for item in restri_A]
    positions_R = [ob.index(item) for item in restri_R]
    random.seed(1)
    nofeasB=boundary(len(ob),act=positions_A,rep=positions_R)
    test1=get_args2(tfexpression_use,ob=ob,target=gene_list[ii]) # get the initial array
    #stdatabase=stpoint(len(ob))
    MSE_all=MSE_table[1]
    stdatabase=stpoints_nums(MSE_all[ii],10)
    MSE_test_term= CV_fit(ob=ob,t_array=test1,stdatabase=stdatabase,boundary_cv=nofeasB)
    means=MSE_test_term.mean(axis=1)
    MSE_test_term['mean']=means
    MSE_test_term['mean'] = MSE_test_term['mean'].astype(float)
    stnumber=MSE_test_term.nsmallest(4, 'mean').index # get the row number of the stPoint database
    ## New stPoint dataset
    stnumber=stnumber.tolist()
    selected_stpoint = [stdatabase[i] for i in stnumber]
    test1=get_args(tfexpression_use,ob=ob,target=gene_list[ii],w=0)
    stPoint=selected_stpoint[0]
    resCl=minimize(fmodelOneHilltestw,stPoint, args=test1, method='L-BFGS-B', jac='2-point', bounds=nofeasB)
    MSE_min=CV_tuning(ob=ob,t_array=test1,nums=num_tuning,stdatabase=selected_stpoint,boundary_cv=nofeasB,cc=cc)

    return selected_stpoint, stdatabase, resCl.x, MSE_min, MSE_test_term




def MSE_top(gene_list,network,expression,log_expression,MSE_table,cc=5,num_tuning=20):
    """
    This function refines the selection of starting points for each gene in a list by selecting 
    the top MSE values and generating new starting points for further optimization. It returns 
    refined starting points, fitted parameters, and minimized MSE values for each gene.

    Parameters:
    ----------
    gene_list : list
        List of target genes for which refined MSE calculations will be performed.
    
    network : pandas.DataFrame
        A dataframe representing the regulatory network. It contains interactions where each 
        row specifies a source gene, a target gene, and their interaction type. The 'Interaction' 
        column should have values like 2 (activation) or 0 (repression).
    
    expression : pandas.DataFrame
        A dataframe of scaled gene expression data where rows represent genes, and columns 
        represent samples or conditions.
    
    log_expression : pandas.DataFrame
        A dataframe of scaled gene expression data in log scale, where rows represent genes, 
        and columns represent samples or conditions.
    
    MSE_table : list
        A list where each element corresponds to a gene and contains a table of MSE values 
        for previously calculated starting points.
    
    cc : int, optional, default 5
        The number of cross-validation folds to use for tuning the model.
    
    num_tuning : int, optional, default 20
        The number of tuning iterations to perform when minimizing MSE.

    Returns:
    --------
    list :
        A list containing:
        - selected_total: A list of selected starting points for each gene.
        - all_st2: A list of newly generated starting points for each gene.
        - res_par: A list of optimized parameter values obtained by minimizing MSE for each gene.
        - MSE_min_total: A list of minimized MSE values for each gene after further tuning.
        - MSE_st2: A list of updated dataframes for each gene, containing MSE values for 
          the newly generated starting points and their corresponding means.
    """
    selected_total=[]
    all_st2=[]## *10
    res_par=[]
    MSE_min_total=[]
    MSE_st2=[]
    for ii in range(len(gene_list)):
        # work flow
        #### try to write down a work flow
        # the number of gene # K fold
        choose=gene_list[ii]
        print(choose)
        ob = list(network[network['Target']==choose]['Source'])
        if choose in ob:
            ob.remove(choose)
        tfexpression_use=expression.copy()
        tfexpression_use.loc[choose]=log_expression.loc[choose]
        #### 1. choose the regulators gene   ### get a series of stPoints
        restri_A = list(network[(network['Target'] == choose) & (network['Interaction'] == 2)]['Source'])
        restri_R = list(network[(network['Target'] == choose) & (network['Interaction'] == 0)]['Source'])
        positions_A = [ob.index(item) for item in restri_A]
        positions_R = [ob.index(item) for item in restri_R]
        random.seed(1)
        nofeasB=boundary(len(ob),act=positions_A,rep=positions_R)
        test1=get_args2(tfexpression_use,ob=ob,target=gene_list[ii]) # get the initial array
        #stdatabase=stpoint(len(ob))
        MSE_all=MSE_table[1]
        stdatabase=stpoints_nums(MSE_all[ii],10)
        all_st2.append(stdatabase)
        MSE_test_term= CV_fit(ob=ob,t_array=test1,stdatabase=stdatabase,boundary_cv=nofeasB)
        means=MSE_test_term.mean(axis=1)
        MSE_test_term['mean']=means
        MSE_test_term['mean'] = MSE_test_term['mean'].astype(float)
        MSE_st2.append(MSE_test_term)
        stnumber=MSE_test_term.nsmallest(4, 'mean').index # get the row number of the stPoint database
        ## New stPoint dataset
        stnumber=stnumber.tolist()
        selected_stpoint = [stdatabase[i] for i in stnumber]
        selected_total.append(selected_stpoint)
        test1=get_args(tfexpression_use,ob=ob,target=gene_list[ii],w=0)
        stPoint=selected_stpoint[0]
        resCl=minimize(fmodelOneHilltestw,stPoint, args=test1, method='L-BFGS-B', jac='2-point', bounds=nofeasB)
        res_par.append(resCl.x)
        MSE_min=CV_tuning(ob=ob,t_array=test1,nums=num_tuning,stdatabase=selected_stpoint,boundary_cv=nofeasB,cc=cc)
        MSE_min_total.append(MSE_min)
    for ct in range(len(gene_list)):
        MSE_st2[ct]['SP'] = all_st2[ct]
    return [selected_total,all_st2,res_par,MSE_min_total,MSE_st2]



def tuning(MSE_table2,nums,gene_list):
    """
    This function performs parameter tuning by calculating optimal weight parameters (w) 
    based on the MSE values obtained for a list of genes. It uses spline interpolation 
    to fit the MSE data and determine the optimal weights for further optimization.

    Parameters:
    ----------
    MSE_table2 : list
        A list where each element corresponds to a gene and contains MSE values 
        for tuning iterations. The 4th element of this list (index 3) should contain 
        the minimized MSE values (`MSE_min_total`) for each gene.
    
    nums : int
        The number of weight tuning steps to be used, i.e., how many different 
        weight values to try when optimizing, which should same to num_tuning in function MSE_top
    
    gene_list : list
        List of genes for which the tuning parameters will be calculated.

    Returns:
    --------
    w2_tot : list
        A list containing the optimal weight parameters (`w`) for each gene in `gene_list`.

    """
    w2_tot=[]
    logwstep=np.linspace(start=-15,stop=2, num=nums-1,dtype='float64')
    wstep=[0]*nums
    wstep[1:nums]=np.exp(logwstep)
    for ii in range(len(gene_list)):
        MSE_min_total=MSE_table2[3]
        y = np.log(MSE_min_total[ii])#par_error
        x=range(0,len(y))
        knot_numbers = 2
        x_new = np.linspace(0, 1, knot_numbers+2)[1:-1]
        q_knots = np.quantile(x, x_new)
        t,c,k = interpolate.splrep(x, y, t=q_knots, s=1)
        yfit = interpolate.BSpline(t,c,k)(x)
        w2=wstep[wfind(yfit=yfit)]
        w2_tot.append(w2)
    return w2_tot


def MSE_bin(MSE_table2,w_tot,gene_list,network,expression,log_expression):
    """
    This function refines the model by applying binning to the previously computed MSE results 
    and further optimizes the parameter space for each gene in the provided gene list. 
    It performs iterative optimization for different starting points and logs the results.

    Parameters:
    ----------
    MSE_table2 : list
        A list of data where each element corresponds to a gene and contains MSE values 
        for previously calculated starting points. The 5th element (index 4) should contain 
        the MSE table (`MSE_st2`) to be updated.
    
    w_tot : list
        A list of tuning weight steps (`w2`) for each gene, previously calculated by the `tuning` function.
    
    gene_list : list
        List of genes for which further MSE refinement will be performed.
    
    network : pandas.DataFrame
        A dataframe representing the regulatory network. It contains interactions where each 
        row specifies a source gene, a target gene, and their interaction type.
    
    expression : pandas.DataFrame
        A dataframe of scaled gene expression data where rows represent genes, and columns 
        represent samples or conditions.
    
    log_expression : pandas.DataFrame
        A dataframe of scaled gene expression data in log scale, where rows represent genes, 
        and columns represent samples or conditions.

    Returns:
    --------
    MSE_table2 : list
        The updated `MSE_table2` with an additional element containing the selected 
        optimized parameters (`res_pick_all`) for each gene. It also updates the MSE table 
        with bin information (`stbinss`).

    """
    MSE_st2=copy.deepcopy(MSE_table2[4])
    res_w_all=[]
    lam_w_all=[]
    for ii in range(len(gene_list)):
        w2=w_tot[ii]
        choose=gene_list[ii]
        print(choose)
        ob = list(network[network['Target']==choose]['Source'])
        if choose in ob:
            ob.remove(choose)
        tfexpression_use=expression.copy()
        tfexpression_use.loc[choose]=log_expression.loc[choose]
        restri_A = list(network[(network['Target'] == choose) & (network['Interaction'] == 2)]['Source'])
        restri_R = list(network[(network['Target'] == choose) & (network['Interaction'] == 0)]['Source'])
        positions_A = [ob.index(item) for item in restri_A]
        positions_R = [ob.index(item) for item in restri_R]
        nofeasB=boundary(len(ob),act=positions_A,rep=positions_R)
        test1=get_args(tfexpression_use,ob=ob,target=gene_list[ii],w=w2)
        res_w=[]
        lam_w=[]
        for i in range(len(MSE_st2[ii])):
            stPoint=MSE_st2[ii].sort_values(by='mean').iloc[i,6]
            resCl=minimize(fmodelOneHilltestw,stPoint, args=test1, method='L-BFGS-B', jac='2-point', bounds=nofeasB)
            res_w.append(resCl.x)
            lam_w.append(lamlist(resCl.x))
        res_w_all.append(res_w)
        lam_w_all.append(lam_w)
    stbinss=[]
    for ii in range(len(gene_list)):
        stbins=[]
        for i in range(len(res_w_all[ii])):
            stbins.append(checkstPoints(res_w_all[ii][i]))
        stbinss.append(stbins)
    res_pick_all=[]
    for ii in range(len(gene_list)):
        res_pick=[]
        for i in range(10):
            res_pick.append(res_w_all[ii][i])
        res_pick_all.append(res_pick)
    for ct in range(len(gene_list)):
            #MSE_st2[ct]['SP'] = all_st2[ct]
        MSE_st2[ct]['bin']=stbinss[ct]
    MSE_table2[4]=MSE_st2
    MSE_table2.append(res_pick_all)
    return MSE_table2


def consis_test(gene_list,MSE_table):
    """
    This function tests the consistency of starting points across different genes by 
    calculating an ROC-like score for the top MSE values. It returns the consistency 
    scores (`roc`) ( 1 is 100% starting points fitted activation) and a list of starting points for each gene.

    Parameters:
    ----------
    gene_list : list
        List of genes for which the consistency test will be performed.
    
    MSE_table : list
        A list of data where each element corresponds to a gene and contains MSE values 
        for previously calculated starting points. The 5th element (index 4) should contain 
        the MSE table (`MSE_st2`), which holds the sorted MSE values and associated starting points.
        Mainly obtained from the results MSE_bin

    Returns:
    --------
    list :
        A list containing two elements:
        - roc: A list where each element is a consistency score (ROC-like) for each gene 
          based on the top MSE values.
        - st_roc: A list of starting points for each gene, where each element corresponds 
          to the top-ranked starting points based on the MSE values.

    """
    MSE_st2=MSE_table[4]
    roc=[]
    for ii in range(len(gene_list)):
        lenall=[0]*len(MSE_st2[ii].sort_values(by='mean').iloc[0,7])
        if len(MSE_st2[ii])>200:
            for i in range(100):
                stb=MSE_st2[ii].sort_values(by='mean').iloc[i,7]
                lenall=sum_lists(lenall,stb,1)
            roc_new = [x/100 for x in lenall]
            roc.append(roc_new)
        else:
            for i in range(int(len(MSE_st2[ii])/2)):
                itt=len(MSE_st2[ii])/2
                stb=MSE_st2[ii].sort_values(by='mean').iloc[i,7]
                lenall=sum_lists(lenall,stb,1)
            roc_new = [x/itt for x in lenall]
            roc.append(roc_new )
    ### get smaller starting points lists
    st_roc=[]
    for ii in range(len(gene_list)):
        st_new=[]
        if len(MSE_st2[ii])>200:
            for i in range(100):
                stb=MSE_st2[ii].sort_values(by='mean').iloc[i,6]
                st_new.append(stb)
        else:
            for i in range(int(len(MSE_st2[ii])/2)):
                stb=MSE_st2[ii].sort_values(by='mean').iloc[i,6]
                st_new.append(stb)
        st_roc.append(st_new)
    return [roc,st_roc]


def interactions_test(gene_list,network,consis,MSE_table,expression,log_expression,tolerance=[0.01,0.99]):
    """
    This function performs a stability test for different combinations of regulators (sources) 
    of target genes. It evaluates how removing specific combinations of regulators affects 
    the Mean Squared Error (MSE) of the model and fits the results to assess the impact of 
    these interactions on the target gene.

    Parameters:
    ----------
    gene_list : list
        List of target genes for which the stability test will be performed.
    
    network : pandas.DataFrame
        A dataframe representing the regulatory network. It contains interactions where each 
        row specifies a source gene, a target gene, and their interaction type.
    
    consis : list
        The result from the `consis_test` function. It contains the ROC-like consistency scores 
        and a list of starting points for each gene.
    
    MSE_table : list
        A list of data where each element corresponds to a gene and contains MSE values 
        and other related information, such as previously computed parameter sets. The 6th 
        element (index 5) should contain the optimized parameters for each gene (`res_pick_all`).
    
    expression : pandas.DataFrame
        A dataframe of scaled gene expression data where rows represent genes, and columns 
        represent samples or conditions.
    
    log_expression : pandas.DataFrame
        A dataframe of log-transformed gene expression data where rows represent genes, 
        and columns represent samples or conditions.

    tolerance : list, optional, default [0.01, 0.99]
        A list of two float values that define the tolerance range for filtering interactions 
        based on their ROC-like consistency score (actually the percentages of starting points
        can hit another interaction type after fitting).

    Returns:
    --------
    list :
        A list containing:
        - combs_all: A list of tested regulator combinations for each gene.
        - MSE_ratio_a: A list of MSE ratios for different regulator combinations, showing 
          how the removal of specific regulators affects the MSE.
        - res_fit_all: A list of optimized parameter sets for each gene after testing 
          different combinations of regulators.
        - MSE_value_all: A list of raw MSE values for each combination of regulators.
        - ob_new_all: A list of remaining regulators after removing specific combinations.

    """

    roc=consis[0]
    res_pick_all=MSE_table[5]
    combs_all=[]
    MSE_ratio_a=[]
    res_fit_all=[]
    ob_new_all = []
    MSE_value_all=[]
    for ii in range(len(gene_list)):
        choose=gene_list[ii]
        print(choose)
        ob = list(network[network['Target']==choose]['Source'])
        if choose in ob:
            ob.remove(choose)
        comb=get_combinations(search_list(roc[ii],ob,tolerance[0],tolerance[1]))
        combs_all.append(comb)
        MSE1=1
        MSE_ratio_dif=[]
        MSE_value=[]
        res_fit_part=[]
        for ii2 in range(len(res_pick_all[ii])):
            test1=get_args2(expression,ob=ob,target=gene_list[ii])
            arrays2=test1[:-1]
            y_pred1=fit_results(res_pick_all[ii][ii2], *arrays2)
            realdata=np.array(log_expression.loc[choose].iloc[0:],dtype='float64')
            MSE_test=mean_squared_error(realdata,y_pred1)
            MSE1=float(MSE_test)
            MSE_ratio=[]
            MSE_v=[]
            res_dif_fit=[]
            for dt in range(len(comb)):
                mask = np.isin(ob,comb[dt])
                idx = np.where(mask)[0]
                ob2 = [elem for i, elem in enumerate(ob) if i not in idx]
                ob_new_all.append(ob2)
                idx_par=[]
                for i in idx:
                    idx_par.append(3*i)
                    idx_par.append(3*i+1)
                    idx_par.append(3*i+2)
                parameters_new = [elem for i, elem in enumerate(res_pick_all[ii][ii2]) if i not in idx_par]
                if len(ob2)==0:
                    MSE_ratio.append(100000)
                    MSE_v.append(1)
                    res_dif_fit.append('Null')
                    break
                else:
                    restri_A = list(network[(network['Target'] == choose) & (network['Interaction'] == 2) &(network['Source'].isin(ob2)) ]['Source'])
                    restri_R = list(network[(network['Target'] == choose) & (network['Interaction'] == 0)&(network['Source'].isin(ob2))]['Source'])
                    positions_A = [ob2.index(item) for item in restri_A]
                    positions_R = [ob2.index(item) for item in restri_R]
                    nofeasB=boundary(len(ob2),act=positions_A,rep=positions_R)
                    tfexpression_use=expression.copy()
                    tfexpression_use.loc[choose]=log_expression.loc[choose]
                    test1=get_args(tfexpression_use,ob=ob2,target=choose,w=0)
                    resCl=minimize(fmodelOneHilltestw,parameters_new, args=test1, method='L-BFGS-B', jac='2-point', bounds=nofeasB)
                    res_dif_fit.append(resCl.x)
                    par2=resCl.x
                    test2=get_args2(expression,ob=ob2,target=gene_list[ii])
                    arrays3=test2[:-1]
                    y_pred2=fit_results(par2,*arrays3)
                    MSE2=mean_squared_error(realdata,y_pred2)
                    MSE_v.append(MSE2)
                    MSE_ratio.append(MSE2/MSE1)
            res_fit_part.append(res_dif_fit)
            MSE_ratio_dif.append(MSE_ratio)
            MSE_value.append(MSE_v)
        MSE_ratio_a.append(MSE_ratio_dif)
        MSE_value_all.append(MSE_value)
        res_fit_all.append(res_fit_part)
    return [combs_all,MSE_ratio_a,res_fit_all,MSE_value_all,ob_new_all]






def MSE_delete_fitting(gene_list,network,expression,log_expression,del_int,cc=5,num_tuning=20):
    """
    This function refines the MSE calculation by removing certain interactions (regulators) 
    from the regulatory network for each target gene, fitting new models, and recalculating 
    MSE values. It returns the optimized results after deleting the specified interactions.

    Parameters:
    ----------
    gene_list : list
        List of target genes for which MSE recalculations will be performed after deleting 
        certain regulators.
    
    network : pandas.DataFrame
        A dataframe representing the regulatory network. It contains interactions where each 
        row specifies a source gene, a target gene, and their interaction type.
    
    expression : pandas.DataFrame
        A dataframe of scaled gene expression data where rows represent genes, and columns 
        represent samples or conditions.
    
    log_expression : pandas.DataFrame
        A dataframe of log-transformed gene expression data where rows represent genes, 
        and columns represent samples or conditions.
    
    del_int : list
        A list where:
        - `del_int[0]` contains lists of regulators (sources) to be deleted for each target gene.
        - `del_int[1]` contains lists of starting points (`st_new_cut`) for each target gene, 
          used in the fitting process.
    
    cc : int, optional, default 5
        The number of cross-validation folds to use when tuning the model.
    
    num_tuning : int, optional, default 20
        The number of tuning iterations to perform when minimizing MSE.

    Returns:
    --------
    list :
        A list containing:
        - selected_total_sec: A list of selected starting points after deletion for each gene.
        - ob_newall: A list of remaining regulators after deletion for each gene.
        - res_par_sec: A list of optimized parameter values obtained by minimizing MSE 
          for each gene after the deletion of specified interactions.
        - MSE_min_total_sec: A list of minimized MSE values for each gene after tuning 
          the model with the updated regulators.
        - MSE_st2_sec: A list of updated dataframes containing MSE values for each gene, 
          along with the new starting points after the deletion of interactions.

    """
    obs=[]
    for ii in range(len(gene_list)):
        choose=gene_list[ii]

        ob = list(network[network['Target']==choose]['Source'])
        if choose in ob:
            ob.remove(choose)
        obs.append(ob)
    st_new_cut=del_int[1]
    delt=del_int[0]
    selected_total_sec=[]
    res_par_sec=[]
    MSE_min_total_sec=[]
    MSE_st2_sec=[]
    ob_newall=[]
    for ii in range(len(gene_list)):
        # work flow
        choose=gene_list[ii]
        mask = np.isin(obs[ii],delt[ii])
        idx = np.where(mask)[0]
        ob = [elem for i, elem in enumerate(obs[ii]) if i not in idx]
        ob_newall.append(ob)
        print(choose)
        tfexpression_use=expression.copy()
        tfexpression_use.loc[choose]=log_expression.loc[choose]
#### 1. choose the regulators gene   ### get a series of stPoints
        random.seed(1)
        restri_A = list(network[(network['Target'] == choose) & (network['Interaction'] == 2) &(network['Source'].isin(ob)) ]['Source'])
        restri_R = list(network[(network['Target'] == choose) & (network['Interaction'] == 0)&(network['Source'].isin(ob))]['Source'])
        positions_A = [ob.index(item) for item in restri_A]
        positions_R = [ob.index(item) for item in restri_R]
        nofeasB=boundary(len(ob),act=positions_A,rep=positions_R)
        test1=get_args2(tfexpression_use,ob=ob,target=gene_list[ii]) # get the initial array
        stdatabase=st_new_cut[ii]
    #stdatabase=stpoints_nums(MSE_st[ii],10)
        MSE_test_term= CV_fit(ob=ob,t_array=test1,stdatabase=stdatabase,boundary_cv=nofeasB)
        means=MSE_test_term.mean(axis=1)
        MSE_test_term['mean']=means
        MSE_test_term['mean'] = MSE_test_term['mean'].astype(float)
        MSE_st2_sec.append(MSE_test_term)
        stnumber=MSE_test_term.nsmallest(4, 'mean').index # get the row number of the stPoint database
    ## New stPoint dataset
        stnumber=stnumber.tolist()
        selected_stpoint = [stdatabase[i] for i in stnumber]
        selected_total_sec.append(selected_stpoint)
        test1=get_args(tfexpression_use,ob=ob,target=gene_list[ii],w=0)
        stPoint=selected_stpoint[0]
        resCl=minimize(fmodelOneHilltestw,stPoint, args=test1, method='L-BFGS-B', jac='2-point', bounds=nofeasB)
        res_par_sec.append(resCl.x)
        MSE_min=CV_tuning(ob=ob,t_array=test1,nums=num_tuning,stdatabase=selected_stpoint,boundary_cv=nofeasB,cc=cc)
        MSE_min_total_sec.append(MSE_min)
    for ct in range(len(gene_list)):
        MSE_st2_sec[ct]['SP'] = st_new_cut[ct]
    return [selected_total_sec,ob_newall,res_par_sec,MSE_min_total_sec,MSE_st2_sec]



def parameter_fitting(MSE_del_table,w_tot,gene_list,del_int,network,expression,log_expression):
    """
    This function performs parameter fitting after deleting specific interactions (regulators) 
    from the regulatory network for each target gene. It refines the model by fitting 
    parameters based on the updated network and computes the optimal parameters for each 
    gene after interaction deletion.

    Parameters:
    ----------
    MSE_del_table : list
        A list of data where each element corresponds to a gene and contains MSE values 
        and other related information, such as previously computed parameter sets. 
        The 5th element (index 4) should contain the MSE table (`MSE_st2_sec`).
    
    w_tot : list
        A list of weight paramaters (`w2`) for each gene, used for adjusting the contribution 
        of regulators to the model.

    gene_list : list
        List of target genes for which parameter fitting will be performed.

    del_int : list
        A list where:
        - `del_int[0]` contains lists of regulators (sources) to be deleted for each target gene.
        - `del_int[1]` contains lists of starting points for each gene, used in the fitting process.

    network : pandas.DataFrame
        A dataframe representing the regulatory network. It contains interactions where each 
        row specifies a source gene, a target gene, and their interaction type.

    expression : pandas.DataFrame
        A dataframe of scaled gene expression data where rows represent genes, and columns 
        represent samples or conditions.

    log_expression : pandas.DataFrame
        A dataframe of log-transformed gene expression data where rows represent genes, 
        and columns represent samples or conditions.

    Returns:
    --------
    list :
        A list containing:
        - res_final: A list of the final optimized parameter values for each gene after 
          interaction deletion.
        - ob: The list of remaining regulators (sources) for each gene after the deletion 
          of specific interactions.

    """
    delt=del_int[0]
    MSE_st2_sec=MSE_del_table[4]
    obs=[]
    for ii in range(len(gene_list)):
        choose=gene_list[ii]
        ob = list(network[network['Target']==choose]['Source'])
        if choose in ob:
            ob.remove(choose)
        obs.append(ob)
    res_w_all2=[]
    lam_w_all2=[]
    for ii in range(len(gene_list)):
        w2=w_tot[ii]
        choose=gene_list[ii]
        mask = np.isin(obs[ii],delt[ii])
        idx = np.where(mask)[0]
        ob = [elem for i, elem in enumerate(obs[ii]) if i not in idx]
        tfexpression_use=expression.copy()
        tfexpression_use.loc[choose]=log_expression.loc[choose]
        test1=get_args(tfexpression_use,ob=ob,target=gene_list[ii],w=w2)
        restri_A = list(network[(network['Target'] == choose) & (network['Interaction'] == 2) &(network['Source'].isin(ob)) ]['Source'])
        restri_R = list(network[(network['Target'] == choose) & (network['Interaction'] == 0)&(network['Source'].isin(ob))]['Source'])
        positions_A = [ob.index(item) for item in restri_A]
        positions_R = [ob.index(item) for item in restri_R]
        nofeasB=boundary(len(ob),act=positions_A,rep=positions_R)
        res_w=[]
        lam_w=[]
        for i in range(len(MSE_st2_sec[ii])):
            stPoint=MSE_st2_sec[ii].sort_values(by='mean').iloc[i,6]
            resCl=minimize(fmodelOneHilltestw,stPoint, args=test1, method='L-BFGS-B', jac='2-point', bounds=nofeasB)
            res_w.append(resCl.x)
            lam_w.append(lamlist(resCl.x))
        res_w_all2.append(res_w)
        lam_w_all2.append(lam_w)
    stbinss=[]
    for ii in range(len(gene_list)):
        stbins=[]
        for i in range(len(res_w_all2[ii])):
            stbins.append(checkstPoints(res_w_all2[ii][i]))
        stbinss.append(stbins)
    for ct in range(len(gene_list)):
        #MSE_st2_sec[ct]['SP'] = st_new_cut[ct]
        MSE_st2_sec[ct]['bin']=stbinss[ct]
    roc2=[]
    for ii in range(len(gene_list)):
        lenall=[0]*len(MSE_st2_sec[ii].sort_values(by='mean').iloc[0,7])
        if len(MSE_st2_sec[ii])>200:
            for i in range(100):
                stb=MSE_st2_sec[ii].sort_values(by='mean').iloc[i,7]
                lenall=sum_lists(lenall,stb,1)
            roc_new = [x/100 for x in lenall]
            roc2.append(roc_new)
        else:
            for i in range(int(len(MSE_st2_sec[ii])/2)):
                itt=len(MSE_st2_sec[ii])/2
                stb=MSE_st2_sec[ii].sort_values(by='mean').iloc[i,7]
                lenall=sum_lists(lenall,stb,1)
            roc_new = [x/itt for x in lenall]
            roc2.append(roc_new )
    stpoints_final=[]
    for i in range(len(gene_list)):
        st1=spresults_find(MSE_st2_sec[i],roc2[i])
        stpoints_final.append(st1)
    res_final=[]
    for ii in range(len(gene_list)):
        w2=w_tot[ii]
        choose=gene_list[ii]
        print(choose)
        mask = np.isin(obs[ii],delt[ii])
        idx = np.where(mask)[0]
        ob = [elem for i, elem in enumerate(obs[ii]) if i not in idx]
        print(ob)
        tfexpression_use=expression.copy()
        tfexpression_use.loc[choose]=log_expression.loc[choose]
        test1=get_args(tfexpression_use,ob=ob,target=gene_list[ii],w=w2)
        restri_A = list(network[(network['Target'] == choose) & (network['Interaction'] == 2) &(network['Source'].isin(ob)) ]['Source'])
        restri_R = list(network[(network['Target'] == choose) & (network['Interaction'] == 0)&(network['Source'].isin(ob))]['Source'])
        positions_A = [ob.index(item) for item in restri_A]
        positions_R = [ob.index(item) for item in restri_R]
        nofeasB=boundary(len(ob),act=positions_A,rep=positions_R)
        stPoint=stpoints_final[ii]
        resCl=minimize(fmodelOneHilltestw,stPoint, args=test1, method='L-BFGS-B', jac='2-point', bounds=nofeasB)
        res_final.append(resCl.x)
    return [res_final,ob]


def delete_int2(gene_list,network,int_results,consis_res,MSE_cut):
    MSE_value_all=int_results[3]
    MSE_ratio_a=int_results[1]
    combs_all=int_results[0]
    st_roc=consis_res[1]
    MSE_mean_s_all=[]
    for i2 in range(len(gene_list)):
        MSE_mean_s=[0]*len(MSE_ratio_a[i2][0])
        for i in range(len(MSE_ratio_a[i2])):
            MSE_mean_s=sum_lists(MSE_mean_s,MSE_ratio_a[i2][i],1)
        MSE_mean_s_all.append(MSE_mean_s)
    MSE_del=[]
    for i in range(len(gene_list)):
        data = {'Gene': combs_all[i], 'MSEtot': MSE_mean_s_all[i]}
        df = pd.DataFrame(data)
        MSE_del.append(df)
    MSE_mean_value_all=[]
    for i2 in range(len(gene_list)):
        MSE_mean_value=[0]*len(MSE_value_all[i2][0])
        for i in range(len(MSE_value_all[i2])):
            MSE_mean_value=sum_lists(MSE_mean_value,MSE_value_all[i2][i],1)
        MSE_mean_value_all.append(MSE_mean_value)
#create a df
    MSE_value_del=[]
    for i in range(len(gene_list)):
        data = {'Gene': combs_all[i], 'MSEtot': MSE_mean_s_all[i]}
        df = pd.DataFrame(data)
        MSE_value_del.append(df)
# finding the delete gene
    delt_value=[]
    for i in range(len(gene_list)):
        if len(MSE_value_del)>0:
            df=MSE_value_del[i][MSE_value_del[i]['MSEtot'] == MSE_cut].sort_values(by='MSEtot')
            if len(df)>0:
                max_len = max([len(row) for row in df['Gene']])
                rows_with_max_len = [row for row in df['Gene'] if len(row) == max_len]
                delt_value.append(rows_with_max_len[0])
            else:
                delt_value.append([])
        else:
            delt_value.append([])
    obs=[]
    for ii in range(len(gene_list)):
        choose=gene_list[ii]
        ob = list(network[network['Target']==choose]['Source'])
        if choose in ob:
            ob.remove(choose)
        obs.append(ob)
## get new st points
    st_new_cut=[]
    for ii in range(len(gene_list)):
        mask = np.isin(obs[ii],delt_value[ii])
        idx = np.where(mask)[0]
        ob2 = [elem for i, elem in enumerate(obs[ii]) if i not in idx]
        idx_par=[]
        for i in idx:
            idx_par.append(3*i)
            idx_par.append(3*i+1)
            idx_par.append(3*i+2)
        st_new=[]
        for ii2 in range(len(st_roc[ii])):
            parameters_new = [elem for i, elem in enumerate(st_roc[ii][ii2]) if i not in idx_par]
            st_new.append(parameters_new)
        st_new_cut.append(st_new)
    return [delt_value,st_new_cut,df,MSE_value_del]

def delete_int_mod(gene_list, network, int_results, consis_res, delt_value):
    MSE_value_all=int_results[3]
    MSE_ratio_a=int_results[1]
    combs_all=int_results[0]
    st_roc=consis_res[1]

    obs=[]
    for ii in range(len(gene_list)):
        choose=gene_list[ii]
        ob = list(network[network['Target']==choose]['Source'])
        if choose in ob:
            ob.remove(choose)
        obs.append(ob)
    ## get new st points
    ob_newall=[]
    st_new_cut=[]
    for ii in range(len(gene_list)):
        mask = np.isin(obs[ii],delt_value[ii])
        idx = np.where(mask)[0]
        ob2 = [elem for i, elem in enumerate(obs[ii]) if i not in idx]
        ob_newall.append(ob2)
        idx_par=[]
        for i in idx:
            idx_par.append(3*i)
            idx_par.append(3*i+1)
            idx_par.append(3*i+2)
        st_new=[]
        for ii2 in range(len(st_roc[ii])):
            parameters_new = [elem for i, elem in enumerate(st_roc[ii][ii2]) if i not in idx_par]
            st_new.append(parameters_new)
        st_new_cut.append(st_new)
    return [delt_value,st_new_cut,ob_newall]

def delete_int(gene_list,network,int_results,consis_res,MSE_cut):
    """
    This function identifies and deletes interactions (regulators) for each gene in the given list 
    based on MSE (Mean Squared Error) values. It calculates the average MSE for each interaction, 
    identifies the interaction to delete (if its MSE is below a certain threshold), and generates 
    new starting points for further modeling.

    Parameters:
    ----------
    gene_list : list
        List of target genes for which interaction deletion will be performed.

    network : pandas.DataFrame
        A dataframe representing the regulatory network. It contains interactions where each 
        row specifies a source gene, a target gene, and their interaction type.

    int_results : list
        A list of results from a previous interaction test, where:
        - `int_results[0]`: A list of regulator combinations (`combs_all`) for each gene.
        - `int_results[1]`: A list of MSE ratios (`MSE_ratio_a`) for each gene.
        - `int_results[3]`: A list of raw MSE values (`MSE_value_all`) for each gene.

    consis_res : list
        A list containing the ROC-like consistency scores and starting points (`st_roc`) 
        from a previous consistency test for each gene.

    MSE_cut : float
        The threshold MSE value used to determine which interactions (regulators) should be deleted. 
        Interactions with MSE below this value will be considered for deletion.

    Returns:
    --------
    list :
        A list containing:
        - delt_value: A list of regulators (interactions) to delete for each gene.
        - st_new_cut: A list of new starting points after deleting the identified interactions.
        - ob_newall: A list of remaining regulators (interactions) after deletion for each gene.
        - MSE_value_del: A list of dataframes showing MSE values for each regulator combination 
          for each gene, after interactions have been deleted.

    """
    MSE_value_all=int_results[3]
    MSE_ratio_a=int_results[1]
    combs_all=int_results[0]
    st_roc=consis_res[1]
    MSE_mean_s_all=[]
    for i2 in range(len(gene_list)):
        MSE_mean_s=[0]*len(MSE_ratio_a[i2][0])
        for i in range(len(MSE_ratio_a[i2])):
            MSE_mean_s=sum_lists(MSE_mean_s,MSE_ratio_a[i2][i],1)
        MSE_mean_s_all.append(MSE_mean_s)
    MSE_del=[]
    for i in range(len(gene_list)):
        data = {'Gene': combs_all[i], 'MSEtot': MSE_mean_s_all[i]}
        df = pd.DataFrame(data)
        MSE_del.append(df)
    MSE_mean_value_all=[]
    for i2 in range(len(gene_list)):
        MSE_mean_value=[0]*len(MSE_value_all[i2][0])
        for i in range(len(MSE_value_all[i2])):
            MSE_mean_value=sum_lists(MSE_mean_value,MSE_value_all[i2][i],1)
        MSE_mean_value_all.append(MSE_mean_value)
#create a df
    MSE_value_del=[]
    for i in range(len(gene_list)):
        data = {'Gene': combs_all[i], 'MSEtot': MSE_mean_s_all[i]}
        df = pd.DataFrame(data)
        MSE_value_del.append(df)


# finding the delete gene
    delt_value=[]
    for i in range(len(gene_list)):
        if len(MSE_value_del)>0:
            df=MSE_value_del[i][MSE_value_del[i]['MSEtot'] < MSE_cut].sort_values(by='MSEtot')
            if len(df)>0:
                max_len = max([len(row) for row in df['Gene']])
                rows_with_max_len = [row for row in df['Gene'] if len(row) == max_len]
                delt_value.append(rows_with_max_len[0])
            else:
                delt_value.append([])
        else:
            delt_value.append([])

    obs=[]
    for ii in range(len(gene_list)):
        choose=gene_list[ii]
        ob = list(network[network['Target']==choose]['Source'])
        if choose in ob:
            ob.remove(choose)
        obs.append(ob)
## get new st points
    ob_newall=[]
    st_new_cut=[]
    for ii in range(len(gene_list)):
        mask = np.isin(obs[ii],delt_value[ii])
        idx = np.where(mask)[0]
        ob2 = [elem for i, elem in enumerate(obs[ii]) if i not in idx]
        ob_newall.append(ob2)
        idx_par=[]
        for i in idx:
            idx_par.append(3*i)
            idx_par.append(3*i+1)
            idx_par.append(3*i+2)
        st_new=[]
        for ii2 in range(len(st_roc[ii])):
            parameters_new = [elem for i, elem in enumerate(st_roc[ii][ii2]) if i not in idx_par]
            st_new.append(parameters_new)
        st_new_cut.append(st_new)
    return [delt_value,st_new_cut,ob_newall,MSE_value_del]

##def score (cancer_row, naive_row, max_overlap, k):
##    tfs_cancer = cancer_row[0]
##    tfs_naive = naive_row[0]
##    mse = (cancer_row[1] + naive_row[1]) / 2
##    if max_overlap == 0:
##        percent = 1
##    else:
##        percent = len(list(set(tfs_cancer) & set(tfs_naive))) / max_overlap
##    if percent == 0:
##        score = 123456789
##    else:
##        score = mse / ((percent)**k)
##    return score
##
##def find_best(cancer_list, naive_list, cut_val, k):
##    df_tot = pd.DataFrame(index=range(len(cancer_list)),
##                          columns=range(len(naive_list)))
##    unique_naive = {item for lst in naive_list.iloc[:,0] for item in lst}
##    unique_cancer = {item for lst in cancer_list.iloc[:,0] for item in lst}
##    max_overlap = len(unique_naive & unique_cancer)
##    for i in range(len(cancer_list)):
##        for j in range(len(naive_list)):
##            val = score(cancer_list.iloc[i], naive_list.iloc[j], max_overlap, k)
##            if val > cut_val:
##                df_tot.iloc[i][j] = 9999
##            else:
##                df_tot.iloc[i][j] = val
##    df_lens = pd.DataFrame(index=range(len(cancer_list)),
##                          columns=range(len(naive_list)))
##    for i in range(len(cancer_list)):
##        for j in range(len(naive_list)):
##            if df_tot.iloc[i][j] == 9999:
##                df_lens.iloc[i][j] = 999
##            else:
##                df_lens.iloc[i][j] = 2 * max(len(cancer_list.iloc[i][0]),len(naive_list.iloc[j][0])) + min(len(cancer_list.iloc[i][0]),len(naive_list.iloc[j][0]))
##    len_min = df_lens.min().min()
##    for i in range(len(cancer_list)):
##        for j in range(len(naive_list)):
##            if df_lens.iloc[i][j] - 1 > len_min:
##                df_tot.iloc[i][j] = 9998
##    cancer_min, naive_min = df_tot.stack().idxmin()
##    cancer_keep = cancer_list.iloc[cancer_min]
##    naive_keep = naive_list.iloc[naive_min]
##    return cancer_keep, naive_keep, df_tot.min().min(),df_tot, df_lens
##
##def combined_del (del_naive, genes_naive, network_naive, del_cancer, genes_cancer, network_cancer, cv, k):
##    warnings.filterwarnings("ignore", category=FutureWarning)
##    df_naive = del_naive[3]
##    df_cancer = del_cancer[3]
##    genes_both = list(set(genes_naive) & set(genes_cancer))
##
##    naive_both = []
##    naive_order = []
##    for i in range(len(df_naive)):
##        if genes_naive[i] in genes_both:
##            naive_order.append(genes_naive[i])
##            df_temp = copy.deepcopy(df_naive[i])
##            total_genes = df_temp.iloc[-1][0]
##            df_temp.iloc[:, 0] = df_temp.iloc[:, 0].astype(object)
##            for j in range(len(df_temp)):
##                keeped = [gene for gene in total_genes if gene not in df_temp.iloc[j][0]]
##                df_temp.at[df_temp.index[j], df_temp.columns[0]] = keeped
##            naive_both.append(df_temp)
##
##    cancer_both = []
##    cancer_order = []
##    for i in range(len(df_cancer)):
##        if genes_cancer[i] in genes_both:
##            cancer_order.append(genes_cancer[i])
##            df_temp = copy.deepcopy(df_cancer[i])
##            total_genes = df_temp.iloc[-1][0]
##            df_temp.iloc[:, 0] = df_temp.iloc[:, 0].astype(object)
##            for j in range(len(df_temp)):
##                keeped = [gene for gene in total_genes if gene not in df_temp.iloc[j][0]]
##                df_temp.at[df_temp.index[j], df_temp.columns[0]] = keeped
##            cancer_both.append(df_temp)
##
##    naive_map1 = dict(zip(naive_order, naive_both))
##    cancer_map1 = dict(zip(cancer_order, cancer_both))
##    cut_value = 2 * cv
##    naive_both2 = []
##    for gene in naive_order:
##        df_temp = naive_map1[gene]
##        filtered = df_temp[df_temp.iloc[:, 1] < cut_value]
##        if len(filtered) > 0:
##            naive_both2.append(filtered)
##        else:
##            naive_both2.append(df_temp)
##    cancer_both2 = []
##    for gene in cancer_order:
##        df_temp = cancer_map1[gene]
##        filtered = df_temp[df_temp.iloc[:, 1] < cut_value]
##        if len(filtered) > 0:
##            cancer_both2.append(filtered)
##        else:
##            cancer_both2.append(df_temp)
##
##    naive_map2 = dict(zip(naive_order, naive_both2))
##    cancer_map2 = dict(zip(cancer_order, cancer_both2))
##
##    del_naive_update = copy.deepcopy(del_naive)
##    del_cancer_update = copy.deepcopy(del_cancer)
##
##    for gene in genes_both:
##        cancer_keep, naive_keep, min_MSE, df_tot, df_lens = find_best(cancer_map2[gene], naive_map2[gene], cut_value, k)
##        unique_naive = list({item for lst in naive_map1[gene].iloc[:, 0] for item in lst})
##        unique_cancer = list({item for lst in cancer_map1[gene].iloc[:, 0] for item in lst})
##        naive_tar = [a for a in unique_naive if a not in naive_keep[0]]
##        cancer_tar = [b for b in unique_cancer if b not in cancer_keep[0]]
##        if min_MSE < cut_value:
##            naive_gene_pos = genes_naive.index(gene)
##            cancer_gene_pos = genes_cancer.index(gene)
##            del_naive_update[0][naive_gene_pos] = naive_tar
##            del_cancer_update[0][cancer_gene_pos] = cancer_tar
##
##    return del_naive_update, del_cancer_update



def fmodelwithoutinput(x, *args):
    y = args[0]
    N = y.size

    # Sum up the models for all inputs and calculate the MSE
    mod_pred = x[0]  ## log2(x[-1]) to x[-1]
    resid = (mod_pred - y)**2
    func_val = sum(resid)
    mse = (1/N * func_val)

    return mse



def get_args_withoutinput(gene,target):# gene is a dataframe_tf12  ob is the name of the regulators
    argsss=[]
    resu2=np.array(gene.loc[target].iloc[ 0:],dtype='float64')
    argsss.append(resu2.flatten())
    argsss_tuple = tuple(argsss)
    return argsss_tuple





def Fit_withoutinput(gene_list, res_final, ob_newall,gene_position, tfexpression_logtarget):
    # Make copies of the input lists to avoid modifying the originals
    ob_newall2 = ob_newall.copy()
    res_final2 = res_final.copy()
    gene_position2=gene_position.copy()
    # Loop through the genes starting from the position of res_final length
    for gene in gene_list[len(res_final):]:
        print(gene)
        # Get the log expression data for the current gene
        logexpress = get_args_withoutinput(tfexpression_logtarget, gene)

        # Minimize the function and append results
        res = minimize(fmodelwithoutinput, [2], args=logexpress, method='L-BFGS-B', jac='2-point', bounds=[[-6, 6]])
        ob_newall2.append([])  # Append an empty list to ob_newall2
        res_final2.append(res.x)  # Append the result to res_final2
        gene_position2.append([])
    # Return the updated ob_newall2 and res_final2
    return ob_newall2, res_final2,gene_position2


def model_dynamic(pars, u):
    return (pars[0] + (1-pars[0])/(1+(u/pars[2])**(5*(pars[1])+1)))



def multiply_list(lst):
    result = 1
    for num in lst:
        result *= num
    return result



def inputs_func(n,dt,t_tot,tf_df,gene_list):
    input_n=n
    dt=dt
    t_tot=t_tot
    t_all = np.arange(0, t_tot+dt, dt)
    n_all=len(t_all)
    choose=gene_list[input_n]
    gene_s=np.array(tf_df.loc[choose].iloc[0:tf_df.shape[1]],dtype='float64')
    X_tot = np.array([gene_s[0]])  # initialize X_tot with first element of gene
    for ind in range(len(gene_s)-1):
        x_interp = np.linspace(gene_s[ind], gene_s[ind+1], int((n_all-1)/(len(gene_s) - 1)))
        X_tot = np.concatenate((X_tot, x_interp))
    X_tot2=np.insert(X_tot, 0, 0)
    X_tot1=np.append(X_tot,0)
    X_tot3=(X_tot1-X_tot2)
    input_tot=X_tot3[1:len(X_tot)+1]
    input_tot=np.insert(input_tot,0,0)
    input_tot=input_tot[0:len(input_tot)-1]
    return(input_tot)

###



def scale_array(arr):
    # Suppress all warnings during the execution of this block
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # Find the minimum and maximum values of the array
        min_val = np.min(arr)
        max_val = np.max(arr)

        # Scale the array between 0 and 1
        scaled_arr = (arr - min_val) / (max_val - min_val)

    return scaled_arr



def scale_array2(arr,arr2):
    # Find the minimum and maximum values of the array
    min_val = np.min(arr2)
    max_val = np.max(arr2)

    # Scale the array between 0 and 1
    scaled_arr = (arr - min_val) / (max_val - min_val)
    return scaled_arr


def cutnarray(array,n):
    l=int((len(array)-1)/n)
    cutresults=[]
    for i in range(n):
        cutresults.append(array[l*i:l*(i+1)+1])
    return(cutresults)


def dynamic1_1gene(dt,t_tot,exp_array,n,parameters_all,gene_position,save_num):

    dt_use=dt
    t_tot=t_tot
    t_all = np.arange(0, t_tot+dt, dt)
    n_all=len(t_all)
    colnum = np.linspace(1, n_all, save_num)
    #df = pd.DataFrame(index=range(len(exp_array)), columns=range(n_all))
    df = pd.DataFrame(index=range(len(exp_array)), columns=range(len(colnum)))
    n_value=exp_array[n]
    exp=exp_array
    for id in range(1,n_all+1):
        dxdt=dynamic_tot(last_step=exp,para_total=parameters_all,gene_position=gene_position,dt=dt_use)
        x=dxdt+exp
        x[n]=n_value
        exp=x.copy()
        if id in colnum:
            inde=np.where(colnum == id)[0][0]
            df.loc[:,inde]=x
    return(df)


def dynamic2_1gene(dt,t_tot,exp_array,n,parameters_all,gene_position,n_array,save_num):
    dt_use=dt
    t_tot=t_tot
    t_all = np.arange(0, t_tot+dt, dt)
    n_all=len(n_array)
    colnum = np.linspace(1, n_all, save_num)
    #df = pd.DataFrame(index=range(len(exp_array)), columns=range(n_all))
    df = pd.DataFrame(index=range(len(exp_array)), columns=range(len(colnum)))
    df.loc[:,0]=exp_array
    n_value=exp_array[n]
    exp=exp_array
    for id in range(1,n_all+1):
        dxdt=dynamic_tot(last_step=exp,para_total=parameters_all,gene_position=gene_position,dt=dt_use)
        x=dxdt+exp
        x[n]=exp[n]+n_array[id-1]
        exp=x.copy()
        if id in colnum:
            inde=np.where(colnum == id)[0][0]
            df.loc[:,inde]=x
    return(df)




def dynamic2(dt,t_tot,exp_array,n1,n2,parameters_all,gene_position,n1_array,n2_array,save_num):
    dt_use=dt
    t_tot=t_tot
    t_all = np.arange(0, t_tot+dt, dt)
    n_all=len(n1_array)
    colnum = np.linspace(0, n_all-1, save_num)
    df = pd.DataFrame(index=range(len(exp_array)), columns=range(len(colnum)))
    df.loc[:,0]=exp_array
    n1_value=exp_array[n1]
    n2_value=exp_array[n2]
    exp=exp_array
    for id in range(1,n_all+1):
        k1=np.array(dynamic_tot(last_step=exp,para_total=parameters_all,gene_position=gene_position,dt=dt_use))
        k2=np.array(dynamic_tot(last_step=exp+0.5*k1,para_total=parameters_all,gene_position=gene_position,dt=dt_use))
        k3=np.array(dynamic_tot(last_step=exp+0.5*k1+0.5*k2,para_total=parameters_all,gene_position=gene_position,dt=dt_use))
        k4=np.array(dynamic_tot(last_step=exp+k3,para_total=parameters_all,gene_position=gene_position,dt=dt_use))
        dxdt=(k1 + 2 * k2 + 2 * k3 + k4) / 6
        x=dxdt+exp
        x[n2]=exp[n2]+n2_array[id-1]
        x[n1]=exp[n1]+n1_array[id-1]
        exp=x.copy()
        if id in colnum:
            inde=np.where(colnum == id)[0][0]
            df.loc[:,inde]=x
    return(df)

def dynamic_singlegene(last_step, para, gene_position, gene_n, dt):
    """
    Calculate the dynamics of a single gene.
    
    Parameters:
    - last_step: array, values from the previous step
    - para: array, parameters used in the model
    - gene_position: index of the regulator genes in the gene list
    - gene_n: index of the target gene in the gene list
    - dt: time step size
    
    Returns:
    - dxdt: the change in gene expression for the current time step
    """
    gn = gene_position
    xi = last_step[gene_n]
    k = 0.1
    g = (2**(para[-1])) * k  # Growth rate scaling factor

    # If the length of para is 1, use the simplified equation
    if len(para) == 1:
        dxdt = (g - k * xi) * dt
    else:
        dxdt_parameter = []
        new_laststep = last_step[gn]

        # Compute the model dynamics for each regulator gene
        for ia, ib in enumerate(new_laststep):
            xm = para[ia*3:(ia+1)*3]
            dxdt_parameter.append(model_dynamic(xm, ib))

        # Combine the results to get the final dxdt
        dxdt = (multiply_list(dxdt_parameter) * g - k * xi) * dt

    return dxdt




def dynamic1(dt,t_tot,exp_array,n1,n2,parameters_all,gene_position,save_num):
    dt_use=dt
    t_tot=t_tot
    t_all = np.arange(0, t_tot+dt, dt)
    n_all=len(t_all)
    colnum = np.linspace(0, n_all-1, save_num)
    df = pd.DataFrame(index=range(len(exp_array)), columns=range(len(colnum)))
    df.loc[:,0]=exp_array
    n1_value=exp_array[n1]
    n2_value=exp_array[n2]
    exp=exp_array
    for id in range(1,n_all+1):
        k1=np.array(dynamic_tot(last_step=exp,para_total=parameters_all,gene_position=gene_position,dt=dt_use))
        k2=np.array(dynamic_tot(last_step=exp+0.5*k1,para_total=parameters_all,gene_position=gene_position,dt=dt_use))
        k3=np.array(dynamic_tot(last_step=exp+0.5*k1+0.5*k2,para_total=parameters_all,gene_position=gene_position,dt=dt_use))
        k4=np.array(dynamic_tot(last_step=exp+k3,para_total=parameters_all,gene_position=gene_position,dt=dt_use))
        dxdt=(k1 + 2 * k2 + 2 * k3 + k4) / 6
        x=dxdt+exp
        x[n2]=n2_value
        x[n1]=n1_value
        exp=x.copy()
        if id in colnum:
            inde=np.where(colnum == id)[0][0]
            df.loc[:,inde]=x
    return(df)


def dynamic_tot(last_step, para_total, gene_position, dt):
    dxdt = []
    nums = list(range(0, len(last_step)))  # Get range based on length of last_step
    nums2 = list(range(0, len(gene_position)))  # Get range based on length of gene_position

    for i in nums:
        if i in nums2:  # If index is valid for gene_position
            a = dynamic_singlegene(last_step=last_step, para=para_total[i], gene_position=gene_position[i], gene_n=i, dt=dt)
        else:  # If index is not in nums2, append last_step[i] directly
            a = 0
        dxdt.append(a)

    return dxdt


def selectcell(pseudo, tf, step=0.005, iterations=201, pseudotime_col=1, cellname_col='cellnames'):
    """
    This function selects specific cells based on their proximity to evenly spaced pseudotime 
    intervals and returns the corresponding transcription factor (TF) expression data for those cells.

    Parameters:
    ----------
    pseudo : pandas.DataFrame
        A dataframe that contains pseudotime information for cells. One column should hold the pseudotime 
        values, and another column should contain cell names.

    tf : pandas.DataFrame
        A dataframe where rows represent genes (e.g., transcription factors) and columns represent 
        individual cells. The columns should be named based on cell names, which correspond to the cell 
        names in `pseudo`.

    step : float, optional, default 0.005
        The step size used to divide the pseudotime into evenly spaced intervals.

    iterations : int, optional, default 201
        The number of evenly spaced pseudotime intervals to be used. The total range of pseudotime covered 
        will be `iterations * step`.

    pseudotime_col : int, optional, default 1
        The index (0-based) of the column in the `pseudo` dataframe that contains the pseudotime values.

    cellname_col : str, optional, default 'cellnames'
        The name of the column in the `pseudo` dataframe that contains the cell names.

    Returns:
    --------
    tuple :
        - pseudotime: A numpy array of the pseudotime values corresponding to the selected cells.
        - tfpickcell: A dataframe containing the transcription factor expression data for the selected cells.
        - cellpick: A list of the names of the selected cells.

    """

    nrows = []
    for i in range(iterations):
        closest_to_i = abs(pseudo.iloc[:, pseudotime_col] - i * step).idxmin()
        nrows.append(closest_to_i)

    cellpick = list(pseudo.iloc[nrows][cellname_col])
    tfpickcell = tf[cellpick]
    pseudotime = np.array(pseudo.iloc[:, pseudotime_col])

    return pseudotime, tfpickcell,cellpick


def process_tf(tf, tfpickcell, pseudo, cellpick, pseudotime_col=1):
    """
    This function processes transcription factor (TF) expression data for selected cells by applying 
    transformations, scaling the expression values, and returning the processed data along with 
    corresponding pseudotime values.

    Parameters:
    ----------
    tf : pandas.DataFrame
        A dataframe where rows represent genes (e.g., transcription factors) and columns represent 
        individual cells. The columns should be named based on cell names, which correspond to 
        the cell names in `pseudo`.

    tfpickcell : pandas.DataFrame
        A dataframe containing the transcription factor expression data for the selected cells. 
        This data is typically a subset of the `tf` dataframe, corresponding to the `cellpick` list.

    pseudo : pandas.DataFrame
        A dataframe that contains pseudotime information for cells. One column should hold the pseudotime 
        values, and another column should contain cell names.

    cellpick : list
        A list of the names of the selected cells for which the transcription factor expression data 
        is being processed.

    pseudotime_col : int, optional, default 1
        The index (0-based) of the column in the `pseudo` dataframe that contains the pseudotime values.

    Returns:
    --------
    tuple :
        - tfexpression: The transcription factor expression data after applying the power-of-two transformation.
        - tfexpression_logtarget: The scaled and log2-transformed transcription factor expression data.
        - pseudotime_pick: The pseudotime values corresponding to the selected cells.

    """
    # Step 1: Apply power of two transformation to tf12pickcell
    tfexpression = tfpickcell.apply(lambda x: x.map(power_of_two))

    # Step 2: Get pseudotime for picked cells
    col_names = tf.columns
    tot_cell = col_names[0:len(col_names)]
    pick_pseudotime = [tot_cell.get_loc(col_name) for col_name in cellpick]
    pseudotime = np.array(pseudo.iloc[:, pseudotime_col])
    pseudotime_pick = pseudotime[pick_pseudotime]

    # Step 3: Scale gene expression to max of one for each gene
    tfexpression_scaled = tfexpression.div(tfexpression.max(axis=1), axis=0)

    # Step 4: Apply log2 transformation to the scaled data
    tfexpression_logtarget = tfexpression_scaled.apply(lambda x: x.map(transtolog2))

    return tfexpression,tfexpression_logtarget, pseudotime_pick



def save_combined_network(gene_list, ob_newall, res_final, output):
    """
    This function builds a combined network of gene interactions based on the final fitted 
    interaction parameters for each gene in `gene_list` and saves it to a CSV file.

    Parameters:
    ----------
    gene_list : list
        A list of target genes for which the network interactions are being built.
    
    ob_newall : list
        A list where each element corresponds to the remaining regulators (sources) 
        for each gene in `gene_list` after some interactions have been deleted.

    res_final : list
        A list where each element contains the final fitted interaction parameters 
        for the regulators of the corresponding gene in `gene_list`. Each element is 
        a list of interaction parameters for each regulator.

    output : str
        The file path where the combined network will be saved as a CSV file.

    Returns:
    --------
    combined_network : pandas.DataFrame
        A dataframe representing the combined network of interactions for all genes in `gene_list`. 
        It contains the following columns:
        - 'Source': The regulator gene.
        - 'Target': The target gene.
        - 'Interaction': The interaction parameter between the regulator and the target gene.
        - 'n': The n parameter of the interaction.
        - 'thershold': The threshold parameter of the interaction.

    """
    combined_network = pd.DataFrame()

    # Loop through unique values to build the combined network DataFrame
    for k in range(len(gene_list)):
        choose = gene_list[k]
        ob = ob_newall[k]
        inter = res_final[k]
        inter2 = inter[:-1]  # Skip the last interaction
        list2 = influence(inter2)

        first = [inner_list[0] for inner_list in list2]
        second = [inner_list[1] for inner_list in list2]
        third = [inner_list[2] for inner_list in list2]

        target = [choose] * len(ob)
        df0 = pd.DataFrame(list(zip(ob, target, first, second, third)),
                           columns=['Source', 'Target', 'Interaction', 'n', 'thershold'])

        combined_network = pd.concat([combined_network, df0], ignore_index=True)

    combined_network.to_csv(output, index=False)

    return combined_network




def plot_sequence_vs_edges(gene_list, network, int_results, consis_res, start=10, stop=1001, step=1, figsize=(10, 6), plot_type='scatter'):
    """
    Plots the relationship between a sequence of MSE cut-off values and the number of edges in 
    a gene interaction network after applying a deletion function that removes interactions based on MSE_cut.

    Parameters:
    ----------
    gene_list : list
        List of unique genes for which the interactions are being evaluated.

    network : pandas.DataFrame
        The network dataframe that contains interactions between genes (source and target).

    int_results : pandas.DataFrame
        A dataframe containing interaction results. This should include information about 
        interaction strength and MSE for each gene interaction.

    consis_res : pandas.DataFrame
        A dataframe containing consistency results. This data is used to evaluate the stability 
        of interactions and is considered when applying the deletion function.

    start : float, optional, default 10
        The start value for the sequence of MSE cut-off values used in the deletion function.

    stop : float, optional, default 1001
        The stop value for the sequence of MSE cut-off values used in the deletion function.

    step : float, optional, default 1
        The step size for the sequence of MSE cut-off values used in the deletion function.

    figsize : tuple, optional, default (10, 6)
        The size of the plot (width, height).

    plot_type : str, optional, default 'scatter'
        The type of plot to generate: 'scatter' for a scatter plot, or 'line' for a line plot.

    Returns:
    --------
    pd.DataFrame :
        A dataframe containing the MSE cut-off sequence (as 'Sequence') and the corresponding number 
        of edges remaining in the network (as 'n_edges') after applying the deletion function.

    """
    # Generate the sequence
    sequence = np.arange(start, stop, step)
    lines_MSEcut = pd.DataFrame(sequence, columns=['Sequence'])
    tot_gene = []

    # Loop through sequence and apply delete_int function
    for cut_choose in sequence:
        obdele = delete_int(gene_list=gene_list, network=network, int_results=int_results, consis_res=consis_res, MSE_cut=cut_choose)[2]
        total_strings = sum(len(sublist) for sublist in obdele)
        tot_gene.append(total_strings)

    # Store the results in the DataFrame
    lines_MSEcut['n_edges'] = tot_gene

    # Plotting
    plt.figure(figsize=figsize)
    if plot_type == 'scatter':
        plt.scatter(lines_MSEcut['Sequence'], lines_MSEcut['n_edges'])
    elif plot_type == 'line':
        plt.plot(lines_MSEcut['Sequence'], lines_MSEcut['n_edges'])

    # Plot customization
    plt.title('Sequence vs n_edges')
    plt.xlabel('Sequence')
    plt.ylabel('n_edges')
    plt.grid(True)
    plt.show()

    return lines_MSEcut




def plot_network_results(gene_list, network, tfexpression, ob_newall, res_final,
                         tfexpression_logtarget, pseudotime_pick, output_file='Network_fit.pdf'):
    """
    This function generates and saves a series of subplots that compare model-predicted gene 
    expression to actual gene expression over pseudotime for a list of genes. The plots are 
    saved to a PDF file, with each subplot corresponding to one gene.

    Parameters:
    ----------
    gene_list : list
        A list of genes for which the model predictions and actual expression data will be plotted.

    network : pandas.DataFrame
        A dataframe representing the regulatory network. It contains interactions where each 
        row specifies a source gene, a target gene, and their interaction type.

    tfexpression : pandas.DataFrame
        A dataframe where rows represent transcription factors or genes, and columns represent 
        cell pseudotime values or conditions. It contains the expression data to be used for 
        model fitting and comparisons.

    ob_newall : list
        A list of remaining regulators (sources) for each gene after certain interactions 
        have been deleted. Each element corresponds to a gene in `gene_list`.

    res_final : list
        A list containing the final fitted interaction parameters for the remaining regulators 
        of each gene in `gene_list`.

    tfexpression_logtarget : pandas.DataFrame
        A dataframe similar to `tfexpression`, but where the expression values are log-transformed.

    pseudotime_pick : numpy.ndarray
        A numpy array of pseudotime values corresponding to the selected cells. This will be 
        used for plotting the x-axis.

    output_file : str, optional, default 'Network_fit.pdf'
        The file path where the PDF containing the subplots will be saved.

    Returns:
    --------
        The function generates subplots for each gene and saves them in a PDF file. 
        No data is returned.

    """
    from math import ceil

    # Open a PDF file to save the plots
    pdf_pages = PdfPages(output_file)

    # Calculate the number of rows needed (4 plots per row)
    num_plots = len(gene_list)
    num_cols = 4
    num_rows = ceil(num_plots / num_cols)

    # Adjust the figure size dynamically based on the number of rows
    fig_width = 20
    fig_height_per_row = 5  # Set the height per row
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(fig_width, fig_height_per_row * num_rows))

    # Ensure axs is always a 2D array even if there's only one row
    if num_rows == 1:
        axs = np.expand_dims(axs, axis=0)

    # Defining improved visual styles
    data_color = '#377eb8'  # Subdued, publication-friendly blue for data
    model_color = '#e41a1c'  # Vivid, contrast-rich red for model predictions
    line_width = 2

    # Iterate over each gene in the gene list to generate subplots
    for ii, gene in enumerate(gene_list):
        row = ii // num_cols  # Determine the row of the subplot
        col = ii % num_cols   # Determine the column of the subplot
        ax = axs[row, col]

        # Get the sources for the current target gene
        ob = list(network[network['Target'] == gene]['Source'])
        if gene in ob:
            ob.remove(gene)  # Remove self-regulation if it exists

        # Process expression data and fitting for the current gene
        test1 = get_args2(tfexpression, ob=ob, target=gene)
        arrays2 = test1[:-1]
        test2 = get_args2(tfexpression, ob=ob_newall[ii], target=gene)
        arrays3 = test2[:-1]
        y_pred2 = fit_results(res_final[ii], *arrays3)

        # Get the original expression data for the target gene
        orgin = np.array(tfexpression_logtarget.loc[[gene]].iloc[:, 0:tfexpression_logtarget.shape[1]], dtype='float64')

        # Plot the model predictions and actual data
        ax.plot(pseudotime_pick, 2**y_pred2, color=model_color, linewidth=line_width, label='Fitting')
        ax.plot(pseudotime_pick, 2**orgin[0], color=data_color, linewidth=line_width, linestyle='--', label='Data')

        # Customize the appearance of the plot
        ax.grid(False)  # Disable grid
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Set title and axis labels
        ax.set_title(gene, fontsize=40, fontweight='bold')
        ax.set_xlabel('Pseudotime', fontsize=30)
        ax.set_ylabel('Log-expression', fontsize=30)

        # Customize tick parameters
        ax.tick_params(axis='both', which='major', labelsize=20)

    # Turn off unused subplots
    for jj in range(ii + 1, num_rows * num_cols):
        fig.delaxes(axs[jj // num_cols, jj % num_cols])

    # Ensure no layout overlap and save the figure
    plt.tight_layout()
    pdf_pages.savefig(fig)
    pdf_pages.close()



def process_gene_expression(tfexpression_logtarget, tfexpression, gene_list, ob_newall):
    """
    Processes the gene expression data by filtering, scaling, and reindexing based on a list 
    of unique genes and their observed values.

    Parameters:
    ----------
    tfexpression_logtarget : pd.DataFrame
        Log-transformed gene expression data. Rows represent genes (and other associated indices), 
        and columns represent different cells or conditions.

    tfexpression : pd.DataFrame
        Target gene expression data before transformation. This is the raw or scaled expression 
        data that will be processed and reindexed.

    gene_list : pd.Index or list
        A list or index of unique gene names for which the gene expression data will be filtered 
        and reindexed.

    ob_newall : list of lists
        A list of observed gene values for each gene in `gene_list`. Each entry corresponds to 
        the list of observed regulators for the corresponding gene in `gene_list`.

    Returns:
    --------
    tuple :
        - tfexpression_logtarget2 : pd.DataFrame
            The log-transformed gene expression data, filtered and reordered such that genes 
            in `gene_list` appear first.
        
        - tfexpression_target2 : pd.DataFrame
            The scaled and reindexed target gene expression data, reindexed based on the reordered 
            `tfexpression_logtarget2`.

        - gene_position : list of lists
            A list of positions of the observed regulators in `ob_newall`, corresponding to the 
            indices of the reordered target gene list.

        - gene_list2 : np.ndarray
            An array of reordered gene names based on the processed `tfexpression_logtarget2` dataframe.

    """
    # Step 1: Filter the logtarget DataFrame using unique values
    tfexpression_target=tfexpression.copy()
    for i in range(tfexpression.shape[0]): # scale the max for gene expression is one
        arr_target=np.array(tfexpression_target.iloc[i])
        min_val = np.min(arr_target)
        scale_target = arr_target.copy()
        scale_target = scale_target/max(scale_target)
        tfexpression_target.iloc[i]=scale_target




    rows_in_gene_list = tfexpression_logtarget.index.get_level_values(0).isin(gene_list)
    # Get rows that are in the gene_list and reorder them
    tf_in_list = tfexpression_logtarget.loc[rows_in_gene_list].reindex(gene_list, level=0)
    # Get rows that are not in the gene_list
    tf_not_in_list = tfexpression_logtarget.loc[~rows_in_gene_list]
    # Concatenate the rows, with those in gene_list first and the rest last
    tfexpression_logtarget2 = pd.concat([tf_in_list, tf_not_in_list])
    gene_list2 = tfexpression_logtarget2.index[:len(gene_list)]
    tot_gene  = tfexpression_logtarget2.axes[0].tolist()
    gene_position = []
    index_values = []
    # Step 2: Reindex based on observed values in ob_newall
    for ii in range(len(gene_list)):
        choose = gene_list[ii]
        ob = ob_newall[ii]
        indices = [tot_gene.index(x) for x in ob]  # Get indices of observed genes
        gene_position.append(indices)
        index_values.append((choose,))  # Create tuples for reindexing

    # Step 3: Reindex the target gene expression data
    tfexpression_target2 = tfexpression_target.copy()
    tfexpression_target2 = tfexpression_target2.reindex(tfexpression_logtarget2.index.tolist())

    return tfexpression_logtarget2, tfexpression_target2, gene_position,gene_list2


def two_gene_driving(dt,t_tot,gene_list, drivers, tfexpression_target2, res_final, gene_position, tfexpression, ob_all, network, cached, cache):
    # Identify genes not in drivers and their positions
    #gene_list2=np.array([tup[0] for tup in tfexpression_target2.index.tolist()],dtype=object)
    #genes_not_in_drivers = np.setdiff1d(gene_list2, drivers)
    #positions = list(np.where(np.isin(gene_list2, genes_not_in_drivers))[0])
    #positions = [int(x) for x in positions]
    # Generate all combinations of genes
    all_combinations = []
    for i in range(len(gene_list)):
        for j in range(i + 1, len(gene_list)):
            all_combinations.append([i, j])
    gene_combinations = []
    for ts in all_combinations:
        gene_combinations.append(list(gene_list[ts]))

    # Initialize lists for storing MSE values
    MSE_twolines_all_f = []
    MSE_twolines_all_b = []
    MSE_twolines_all_dif = []

    # Iterate through all combinations and compute MSE
    for t_c, (n_1, n_2) in enumerate(all_combinations):
        genes = gene_combinations[t_c]
        print(genes)
        cache_key = tuple(sorted(genes))
        if (cached and cache_key in cache):
            df3, df4 = cache[cache_key]
        else:
            # Call the functions as per the logic given
            n1_array = inputs_func(n=n_1, dt=dt, t_tot=t_tot, tf_df=tfexpression_target2, gene_list=gene_list)
            n2_array = inputs_func(n=n_2, dt=dt, t_tot=t_tot, tf_df=tfexpression_target2, gene_list=gene_list)
            exp_array = np.array(tfexpression_target2.iloc[:, 0])

            df2 = dynamic1(dt=dt, t_tot=t_tot, exp_array=exp_array, n1=n_1, n2=n_2, parameters_all=res_final, gene_position=gene_position, save_num=201)
            exp_array = np.array(df2.iloc[:, -1])

            df3 = dynamic2(dt=dt, t_tot=t_tot, exp_array=exp_array, n1=n_1, n2=n_2, parameters_all=res_final, gene_position=gene_position,
                       n1_array=n1_array, n2_array=n2_array, save_num=201)

            n1_array_reverse = -np.flip(n1_array)
            n2_array_reverse = -np.flip(n2_array)
            exp_array = np.array(df3.iloc[:, -1])

            df4 = dynamic2(dt=dt, t_tot=t_tot, exp_array=exp_array, n1=n_1, n2=n_2, parameters_all=res_final, gene_position=gene_position,
                       n1_array=n1_array_reverse, n2_array=n2_array_reverse, save_num=201)

        MSE_twolines_f = []
        MSE_twolines_b = []
        MSE_twolines_dif = []

        for ii, choose in enumerate(gene_list):
            ob = list(network[network['Target'] == choose]['Source'])
            if choose in ob:
                ob.remove(choose)
            test2 = get_args2(tfexpression, ob=ob_all[ii], target=gene_list[ii])
            arrays3 = test2[:-1]
            y_pred2 = fit_results(res_final[ii], *arrays3)
            test4 = np.array(df3.iloc[ii], dtype='float64')
            test4_2 = scale_array(test4)
            test41 = np.array(df4.iloc[ii], dtype='float64')
            test41_2 = np.flip(scale_array(test41))

            if np.isnan(test4_2).any() or np.isnan(test41_2).any():
                ms, ms2, ms3 = 10000, 10000, 10000
            else:
                y_pred2_scale = scale_array(2 ** y_pred2)
                ms = mean_squared_error(y_pred2_scale, test4_2) if not np.isnan(y_pred2_scale).any() else 0
                ms2 = mean_squared_error(test41_2, test4_2) if not np.isnan(test41_2).any() else 0
                ms3 = mean_squared_error(y_pred2_scale, test41_2) if not np.isnan(y_pred2_scale).any() else 0

            MSE_twolines_f.append(ms)
            MSE_twolines_b.append(ms3)
            MSE_twolines_dif.append(ms2)

        MSE_twolines_all_f.append(MSE_twolines_f)
        MSE_twolines_all_b.append(MSE_twolines_b)
        MSE_twolines_all_dif.append(MSE_twolines_dif)
        print(f'{t_c / len(all_combinations):.2f}')

    # Compute sum of MSE and MSE^2
    sum_MSE_f = [sum(mse_f) for mse_f in MSE_twolines_all_f]
    sum_MSE_sqr_f = [sum(np.array(mse_f) ** 2) for mse_f in MSE_twolines_all_f]

    # Prepare DataFrame for MSE results
    sign_results_f = pd.DataFrame(index=range(len(all_combinations)), columns=range(6))
    sign_results_f.loc[:, 4] = sum_MSE_f
    sign_results_f.loc[:, 5] = sum_MSE_sqr_f

    for i in range(len(gene_combinations)):
        sign_results_f.loc[i, 0] = gene_combinations[i][0]
        sign_results_f.loc[i, 1] = gene_combinations[i][1]
        sign_results_f.loc[i, 2] = all_combinations[i][0]
        sign_results_f.loc[i, 3] = all_combinations[i][1]

    df_2lines_f = sign_results_f.sort_values(by=5)
    df_2lines_f.columns = ['gene1', 'gene2', 'position1', 'position2', 'MSE', 'MSE^2']

    #Repeat for MSE backwards results
    sum_MSE_b = [sum(mse_b) for mse_b in MSE_twolines_all_b]
    sum_MSE_sqr_b = [sum(np.array(mse_b) ** 2) for mse_b in MSE_twolines_all_b]

    # Prepare DataFrame for MSE results
    sign_results_b = pd.DataFrame(index=range(len(all_combinations)), columns=range(6))
    sign_results_b.loc[:, 4] = sum_MSE_b
    sign_results_b.loc[:, 5] = sum_MSE_sqr_b

    for i in range(len(gene_combinations)):
        sign_results_b.loc[i, 0] = gene_combinations[i][0]
        sign_results_b.loc[i, 1] = gene_combinations[i][1]
        sign_results_b.loc[i, 2] = all_combinations[i][0]
        sign_results_b.loc[i, 3] = all_combinations[i][1]

    df_2lines_b = sign_results_b.sort_values(by=5)
    df_2lines_b.columns = ['gene1', 'gene2', 'position1', 'position2', 'MSE', 'MSE^2']

    # Repeat for the MSE differences
    sum_MSE_dif = [sum(mse_dif) for mse_dif in MSE_twolines_all_dif]
    sum_MSE_sqr_dif = [sum(np.array(mse_dif) ** 2) for mse_dif in MSE_twolines_all_dif]

    sign_results_dif = pd.DataFrame(index=range(len(all_combinations)), columns=range(6))
    sign_results_dif.loc[:, 4] = sum_MSE_dif
    sign_results_dif.loc[:, 5] = sum_MSE_sqr_dif

    for i in range(len(gene_combinations)):
        sign_results_dif.loc[i, 0] = gene_combinations[i][0]
        sign_results_dif.loc[i, 1] = gene_combinations[i][1]
        sign_results_dif.loc[i, 2] = all_combinations[i][0]
        sign_results_dif.loc[i, 3] = all_combinations[i][1]

    df_2lines_dif = sign_results_dif.sort_values(by=5)
    df_2lines_dif.columns = ['gene1', 'gene2', 'position1', 'position2', 'MSE', 'MSE^2']

    return df_2lines_f, df_2lines_b, df_2lines_dif




def one_gene_driving(dt,t_tot,gene_list, drivers, tfexpression_target2, res_final, gene_position, network, tfexpression, ob_all, cached, cache):
    #genes_not_in_drivers = np.setdiff1d(gene_list, drivers)
    #positions = np.where(np.isin(gene_list, genes_not_in_drivers))[0]

    # Initialize result lists
    MSE_oneline_all_f = []
    MSE_oneline_all_b = []
    MSE_oneline_all_dif = []

    # Range of genes to loop over, excluding positions found in 'positions'x
    range_gene = [i for i in range(len(gene_list))]
    for t_c in range_gene:
        n_1 = t_c
        if (cached and (gene_list[n_1],) in cache):
            df3, df4 = cache[(gene_list[n_1],)]
        else:
            # Perform calculations for n1_array and other required arrays
            n1_array = inputs_func(n=n_1, dt=dt, t_tot=t_tot, tf_df=tfexpression_target2, gene_list=gene_list)
            n1_array_reverse = -np.flip(n1_array)
            exp_array = np.array(tfexpression_target2.iloc[:, 0])
            # Dynamic gene functions
            df2 = dynamic1_1gene(dt=dt, t_tot=t_tot, exp_array=exp_array, n=n_1, parameters_all=res_final, gene_position=gene_position, save_num=len(tfexpression.iloc[0]))
            exp_array = np.array(df2.iloc[:, -1])
            df3 = dynamic2_1gene(dt=dt, t_tot=t_tot, exp_array=exp_array, n=n_1, parameters_all=res_final, gene_position=gene_position, n_array=n1_array, save_num=len(tfexpression.iloc[0]))
            exp_array = np.array(df3.iloc[:, -1])
            df4 = dynamic2_1gene(dt=dt, t_tot=t_tot, exp_array=exp_array, n=n_1, parameters_all=res_final, gene_position=gene_position, n_array=n1_array_reverse, save_num=len(tfexpression.iloc[0]))

        # MSE Calculations
        MSE_oneline_f = []
        MSE_oneline_b = []
        MSE_oneline_dif = []

        for ii in range(len(gene_list)):
            choose = gene_list[ii]
            ob = list(network[network['Target'] == choose]['Source'])
            if choose in ob:
                ob.remove(choose)

            # Get args for prediction
            test2 = get_args2(tfexpression, ob=ob_all[ii], target=gene_list[ii])
            arrays3 = test2[:-1]
            y_pred2 = fit_results(res_final[ii], *arrays3)

            # Scale and calculate MSE
            test4 = np.array(df3.iloc[ii], dtype='float64')
            test4_2 = scale_array(test4)
            test41 = np.array(df4.iloc[ii], dtype='float64')
            test41_2 = np.flip(scale_array(test41))

            y_pred2_scale = scale_array(2**y_pred2)

            if np.isnan(test4_2).any() or np.isnan(test41_2).any():
                ms, ms2, ms3 = 10000, 10000, 10000
            else:
                y_pred2_scale = scale_array(2 ** y_pred2)
                ms = mean_squared_error(y_pred2_scale, test4_2) if not np.isnan(y_pred2_scale).any() else 0
                ms2 = mean_squared_error(test41_2, test4_2) if not np.isnan(test41_2).any() else 0
                ms3 = mean_squared_error(y_pred2_scale, test41_2) if not np.isnan(y_pred2_scale).any() else 0

            MSE_oneline_f.append(ms)
            MSE_oneline_b.append(ms3)
            MSE_oneline_dif.append(ms2)



        MSE_oneline_all_f.append(MSE_oneline_f)
        MSE_oneline_all_b.append(MSE_oneline_b)
        MSE_oneline_all_dif.append(MSE_oneline_dif)

    # Summing the MSE results
    sum_MSE_f = [sum(mse_f) for mse_f in MSE_oneline_all_f]
    sum_MSE_sqr2_f = [sum(np.array(mse_f)**2) for mse_f in MSE_oneline_all_f]

    sum_MSE_b = [sum(mse_b) for mse_b in MSE_oneline_all_b]
    sum_MSE_sqr2_b = [sum(np.array(mse_b)**2) for mse_b in MSE_oneline_all_b]

    sum_MSE_dif = [sum(mse_dif) for mse_dif in MSE_oneline_all_dif]
    sum_MSE_sqr2_dif = [sum(np.array(mse_dif)**2) for mse_dif in MSE_oneline_all_dif]

    # Create DataFrames
    sign1_results_f = pd.DataFrame({'gene': gene_list[range_gene],"MSE":sum_MSE_f ,'MSE^2': sum_MSE_sqr2_f}).sort_values(by='MSE^2')
    sign1_results_b = pd.DataFrame({'gene': gene_list[range_gene],"MSE":sum_MSE_b ,'MSE^2': sum_MSE_sqr2_b}).sort_values(by='MSE^2')
    sign1_results_dif = pd.DataFrame({'gene': gene_list[range_gene], "MSE":sum_MSE_dif ,'MSE^2': sum_MSE_sqr2_dif}).sort_values(by='MSE^2')


    return sign1_results_f, sign1_results_b, sign1_results_dif



def driving_results(dri_genes, gene_list, tfexpression_target2, network, tfexpression, pseudotime_pick, res_final, gene_position, ob_newall, cached=False,cache=None):
    """
    Simulates dynamic gene expression based on driving genes and plots the forward and backward results.
    The function can handle one or two driving genes and caches the results for future use to improve efficiency.

    Parameters:
    ----------
    dri_genes : list
        List of driving genes used to control the dynamic simulations. Can contain one or two genes.

    gene_list : list
        List of target genes for which the simulation will be performed and plotted.

    tfexpression_target2 : pandas.DataFrame
        Scaled target gene expression data, used to calculate the driving signals.

    network : pandas.DataFrame
        A dataframe representing the regulatory network. It contains interactions where each 
        row specifies a source gene, a target gene, and their interaction type.

    tfexpression : pandas.DataFrame
        Target gene expression data before transformation, used to generate model predictions.

    pseudotime_pick : numpy.ndarray
        A numpy array of pseudotime values corresponding to the selected cells. Used for plotting the x-axis.

    res_final : list
        A list containing the final fitted interaction parameters for the remaining regulators 
        of each gene in `gene_list`.

    gene_position : list
        A list of positions of the observed regulators in `ob_newall`, corresponding to the indices 
        of the reordered target gene list.

    ob_newall : list of lists
        List of observed regulators for each gene in `gene_list`, used to generate the model predictions.

    cached : bool, optional, default False
        Whether to use cached results for faster future simulations.

    cache : dict, optional
        A dictionary used to store cached results based on the driving genes.

    Returns:
    --------
    tuple :
        - df3: The forward simulation results for each gene.
        - df4: The backward simulation results for each gene.
        - cache: Updated cache containing the results for the driving genes.
        - fig: The figure object containing the plots for forward and backward simulations.

    """
    from math import ceil
    if cache is None:
        cache = {}
    cache_key = tuple(dri_genes)

    # Check if the results are already in the cache
    if cached and cache_key in cache:
        df3, df4 = cache[cache_key]
        print("Using cached results.")
    else:
        if len(dri_genes) == 2:
            # Two gene logic
            positions = np.where(np.isin(gene_list, dri_genes))[0]
            n_1 = positions[0]
            n_2 = positions[1]

            # Calculate arrays for both driving genes
            n1_array = inputs_func(n=n_1, dt=1, t_tot=50000, tf_df=tfexpression_target2, gene_list=gene_list)
            n2_array = inputs_func(n=n_2, dt=1, t_tot=50000, tf_df=tfexpression_target2, gene_list=gene_list)

            # Initial expression array
            exp_array = np.array(tfexpression_target2.iloc[:, 0])

            # Perform dynamic simulations for two genes
            df2 = dynamic1(dt=0.1, t_tot=200, exp_array=exp_array, n1=n_1, n2=n_2, parameters_all=res_final, gene_position=gene_position, save_num=len(tfexpression.iloc[0]))
            exp_array = np.array(df2.iloc[:, -1])

            df3 = dynamic2(dt=1, t_tot=50000, exp_array=exp_array, n1=n_1, n2=n_2, parameters_all=res_final, gene_position=gene_position,
                           n1_array=n1_array, n2_array=n2_array, save_num=len(tfexpression.iloc[0]))

            # Reverse arrays for backward simulation
            n1_array_reverse = -np.flip(n1_array)
            n2_array_reverse = -np.flip(n2_array)
            exp_array = np.array(df3.iloc[:, -1])

            df4 = dynamic2(dt=1, t_tot=50000, exp_array=exp_array, n1=n_1, n2=n_2, parameters_all=res_final, gene_position=gene_position,
                           n1_array=n1_array_reverse, n2_array=n2_array_reverse, save_num=len(tfexpression.iloc[0]))

        elif len(dri_genes) == 1:
            # One gene logic
            n_1 = np.where(np.isin(gene_list, dri_genes))[0][0]
            print(f"Processing single gene: {gene_list[n_1]}")

            # Calculate array for the driving gene
            n1_array = inputs_func(n=n_1, dt=1, t_tot=50000, tf_df=tfexpression_target2, gene_list=gene_list)
            n1_array_reverse = -np.flip(n1_array)

            # Initial expression array
            exp_array = np.array(tfexpression_target2.iloc[:, 0])

            # Perform dynamic simulations for one gene
            df2 = dynamic1_1gene(dt=0.01, t_tot=200, exp_array=exp_array, n=n_1, parameters_all=res_final, gene_position=gene_position, save_num=len(tfexpression.iloc[0]))
            exp_array = np.array(df2.iloc[:, -1])

            df3 = dynamic2_1gene(dt=1, t_tot=50000, exp_array=exp_array, n=n_1, parameters_all=res_final, gene_position=gene_position,
                                 n_array=n1_array, save_num=len(tfexpression.iloc[0]))

            exp_array = np.array(df3.iloc[:, -1])

            df4 = dynamic2_1gene(dt=1, t_tot=50000, exp_array=exp_array, n=n_1, parameters_all=res_final, gene_position=gene_position,
                                 n_array=n1_array_reverse, save_num=len(tfexpression.iloc[0]))

        # Cache results for faster future access
        cache[cache_key] = (df3, df4)

    # Dynamically calculate the number of rows needed for the plots
    num_plots = len(gene_list)
    num_cols = 4  # Fixed number of columns
    num_rows = ceil(num_plots / num_cols)  # Calculate required number of rows

    # Adjust the figure size dynamically based on the number of rows
    fig_width = 20
    fig_height_per_row = 5  # Set height per row
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(fig_width, fig_height_per_row * num_rows))

    # Flatten the axs array to handle 1D indexing
    axs = axs.flatten()

    plt.rcParams.update({'font.size': 16, 'font.family': 'sans-serif'})
    colors = ['#d95f02', '#1b9e77', '#7570b3']
    linewidths = [4, 4, 3]

    for ii in range(num_plots):
        ax = axs[ii]
        choose = gene_list[ii]
        ob = list(network[network['Target'] == choose]['Source'])
        if choose in ob:
            ob.remove(choose)

        # Perform the argument fetching and result fitting
        test2 = get_args2(tfexpression, ob=ob_newall[ii], target=gene_list[ii])
        arrays3 = test2[:-1]
        y_pred2 = fit_results(res_final[ii], *arrays3)
        orgin = np.array(tfexpression_target2.loc[[choose]].iloc[:, 0:tfexpression_target2.shape[1]], dtype='float64')

        # Forward simulation
        test01 = np.array(df3.iloc[ii], dtype='float64')
        test02 = scale_array(test01)
        arr1 = np.linspace(0, 1, len(test01))

        # Backward simulation
        test11 = np.array(df4.iloc[ii], dtype='float64')
        test12 = scale_array2(test11, test01)
        arr2 = np.linspace(0, 1, len(test11))

        # Plot results
        ax.plot(arr2[1:], test02[1:], color=colors[0], linewidth=linewidths[0], alpha=0.7, label='Forward signal result')
        ax.plot(arr1[1:], np.flip(test12[1:]), color=colors[1], linewidth=linewidths[1], alpha=0.7, label='Backward signal result')
        ax.plot(pseudotime_pick, scale_array(2**orgin[0]), color=colors[2], linewidth=linewidths[2], linestyle='--', label='Data')

        # Set title and labels
        ax.set_title(choose, fontsize=30, fontweight='bold')
        ax.set_xlabel('Pseudotime', fontsize=20)
        ax.set_ylabel('Log-expression', fontsize=20)
        ax.grid(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Add legend to the first subplot only
        if ii == 0:
            ax.legend(fontsize=15, loc='upper right')

    # Turn off unused subplots
    for jj in range(num_plots, len(axs)):
        fig.delaxes(axs[jj])

    # Adjust layout
    plt.tight_layout()
    plt.show()

    return df3, df4,cache,fig



import os
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy.stats import gaussian_kde
from scipy.stats import trim_mean


# Combined / overlap-aware edge deletion across two networks


def get_top_sign(genes, gene, network, int_test, kept, base_models):
    if len(kept) == 0:
        return {}
    ii = genes.index(gene)
    ob = list(network[network["Target"] == gene]["Source"])
    if gene in ob:
        ob.remove(gene)
    ob2 = [r for r in ob if r in set(kept)]
    deleted = [r for r in ob if r not in set(kept)]

    # case: did not delete anything
    if len(deleted) == 0:
        par = base_models[ii][0]  # index 0, best model
        t = (len(par) - 1) // 3
        lams = par[:3*t:3]
        signs = [1 if int(lam) > 1 else 0 for lam in lams]
        return dict(zip(ob2, signs))

    # find which deletion combination matches
    combs_all = int_test[0][ii]
    dt_idx = None
    for dt in range(len(combs_all)):
        if sorted(combs_all[dt]) == sorted(deleted):
            dt_idx = dt
            break
    res_fit_all = int_test[2][ii]
    par = res_fit_all[0][dt_idx]  # index 0, best model
    t = (len(par) - 1) // 3
    lams = par[:3*t:3]
    signs = [1 if int(lam) > 1 else 0 for lam in lams]
    return dict(zip(ob2, signs))


def overlap_score(mse1, mse2, regulators1, regulators2, pos1, pos2, n, k):
    if len(regulators1) == 0 or len(regulators2) == 0:
        return 0
    regulators1 = set(regulators1)
    regulators2 = set(regulators2)
    jaccard = len(regulators1 & regulators2) / len(regulators1 | regulators2)
    return 1/(mse1 + mse2) * (jaccard**k) * math.e**(-n * (len(regulators1) + len(regulators2)))


def overlap_score_signed(mse1, mse2, regulators1, regulators2, pos1, pos2, n, k):
    if len(regulators1) == 0 or len(regulators2) == 0:
        return 0
    regulators1 = set(regulators1)
    regulators2 = set(regulators2)
    if len(regulators1 | regulators2) == 0:
        jaccard = 0
    else:
        total = 0
        for r in (regulators1 & regulators2):
            p1 = pos1.get(r)
            p2 = pos2.get(r)
            if (p1 > 0.5 and p2 > 0.5) or (p1 < 0.5 and p2 < 0.5):
                total += 1
        jaccard = total / len(regulators1 | regulators2)
    return 1/(mse1 + mse2) * (jaccard**k) * math.e**(-n * (len(regulators1) + len(regulators2)))


def delete_overlapping(delete_1, delete_2, genes_1, genes_2, network_1, network_2,
                       int_test_1, int_test_2, Part3_1, Part3_2, n, k):
    ob_newall1 = [[] for i in range(len(genes_1))]
    ob_newall2 = [[] for i in range(len(genes_2))]
    del_dict_1 = dict(zip(genes_1, delete_1[3]))
    del_dict_2 = dict(zip(genes_2, delete_2[3]))
    genes_overlap = list(set(genes_1) & set(genes_2))
    for gene in genes_overlap:
        del_1 = del_dict_1[gene]
        del_2 = del_dict_2[gene]
        ob1 = list(network_1[network_1['Target']==gene]['Source'])
        ob2 = list(network_2[network_2['Target']==gene]['Source'])
        empty_row = pd.DataFrame({"Gene": [[]], "MSEtot": [10]})
        del_1 = pd.concat([empty_row, del_1], ignore_index=True)
        del_2 = pd.concat([empty_row, del_2], ignore_index=True)
        pos_cache_1 = {}
        pos_cache_2 = {}
        cur_best = -1
        cur_best_1 = []
        cur_best_2 = []
        for i in range(len(del_1)):
            for j in range(len(del_2)):
                del_1_temp = del_1.iloc[i]
                del_2_temp = del_2.iloc[j]
                regulators1 = [w for w in ob1 if w not in del_1_temp["Gene"]]
                key1 = tuple(sorted(regulators1))
                if key1 not in pos_cache_1:
                    pos_cache_1[key1] = get_top_sign(genes_1, gene, network_1, int_test_1, regulators1, Part3_1[5])
                pos1 = pos_cache_1[key1]
                regulators2 = [v for v in ob2 if v not in del_2_temp["Gene"]]
                key2 = tuple(sorted(regulators2))
                if key2 not in pos_cache_2:
                    pos_cache_2[key2] = get_top_sign(genes_2, gene, network_2, int_test_2, regulators2, Part3_2[5])
                pos2 = pos_cache_2[key2]
                score = overlap_score(del_1_temp["MSEtot"], del_2_temp["MSEtot"], regulators1, regulators2, pos1, pos2, n, k)
                if score > cur_best:
                    cur_best = score
                    cur_best_1 = regulators1
                    cur_best_2 = regulators2
        ob_newall1[genes_1.index(gene)] = cur_best_1
        ob_newall2[genes_2.index(gene)] = cur_best_2
    return ob_newall1, ob_newall2


def delete_overlapping_signed(delete_1, delete_2, genes_1, genes_2, network_1, network_2,
                               int_test_1, int_test_2, Part3_1, Part3_2, n, k):
    ob_newall1 = [[] for i in range(len(genes_1))]
    ob_newall2 = [[] for i in range(len(genes_2))]
    del_dict_1 = dict(zip(genes_1, delete_1[3]))
    del_dict_2 = dict(zip(genes_2, delete_2[3]))
    genes_overlap = list(set(genes_1) & set(genes_2))
    for gene in genes_overlap:
        del_1 = del_dict_1[gene]
        del_2 = del_dict_2[gene]
        ob1 = list(network_1[network_1['Target']==gene]['Source'])
        ob2 = list(network_2[network_2['Target']==gene]['Source'])
        empty_row = pd.DataFrame({"Gene": [[]], "MSEtot": [10]})
        del_1 = pd.concat([empty_row, del_1], ignore_index=True)
        del_2 = pd.concat([empty_row, del_2], ignore_index=True)
        pos_cache_1 = {}
        pos_cache_2 = {}
        cur_best = -1
        cur_best_1 = []
        cur_best_2 = []
        for i in range(len(del_1)):
            for j in range(len(del_2)):
                del_1_temp = del_1.iloc[i]
                del_2_temp = del_2.iloc[j]
                regulators1 = [w for w in ob1 if w not in del_1_temp["Gene"]]
                key1 = tuple(sorted(regulators1))
                if key1 not in pos_cache_1:
                    pos_cache_1[key1] = get_top_sign(genes_1, gene, network_1, int_test_1, regulators1, Part3_1[5])
                pos1 = pos_cache_1[key1]
                regulators2 = [v for v in ob2 if v not in del_2_temp["Gene"]]
                key2 = tuple(sorted(regulators2))
                if key2 not in pos_cache_2:
                    pos_cache_2[key2] = get_top_sign(genes_2, gene, network_2, int_test_2, regulators2, Part3_2[5])
                pos2 = pos_cache_2[key2]
                score = overlap_score_signed(del_1_temp["MSEtot"], del_2_temp["MSEtot"], regulators1, regulators2, pos1, pos2, n, k)
                if score > cur_best:
                    cur_best = score
                    cur_best_1 = regulators1
                    cur_best_2 = regulators2
        ob_newall1[genes_1.index(gene)] = cur_best_1
        ob_newall2[genes_2.index(gene)] = cur_best_2
    return ob_newall1, ob_newall2


def delete_non_overlapping(delete, genes, network, ob_newall, n):
    del_dict = dict(zip(genes, delete[3]))
    for i, gene in enumerate(genes):
        if ob_newall[i] != []:
            continue
        del_gene = del_dict[gene]
        ob = list(network[network['Target']==gene]['Source'])
        empty_row = pd.DataFrame({"Gene": [[]], "MSEtot": [10]})
        del_gene = pd.concat([empty_row, del_gene], ignore_index=True)
        cur_best = -1
        cur_best_genes = []
        for j in range(len(del_gene)):
            delete_temp = del_gene.iloc[j]
            regulators = [m for m in ob if m not in delete_temp["Gene"]]
            if len(regulators) == 0:
                score = 0
            else:
                score = (1/delete_temp["MSEtot"]) * math.e**(-2 * n * (len(regulators)))
            if score > cur_best:
                cur_best = score
                cur_best_genes = regulators
        ob_newall[i] = cur_best_genes
    return ob_newall


def make_stpoints(gene_list, network, ob_newall, consis_res):
    st_roc = consis_res[1]
    obs = []
    for g in gene_list:
        ob = list(network[network["Target"] == g]["Source"])
        if g in ob:
            ob.remove(g)
        obs.append(ob)
    delt_value = []
    for i in range(len(gene_list)):
        kept = set(ob_newall[i])
        delt_value.append([r for r in obs[i] if r not in kept])
    st_new_cut = []
    for i in range(len(gene_list)):
        mask = np.isin(obs[i], delt_value[i])
        idx = np.where(mask)[0]
        idx_par = []
        for j in idx:
            idx_par.extend([3*j, 3*j + 1, 3*j + 2])
        idx_par = set(idx_par)
        st_new = []
        for vec in st_roc[i]:
            vec = list(vec)
            st_new.append([val for p, val in enumerate(vec) if p not in idx_par])
        st_new_cut.append(st_new)
    return delt_value, st_new_cut


def run_delete_combined(network_1, genes_1, delete_1, consis_res_1, int_test_1, Part3_1,
                       network_2, genes_2, delete_2, consis_res_2, int_test_2, Part3_2,
                       n_val, k_val):
    ob_newall1, ob_newall2 = delete_overlapping(delete_1, delete_2, genes_1, genes_2, network_1, network_2,
                                                int_test_1, int_test_2, Part3_1, Part3_2, n_val, k_val)
    ob_newall1 = delete_non_overlapping(delete_1, genes_1, network_1, ob_newall1, n_val)
    ob_newall2 = delete_non_overlapping(delete_2, genes_2, network_2, ob_newall2, n_val)

    delt_value1, st_new_cut1 = make_stpoints(genes_1, network_1, ob_newall1, consis_res_1)
    delt_value2, st_new_cut2 = make_stpoints(genes_2, network_2, ob_newall2, consis_res_2)

    return ob_newall1, ob_newall2, delt_value1, delt_value2, st_new_cut1, st_new_cut2


def run_delete_combined_signed(network_1, genes_1, delete_1, consis_res_1, int_test_1, Part3_1,
                       network_2, genes_2, delete_2, consis_res_2, int_test_2, Part3_2,
                       n_val, k_val):
    ob_newall1, ob_newall2 = delete_overlapping_signed(delete_1, delete_2, genes_1, genes_2, network_1, network_2,
                                                int_test_1, int_test_2, Part3_1, Part3_2, n_val, k_val)
    ob_newall1 = delete_non_overlapping(delete_1, genes_1, network_1, ob_newall1, n_val)
    ob_newall2 = delete_non_overlapping(delete_2, genes_2, network_2, ob_newall2, n_val)

    delt_value1, st_new_cut1 = make_stpoints(genes_1, network_1, ob_newall1, consis_res_1)
    delt_value2, st_new_cut2 = make_stpoints(genes_2, network_2, ob_newall2, consis_res_2)

    return ob_newall1, ob_newall2, delt_value1, delt_value2, st_new_cut1, st_new_cut2


def ob_to_network(gene_list, ob_newall):

    rows = []
    for target, regs in zip(gene_list, ob_newall):
        if regs is None:
            continue
        for source in regs:
            rows.append({"Source": source, "Target": target, "Interaction": 1})

    df = pd.DataFrame(rows, columns=["Source", "Target", "Interaction"])
    return df


def calculate_network_mse(gene_list, network, tfexpression, ob_newall, res_final,
                          tfexpression_logtarget, trim=0.1):
    mse_dict = {}
    for ii, gene in enumerate(gene_list):
        ob = list(network[network['Target'] == gene]['Source'])
        test2 = get_args2(tfexpression, ob=ob_newall[ii], target=gene)
        arrays3 = test2[:-1]
        y_pred2 = fit_results(res_final[ii], *arrays3)
        orgin = np.array(tfexpression_logtarget.loc[[gene]].iloc[:, 0:tfexpression_logtarget.shape[1]], dtype='float64')
        mse = np.mean((2**y_pred2 - 2**orgin[0]) ** 2)
        mse_dict[gene] = float(mse)

    trimmed_mse = trim_mean(list(mse_dict.values()), trim)
    return mse_dict, trimmed_mse


# Stochastic force-and-release perturbation simulation


def dynamic_state(dt,t_tot,exp_array,parameters_all,gene_position,save_num,e_mean):
    dt_use=dt
    t_tot=t_tot
    t_all = np.arange(0, t_tot+dt, dt)
    n_all=len(t_all)
    colnum = np.linspace(1, n_all, save_num)
    #df = pd.DataFrame(index=range(len(exp_array)), columns=range(n_all))
    df = pd.DataFrame(index=range(len(exp_array)), columns=range(len(colnum)))
    exp=exp_array
    sde_n=len(exp)
    for id in range(1,n_all+1):
        dxdt=dynamic_tot(last_step=exp,para_total=parameters_all,gene_position=gene_position,dt=dt_use)
        #noise at 0.04
        x=dxdt+exp+0.04*e_mean*np.random.normal(0, dt_use, sde_n)
        x[x < 0] = 0
        exp=x.copy()
        if id in colnum:
            inde=np.where(colnum == id)[0][0]
            df.loc[:,inde]=x
    return(df)


def make_jump(gene_i, start_end, per, dt=1, t_tot=20):
    arr = np.zeros(int(t_tot / dt) + 1)
    arr[0] = (start_end[1][gene_i] - start_end[0][gene_i]) * per
    return arr


def dynamic2_sde(dt,t_tot,exp_array,n1,n2,parameters_all,gene_position,n1_array,n2_array,e_mean,save_num):
    dt_use=dt
    t_tot=t_tot
    t_all = np.arange(0, t_tot+dt, dt)
    n_all=len(t_all)
    colnum = np.linspace(0, n_all-1, save_num)
    df = pd.DataFrame(index=range(len(exp_array)), columns=range(len(colnum)))
    df.loc[:,0]=exp_array
    n1_value=exp_array[n1]
    n2_value=exp_array[n2]
    exp=exp_array
    sde_n=len(exp)
    for id in range(1,n_all+1):
        dxdt=np.array(dynamic_tot(last_step=exp,para_total=parameters_all,gene_position=gene_position,dt=dt_use))
        x=dxdt+exp+0.04*e_mean*np.random.normal(0, dt_use, sde_n)
        x[n2]=exp[n2]+n2_array[id-1]
        x[n1]=exp[n1]+n1_array[id-1]
        x[x < 0] = 0
        exp=x.copy()
        if id in colnum:
            inde=np.where(colnum == id)[0][0]
            df.loc[:,inde]=x
    return(df)


def dynamic2_sde_relax(dt,t_tot,exp_array,n1,n2,parameters_all,gene_position,e_mean,save_num):
    dt_use=dt
    t_tot=t_tot
    t_all = np.arange(0, t_tot+dt, dt)
    n_all=len(t_all)
    colnum = np.linspace(0, n_all-1, save_num)
    df = pd.DataFrame(index=range(len(exp_array)), columns=range(len(colnum)))
    df.loc[:,0]=exp_array
    exp=exp_array
    sde_n=len(exp)
    for id in range(1,n_all+1):
        dxdt=np.array(dynamic_tot(last_step=exp,para_total=parameters_all,gene_position=gene_position,dt=dt_use))
        x=dxdt+exp+0.04*e_mean*np.random.normal(0, dt_use, sde_n)
        x[x < 0] = 0
        exp=x.copy()
        if id in colnum:
            inde=np.where(colnum == id)[0][0]
            df.loc[:,inde]=x
    return(df)


def fit_base_pca(baseline_cloud):
    X = baseline_cloud.T.values

    pca = PCA(n_components=2, random_state=0)
    pcs = pca.fit_transform(X)

    pcs_base = pd.DataFrame(
        pcs,
        index=baseline_cloud.columns,
        columns=["PC1", "PC2"]
    )

    explained = pca.explained_variance_ratio_ * 100

    return pca, pcs_base, explained


def run_pair_release(pair, per, init_states, start_end, gene_list, res_final, gene_position, e_mean,
                     dt=1, force_t_tot=20, release_t_tot=2000, save_num=201, seed=2):
    n1 = gene_list.index(pair[0])
    n2 = gene_list.index(pair[1])

    n1_array = make_jump(n1, start_end, per=per, dt=dt, t_tot=force_t_tot)
    n2_array = make_jump(n2, start_end, per=per, dt=dt, t_tot=force_t_tot)

    n_runs = init_states.shape[0]
    n_genes = len(gene_list)
    out = np.zeros((n_genes, n_runs))

    np.random.seed(seed)

    for i in range(n_runs):
        df_signal = dynamic2_sde(
            dt=dt,
            t_tot=force_t_tot,
            exp_array=init_states[i],
            n1=n1,
            n2=n2,
            parameters_all=res_final,
            gene_position=gene_position,
            n1_array=n1_array,
            n2_array=n2_array,
            e_mean=e_mean,
            save_num=save_num
        )

        signal_final = np.array(df_signal.iloc[:, -1]).astype(float)

        df_release = dynamic2_sde_relax(
            dt=dt,
            t_tot=release_t_tot,
            exp_array=signal_final,
            n1=n1,
            n2=n2,
            parameters_all=res_final,
            gene_position=gene_position,
            e_mean=e_mean,
            save_num=save_num
        )

        out[:, i] = np.array(df_release.iloc[:, -1]).astype(float)

    return pd.DataFrame(
        out,
        index=gene_list,
        columns=[f"run_{i+1}" for i in range(n_runs)]
    )


def plot_pair_release(pair, per, baseline_cloud, start_end, gene_list, res_final, gene_position, e_mean,
                      n_runs=200, dt=1, force_t_tot=20, release_t_tot=2000, seed=2,
                      out_dir=None):

    #sample from previous baseline steady states
    init_states = (baseline_cloud.loc[gene_list].sample(n=n_runs, axis=1, replace=False, random_state=seed).T.values.astype(float))

    pca, pcs_base, explained = fit_base_pca(baseline_cloud)

    cloud = run_pair_release(
        pair=pair,
        per=per,
        init_states=init_states,
        start_end=start_end,
        gene_list=gene_list,
        res_final=res_final,
        gene_position=gene_position,
        e_mean=e_mean,
        dt=dt,
        force_t_tot=force_t_tot,
        release_t_tot=release_t_tot,
        seed=seed
    )

    pcs = pca.transform(cloud.T.values)

    label = f"{pair[0]}, {pair[1]}"
    if pair[0] == pair[1]:
        label = pair[0]

    plt.figure(figsize=(6, 5))

    plt.scatter(
        pcs_base["PC1"],
        pcs_base["PC2"],
        s=12,
        alpha=0.30,
        label="Baseline"
    )

    plt.scatter(
        pcs[:, 0],
        pcs[:, 1],
        s=12,
        alpha=1,
        label=label
    )

    plt.xlabel(f"PC1 ({explained[0]:.1f}%)")
    plt.ylabel(f"PC2 ({explained[1]:.1f}%)")
    plt.title(label)
    plt.legend(fontsize=8)
    plt.tight_layout()

    if out_dir is not None:
        os.makedirs(out_dir, exist_ok=True)
        filename = f"{pair[0]}_{pair[1]}.pdf"
        if pair[0] == pair[1]:
            filename = f"{pair[0]}.pdf"
        plt.savefig(os.path.join(out_dir, filename), bbox_inches="tight")

    plt.show()

    return cloud, pca


def run_one_pair(i, j, per, n_runs, baseline_cloud, start_end, gene_list,
                 res_final, gene_position, e_mean, prefix,
                 dt=1, force_t_tot=20, release_t_tot=2000):

    gene1, gene2 = gene_list[i], gene_list[j]
    seed = 0

    init_states = (baseline_cloud.loc[gene_list]
                   .sample(n=n_runs, axis=1, replace=False, random_state=seed)
                   .T.values.astype(float))

    cloud = run_pair_release(pair=[gene1, gene2], per=per, init_states=init_states, start_end=start_end,
        gene_list=gene_list, res_final=res_final, gene_position=gene_position,
        e_mean=e_mean, dt=dt, force_t_tot=force_t_tot, release_t_tot=release_t_tot,
        seed=seed)

    if gene1 == gene2:
        filename = f"{prefix}_{gene1}.csv"
    else:
        filename = f"{prefix}_{gene1}_{gene2}.csv"

    cloud.to_csv(filename)

    return cloud


def run_one_sweep(pair, per, n_runs, baseline_cloud, start_end, gene_list,
                  res_final, gene_position, e_mean, out_path, seed=0,
                  dt=1, force_t_tot=20, release_t_tot=1000):
    init_states = (baseline_cloud.loc[gene_list]
                   .sample(n=n_runs, axis=1, replace=False, random_state=seed)
                   .T.values.astype(float))
    cloud = run_pair_release(
        pair=pair, per=per, init_states=init_states, start_end=start_end,
        gene_list=gene_list, res_final=res_final, gene_position=gene_position,
        e_mean=e_mean, dt=dt, force_t_tot=force_t_tot,
        release_t_tot=release_t_tot, seed=seed)
    cloud.to_csv(out_path)
    return out_path


def scaled_and_inputs(expression, gene_list, ob_newall):
    scaled    = expression.div(expression.max(axis=1), axis=0)   # genes x pseudotime
    e_mean    = scaled.mean(axis=1).values.astype(float)
    start_end = [scaled.iloc[:, 0].values.astype(float),         # first / last pseudotime cols
                 scaled.iloc[:, -1].values.astype(float)]
    gene_to_i = {g: i for i, g in enumerate(gene_list)}
    gene_pos  = [[gene_to_i[r] for r in regs] for regs in ob_newall]
    return scaled, e_mean, start_end, gene_pos


def make_init_states(start_end, n_runs, seed=1):
    lo = np.minimum(start_end[0], start_end[1])
    hi = np.maximum(start_end[0], start_end[1])
    np.random.seed(seed)
    return np.array([np.random.uniform(lo, hi) for _ in range(n_runs)])


def find_wells(pc1_stat, k):
    centers = np.sort(KMeans(n_clusters=k, n_init=10, random_state=0)
                      .fit(pc1_stat.reshape(-1, 1)).cluster_centers_.ravel())
    kde  = gaussian_kde(pc1_stat)
    grid = np.linspace(pc1_stat.min(), pc1_stat.max(), 2000); dens = kde(grid)
    bounds = []
    for a, b in zip(centers[:-1], centers[1:]):
        m = (grid >= a) & (grid <= b)
        bounds.append(grid[m][np.argmin(dens[m])])
    return centers, np.array(bounds), grid, dens


def report_plateau(steps, m, tol=1e-3, label=""):
    d = np.abs(np.diff(m))
    conv = np.where(d < tol)[0]
    step_c = steps[conv[0] + 1] if len(conv) else np.nan
    print(f"[{label}] mean PC1 final={m[-1]:.4f}  "
          f"|Δ| < {tol} first at step ~{step_c:.0f}  "
          f"(value at production release_t_tot=1000: {np.interp(1000, steps, m):.4f}, "
          f"2000: {np.interp(2000, steps, m):.4f})")


# Driving-MSE comparison figure helpers


def prep_pairs(f2, dif2):
    """
    Pair deletions: merge f2 and dif2 on edge identity.
    Keep only positive MSE values <= 10000.
    """
    f2_tmp = f2[["gene1", "gene2", "position1", "position2", "MSE"]].rename(
        columns={"MSE": "mse_f"}
    )

    dif2_tmp = dif2[["gene1", "gene2", "position1", "position2", "MSE"]].rename(
        columns={"MSE": "mse_dif"}
    )

    out = f2_tmp.merge(
        dif2_tmp,
        on=["gene1", "gene2", "position1", "position2"],
        how="inner"
    )

    out = out[
        (out["mse_f"] > 0) &
        (out["mse_dif"] > 0) &
        (out["mse_f"] <= 10000) &
        (out["mse_dif"] <= 10000)
    ].copy()

    out["kind"] = "Pair"
    return out[["mse_f", "mse_dif", "kind"]]


def prep_singles(f1, dif1):
    """
    Single deletions: merge f1 and dif1 on gene identity.
    Keep only positive MSE values <= 10000.
    """
    f1_tmp = f1[["gene", "MSE"]].rename(columns={"MSE": "mse_f"})
    dif1_tmp = dif1[["gene", "MSE"]].rename(columns={"MSE": "mse_dif"})

    out = f1_tmp.merge(dif1_tmp, on="gene", how="inner")

    out = out[
        (out["mse_f"] > 0) &
        (out["mse_dif"] > 0) &
        (out["mse_f"] <= 10000) &
        (out["mse_dif"] <= 10000)
    ].copy()

    out["kind"] = "Single"
    return out[["mse_f", "mse_dif", "kind"]]


def add_condition(df, condition_name):
    df = df.copy()
    df["condition"] = condition_name
    return df
