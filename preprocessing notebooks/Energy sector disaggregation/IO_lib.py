''' Input-Output library
'''

import numpy as np

############################
# Code for computing Direct requirements tables from BEA (See chapter 12 of the BEA IO manual https://www.bea.gov/sites/default/files/methodologies/IOmanual_092906.pdf)
############################

def derive_IO_no_scrap_adjustment(U, V, q, g, fc = None, how = "industry"):
    ''' probability of an application sent to occupation j is successful
    
        return: Z, A, L
        Z ... Industry by industry direct requirements
        A ... Coefficient matrix
        L ... Leontief inverse
    
        All numpy arrays:
        U ... Use table commodity x industry
        V ... Make table industry x commodity
        q ... total commodity output 
        g ... total industry output
        fc .. final consumption commodity x nr of final consumption types
    '''
    n = g.size
    
    # Market share matrix
    D = V * 1/q
    
    # commodity x industry direct requirement
    B = U * 1/g
    
    # scrap adjustment of industry output
    #gstar = g - h
    #p = g / gstar 
    #W = p * D # increase market share with 
    
    
    if (how == "industry"):
        # intermediate consumption
        #Z = np.dot(W, U) # industry x industry direct requirement
        Z = np.dot(D, U)
        # industry x industry direct requirements coefficients
        #A = np.dot(W, B)
        A = Z * 1/g
        # A = Z * 1/g # Alternative, equivalent
    elif (how == "commodity"):
        # commodity x commodity direct requirement
        #Z = np.dot(U * 1/g, W) * q
        Z = np.dot(U * 1/g, D) * q
        # commodity x commodity direct requirement coefficients
        #A = np.dot(B, W)
        A = np.dot(B, D)
        # A = Z * 1/q # Alternative, equivalent
    else:
        Z = None
        A = None
    
    # leontief inverse
    L = np.linalg.inv(np.identity(n) - A)
    
    if fc is not None:
        f = np.dot(D, fc)
    else:
        f = None
    

    return Z, A, L, f


def derive_IO(U, V, h, q, g, fc = None, how = "industry"):
    ''' probability of an application sent to occupation j is successful
    
        return: Z, A, L
        Z ... Industry by industry direct requirements
        A ... Coefficient matrix
        L ... Leontief inverse
    
        All numpy arrays:
        U ... Use table commodity x industry
        V ... Make table industry x commodity
        h ... Scrap and secondhand goods
        q ... total commodity output 
        g ... total industry output
        fc .. final consumption commodity x nr of final consumption types
    '''
    n = g.size
    
    # Market share matrix
    D = V * 1/q
    
    # commodity x industry direct requirement
    B = U * 1/g
    
    # scrap adjustment of industry output
    gstar = g - h
    p = g / gstar 
    #W = p * D # increase market share with 
    W = (p * D.T).T
    
    
    if (how == "industry"):
        # intermediate consumption
        Z = np.dot(W, U) # industry x industry direct requirement
        # industry x industry direct requirements coefficients
        A = np.dot(W, B)
        # A = Z * 1/g # Alternative, equivalent
    elif (how == "commodity"):
        # commodity x commodity direct requirement
        Z = np.dot(U * 1/g, W) * q
        # commodity x commodity direct requirement coefficients
        A = np.dot(B, W)
        # A = Z * 1/q # Alternative, equivalent
    else:
        Z = None
        A = None
    
    # leontief inverse
    L = np.linalg.inv(np.identity(n) - A)
    
    if fc is not None:
        f = np.dot(D, fc)
    else:
        f = None
    

    return Z, A, L, f