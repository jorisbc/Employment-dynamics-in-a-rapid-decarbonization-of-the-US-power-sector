''' Input-Output library
@Joris Bucker
'''

import numpy as np

############################
# Code for computing Direct requirements tables from BEA (See chapter 12 of the BEA IO manual https://www.bea.gov/sites/default/files/methodologies/IOmanual_092906.pdf)
############################

def derive_IO_domestic_no_scrap_adjustment(U, W, V, q, g, fc = None, how = "industry"):
    '''     
    See domestic requirement creation: https://apps.bea.gov/scb/pdf/2017/03%20March/0317_introducing_domestic_requirement_tables.pdf
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
    
    U_dom = U-W
    U_dom[U_dom<0] = 0 # no negative values

    
    # Market share matrix
    D = V * 1/q
    
    # commodity x industry direct requirement
    B = U_dom * 1/g
    
    # scrap adjustment of industry output
    #gstar = g - h
    #p = g / gstar 
    #W = p * D # increase market share with 
    
    
    if (how == "industry"):
        # intermediate consumption
        #Z = np.dot(W, U_dom) # industry x industry direct requirement
        Z = np.dot(D, U_dom)
        # industry x industry direct requirements coefficients
        #A = np.dot(W, B)
        A = Z * 1/g
        # A = Z * 1/g # Alternative, equivalent
    elif (how == "commodity"):
        # commodity x commodity direct requirement
        #Z = np.dot(U_dom * 1/g, W) * q
        Z = np.dot(U_dom * 1/g, D) * q
        # commodity x commodity direct requirement coefficients
        #A = np.dot(B, W)
        A = np.dot(B, D)
        # A = Z * 1/q # Alternative, equivalent
    else:
        Z = None
        A = None

    # total imports per industry inputs
    Z = np.r_[Z, [W.sum(axis=0)]]
    
    # fraction imports per industry inputs
    imp_frac = W.sum(axis=0) * 1/g
    A = np.r_[A, [imp_frac]]
    A = np.c_[A, np.zeros(A.shape[0])]
    
    # leontief inverse
    L = np.linalg.inv(np.identity(n+1) - A)
    
    if fc is not None:
        f = np.dot(D, fc)
    else:
        f = None
    
    return Z, A[:, :-1], L[:, :-1], f


def derive_IO_no_scrap_adjustment(U, V, q, g, fc = None, how = "industry"):
    '''     
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
    '''     
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

def derive_IO_domestic(U, W, V, h, inon, q, g, fc = None, W_fc = None, how = "industry", final_scrap = True):
    '''     
    See domestic requirement creation: https://apps.bea.gov/scb/pdf/2017/03%20March/0317_introducing_domestic_requirement_tables.pdf
        return: Z, A, L
        Z ... Industry by industry direct requirements
        A ... Coefficient matrix
        L ... Leontief inverse
    
        All numpy arrays:
        U ... Use table commodity x industry
        V ... Make table industry x commodity
        q ... total commodity output 
        h ... Scrap and secondhand goods
        inon  Non-standard imports
        g ... total industry output
        fc .. final consumption commodity x nr of final consumption types
        W_fc  Imports in final consumption commodity x nr of final consumption types
    '''
    n = g.size
    
    U_dom = U-W
    if final_scrap:    # final is scrap
        scrap_factor = ((U_dom.sum(axis=0) - U_dom[-1, :]) / U_dom.sum(axis=0))
        U_dom = U_dom / scrap_factor
        U_dom = U_dom[:-1, :]
    
    U_dom[U_dom<0] = 0 # no negative values
    
    # Market share matrix
    D = V * 1/q
    
    # commodity x industry direct requirement
    B = U_dom * 1/g
    
    # scrap adjustment of industry output
    gstar = g - h
    p = g / gstar 
    SC = p * D # increase market share with scrap
    
    
    if (how == "industry"):
        # intermediate consumption
        Z = np.dot(SC, U_dom) # industry x industry direct requirement
        #Z = np.dot(D, U_dom)
        # industry x industry direct requirements coefficients
        A = np.dot(SC, B)
        #A = Z * 1/g
        # A = Z * 1/g # Alternative, equivalent
    elif (how == "commodity"):
        # commodity x commodity direct requirement
        Z = np.dot(U_dom * 1/g, W) * q
        #Z = np.dot(U_dom * 1/g, D) * q
        # commodity x commodity direct requirement coefficients
        A = np.dot(B, SC)
        #A = np.dot(B, D)
        # A = Z * 1/q # Alternative, equivalent
    else:
        Z = None
        A = None

    # total imports per industry inputs
    Z = np.r_[Z, [W.sum(axis=0) + inon]]
    
    # fraction imports per industry inputs
    # imports of all goods + non-standard imports
    imp_frac = (W.sum(axis=0) + inon) * 1/g
    A = np.r_[A, [imp_frac]]
    A = np.c_[A, np.zeros(A.shape[0])]
    
    # leontief inverse
    L = np.linalg.inv(np.identity(n+1) - A)
    
    if fc is not None:
        fc = fc - W_fc
        if final_scrap:    # final is scrap
            scrap_factor_fc = (fc.sum(axis=0) - fc[-1, :]) / fc.sum(axis=0)
            fc = fc / scrap_factor_fc
            fc = fc[:-1, :]
        f = np.dot(D, fc)
    else:
        f = None
    
    return Z, A[:, :-1], L[:, :-1], f, p

def derive_import_perc(U, W, V, h, inon, q, g, fc = None, W_fc = None, how = "industry"):
    '''     
    See domestic requirement creation: https://apps.bea.gov/scb/pdf/2017/03%20March/0317_introducing_domestic_requirement_tables.pdf
        return: Z, A, L
        Z ... Industry by industry direct requirements
        A ... Coefficient matrix
        L ... Leontief inverse
    
        All numpy arrays:
        U ... Use table commodity x industry
        V ... Make table industry x commodity
        q ... total commodity output 
        h ... Scrap and secondhand goods
        inon  Non-standard imports
        g ... total industry output
        fc .. final consumption commodity x nr of final consumption types
        W_fc  Imports in final consumption commodity x nr of final consumption types
    '''
    n = g.size
    
    # Market share matrix
    D = V * 1/q
    
    # industry by industry
    Z = np.dot(D, U)
    Z_imp = np.dot(D, W)
    
    tot_use_by_ind = Z.sum(axis=1)
    tot_imp_by_ind = Z_imp.sum(axis=1)
    
    import_perc = tot_imp_by_ind / tot_use_by_ind
    
    return import_perc

