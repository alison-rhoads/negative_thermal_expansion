import numpy as np
import h5py, pprint
import matplotlib.pyplot as plt

# constants
hbar = 6.626e-34 / (2 * np.pi)     
kB= 1.38e-23

def get_cte(infile, C, T, V, iso, tetra):
    
    '''
    infile: gruniesen.hdf5
    C: elastic stiffness matrix (Pa)
    T: temperature (K)
    V: volume of primitive cell (m^3)
    iso and tetra are checks for symmetry of system and take True/False as argument
    '''
    
    # Load Data
    with h5py.File(infile) as f:
        #pprint.pprint(list(f.keys()))
        freqs_thz = f["frequency"][()] # THZ
        G = np.array(f["gruneisen_tensor"][()]) # Gruneisen TENSOR
        weights = f["weight"][()]
    
    freqs = freqs_thz * 1e12  # Hz
    omegas = freqs * 2*np.pi     
    #print(f'Gruniesen Tensor Shape: {G.shape[:]}')  
    n_qpoints = G.shape[0]
    n_modes = G.shape[1]
    #print(f'There are {n_qpoints} q-points and {n_modes} modes')
    C_inv = np.linalg.inv(C) #inverse of stiffness matrix

    # Heat Capacity
    A = hbar * omegas / (kB * T)
    c_mod = kB * A**2 * ( (np.exp(A))/(np.exp(A)-1)**2 ) # modal heat cap
    #print('Modal Heat Capacity Shape:', c_mod.shape)
    cmod_weighted = c_mod * weights[:, None] # weighted   
    #print('Weighted Modal Heat Capacity Shape:', cmod_weighted.shape)
    c_sum = np.sum(cmod_weighted)
    w_sum = np.sum(weights)
    c_tot = c_sum / w_sum

    # Gruneisen and CTE Calculations
    g_mods = []
    g_mod_xxs = []
    g_mod_zzs = []
    
    if iso == True:
        gc_array = [] # array of products of modal grun and heat cap
        #print('Proceeding for ISOTROPIC/CUBIC symmetry...')
        for q in range(n_qpoints): 
            for m in range(n_modes):

                c = float(cmod_weighted[q, m]) # element of heat cap matrix (float)
                #print('cs shape', cs)
                G_mod = G[q, m]  # second rank 3x3 matrix 
                #print(G_mod)
                #print('G_mod shape',G_mod.shape)
                g_mod = (1/3) * np.trace(G_mod) # scalar gruneisen parameter for specific (q,m)

                gc = g_mod * c # product of modal grun and heat cap
                gc_array.append(gc) 
                
        g_tot = np.sum(gc_array) / c_sum

        cte_lin = (g_tot * c_tot) / (V*(C[0,0] + 2*C[0,1]))
        cte_vol = 3 * cte_lin

        return cte_lin, cte_vol
        sys.exit(0)     

    if tetra == True:
        #print('Proceeding for TETRAHEDRAL symmetry...')
    
        gc_xx_array = []
        gc_zz_array = []
        gc_yy_array = []
        gc_yz_array = []
        gc_xz_array = []
        gc_xy_array = []
        
        for q in range(n_qpoints):       
            for m in range(n_modes):
                c = float(cmod_weighted[q, m]) # element of heat cap matrix (float)
                G_mod = G[q, m]  # second rank 3x3 matrix
                #print(G_mod)
                #TODO: later have this in voigt vector form rather than component wise
                gc_xx = c * G_mod[0, 0]  # xx
                gc_yy = c * G_mod[1, 1]  # yy
                gc_zz = c * G_mod[2, 2]  # zz
                gc_yz = c * G_mod[1, 2]  # yz
                gc_xz = c * G_mod[0, 2]  # xz
                gc_xy = c * G_mod[0, 1]  # xy
              
                gc_xx_array.append(gc_xx)
                gc_zz_array.append(gc_zz)
                gc_yy_array.append(gc_yy)
                gc_yz_array.append(gc_yz)
                gc_xz_array.append(gc_xz)
                gc_xy_array.append(gc_xy)
                
        g_tot_xx = np.sum(gc_xx_array) / (c_sum)
        g_tot_zz = np.sum(gc_zz_array) / (c_sum)

        gtot_vec = np.array([ 
            np.sum(gc_xx_array) / (c_sum),
            np.sum(gc_yy_array) / (c_sum),
            np.sum(gc_zz_array) / (c_sum),
            np.sum(gc_yz_array) / (c_sum),
            np.sum(gc_xz_array) / (c_sum),
            np.sum(gc_xy_array) / (c_sum)
        ])

        cte_vec = (c_tot / V) * (C_inv @ gtot_vec)
        
        cte_lin_xx = cte_vec[0]
        cte_lin_zz = cte_vec[2]
        cte_vol = 2*cte_lin_xx + cte_lin_zz

        return cte_lin_xx, cte_lin_zz, cte_vol, cte_vec 
        sys.exit(0)
