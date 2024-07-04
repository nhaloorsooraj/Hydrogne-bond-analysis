#folder 3
import MDAnalysis as mda
import numpy as np

R_CUT_OwOw = 0.345
R_CUT_OwHw = 0.25
R_CUT_OwOe = 0.354
THETA_CUT  = round(np.degrees(np.pi / 6.00), 3)

ethend = #end of ethanol
solend = #end of sol



u = mda.Universe("md.gro", "md.trr", convert_units=False)
nwat = u.select_atoms("resname SOL").n_residues
neth = u.select_atoms("resname ETH").n_residues
tframe = len(u.trajectory) # number of frames



print("")
print("")
print("****** INFO ******")
print("")
print("")
print("Number of frames found in trajectory = ", tframe )
print("")
print("Number of atoms = ", len(u.atoms))
print("")
print("Number of Residues = ", len(u.residues))
print("")
residue_types = list(set(u.residues.resnames))
print("Unique residue types =", residue_types)
print ("number of ETH = ", neth)
print ("number of SOL = ", nwat)

print("")
print("-- Box Dimension --")
print("")
print("X= ", u.dimensions[0])
print("Y= ", u.dimensions[1])
print("Z= ", u.dimensions[2])
print("")
print("-- -- -- -- -- -- --")




box_x = u.dimensions[0]
box_y = u.dimensions[1]
box_z = u.dimensions[2]

SOL_OXY_POS = []
SOL_HYD1_POS = []
SOL_HYD2_POS = []
ETH_OXY_POS = []
ETH_HYD_POS = []
with open('hflag.dat', 'w') as f:
    f.write("---Hydrogen bond analysis--- \n")
    f.write("---------------------------------------------------\n")
    f.write("---Molecule----------H1w-----------------H2w---\n")

    
    for ts in u.trajectory[0:tframe+1]: # start 0 ; end lastframe ; step 1 
        SOL_OXY_POS.clear()
        SOL_HYD1_POS.clear()
        SOL_HYD2_POS.clear()
        ETH_OXY_POS.clear()
        ETH_HYD_POS.clear()
        f.write(f"---FRAME:{ts.frame}---\n")

        COORDINATE = ts.positions
        #-------- GETTING COORDINATES OF WATER ATOMS --------
        #--- The index of array will  be the molecule index ---------
        for i in range(ethend+1,solend+1, 4):  # start : soldend+1, end : solend, step : 4
            SOL_OXY_POS.append(COORDINATE[i])
            SOL_HYD1_POS.append(COORDINATE[i + 1])
            SOL_HYD2_POS.append(COORDINATE[i + 2])

        #-------- GETTING COORDINATES OF ETH ATOMS --------
        #--- The index of array will  be the molecule index ---------
        for k in range(1, ethend+1, 9):  # start : 1, end : ethend, step : 9
            ETH_OXY_POS.append(COORDINATE[k + 7])
            ETH_HYD_POS.append(COORDINATE[k + 8])
            
        SOL_OXY_NUM  = len(SOL_OXY_POS) # Number of WATER OXYGEN
        SOL_HYD1_NUM = len(SOL_HYD1_POS) # Number of WATER HYDROGEN1
        SOL_HYD2_NUM = len(SOL_HYD2_POS) # Number of WATER HYDROGEN2

        ETH_OXY_NUM  = len(ETH_OXY_POS) # Number of ETHANOL OXYGEN
        ETH_HYD_NUM = len(ETH_HYD_POS) # Number of ETHANOL HYDROGEN
    #-------------------------------------------------------------------------------------------------------------
    # 1. Calculate the vector distance between i th and j th Oxygen in Water : Vect Oi --> Oj (Vect_0)
    #-------------------------------------------------------------------------------------------------------------        
        for i in range(SOL_OXY_NUM):
            flag_H1w = 0
            flag_H2w = 0
            for j in range (SOL_OXY_NUM):
                if i != j :
                     SOL_OwOw_dx = SOL_OXY_POS[j][0] - SOL_OXY_POS[i][0]
                     SOL_OwOw_dy = SOL_OXY_POS[j][1] - SOL_OXY_POS[i][1]
                     SOL_OwOw_dz = SOL_OXY_POS[j][2] - SOL_OXY_POS[i][2]
    
                     SOL_OwOw_dx  = SOL_OwOw_dx - (round(SOL_OwOw_dx/box_x))*box_x  
                     SOL_OwOw_dy  = SOL_OwOw_dy - (round(SOL_OwOw_dy/box_y))*box_y 
                     SOL_OwOw_dz  = SOL_OwOw_dz - (round(SOL_OwOw_dz/box_z))*box_z  
                      
                     mod_vecW = np.sqrt( SOL_OwOw_dx**2 + SOL_OwOw_dy**2 + SOL_OwOw_dz**2 )
                    # print("modvectW = ",mod_vecW)
    # Checking with cut off distance btw two molecule, create unit vectors   
                     if mod_vecW < R_CUT_OwOw and mod_vecW != 0:
                         SOL_OwOw_dx = SOL_OwOw_dx/mod_vecW 
                         SOL_OwOw_dy = SOL_OwOw_dy/mod_vecW
                         SOL_OwOw_dz = SOL_OwOw_dz/mod_vecW
     
    #-------------------------------------------------------------------------------------------------------------
    # 2. Calculate the vector distance between i th O  and j th H1 in Water : Oi ---> H1i (Vect_1)
    #-------------------------------------------------------------------------------------------------------------
    
                         SOL_H1wOw_dx = SOL_HYD1_POS[i][0] - SOL_OXY_POS[i][0]
                         SOL_H1wOw_dy = SOL_HYD1_POS[i][1] - SOL_OXY_POS[i][1]
                         SOL_H1wOw_dz = SOL_HYD1_POS[i][2] - SOL_OXY_POS[i][2]
              
                         SOL_H1wOw_dx  = SOL_H1wOw_dx - (round(SOL_H1wOw_dx/box_x))*box_x  
                         SOL_H1wOw_dy  = SOL_H1wOw_dy - (round(SOL_H1wOw_dy/box_y))*box_y 
                         SOL_H1wOw_dz  = SOL_H1wOw_dz - (round(SOL_H1wOw_dz/box_z))*box_z   
                         mod_vecW = np.sqrt( SOL_H1wOw_dx**2 + SOL_H1wOw_dy**2 + SOL_H1wOw_dz**2 )
    # Checking with cut off distance, create unit vectors   
    
                         SOL_H1wOw_dx = SOL_H1wOw_dx/mod_vecW
                         SOL_H1wOw_dy = SOL_H1wOw_dy/mod_vecW
                         SOL_H1wOw_dz = SOL_H1wOw_dz/mod_vecW 
    
    #-------------------------------------------------------------------------------------------------------------
    # 3. Calculate the vector distance between i th O and j th H2 in Water : Oi --> H2i (Vect_2)
    #-------------------------------------------------------------------------------------------------------------
    
                         SOL_H2wOw_dx = SOL_HYD2_POS[i][0] - SOL_OXY_POS[i][0]
                         SOL_H2wOw_dy = SOL_HYD2_POS[i][1] - SOL_OXY_POS[i][1]
                         SOL_H2wOw_dz = SOL_HYD2_POS[i][2] - SOL_OXY_POS[i][2]
    
                         SOL_H2wOw_dx  = SOL_H2wOw_dx - (round(SOL_H2wOw_dx/box_x))*box_x  
                         SOL_H2wOw_dy  = SOL_H2wOw_dy - (round(SOL_H2wOw_dy/box_y))*box_y 
                         SOL_H2wOw_dz  = SOL_H2wOw_dz - (round(SOL_H2wOw_dz/box_z))*box_z   
                         mod_vecW = np.sqrt( SOL_H2wOw_dx**2 + SOL_H2wOw_dy**2 + SOL_H2wOw_dz**2 )
    
    # Checking with cut off distance, create unit vectors   
    
                         SOL_H2wOw_dx = SOL_H2wOw_dx/mod_vecW
                         SOL_H2wOw_dy = SOL_H2wOw_dy/mod_vecW
                         SOL_H2wOw_dz = SOL_H2wOw_dz/mod_vecW
    
    #-------------------------------------------------------------------------------------------------------------
    # 4. Calculate the vector distance between i th H1 and j th O  in Water : H1i --> Oj (Vect_4)
    #-------------------------------------------------------------------------------------------------------------
    
                         SOL_OwH1w_dx = SOL_OXY_POS[j][0] - SOL_HYD1_POS[i][0]
                         SOL_OwH1w_dy = SOL_OXY_POS[j][1] - SOL_HYD1_POS[i][1]
                         SOL_OwH1w_dz = SOL_OXY_POS[j][2] - SOL_HYD1_POS[i][2]
    
                         SOL_OwH1w_dx  = SOL_OwH1w_dx - (round(SOL_OwH1w_dx/box_x))*box_x  
                         SOL_OwH1w_dy  = SOL_OwH1w_dy - (round(SOL_OwH1w_dy/box_y))*box_y 
                         SOL_OwH1w_dz  = SOL_OwH1w_dz - (round(SOL_OwH1w_dz/box_z))*box_z   
                         mod_vecW = np.sqrt( SOL_OwH1w_dx**2 + SOL_OwH1w_dy**2 + SOL_OwH1w_dz**2 )
    
    # Checking with cut off distance, create unit vectors   
                         SOL_OwH1w_dx = SOL_OwH1w_dx/mod_vecW
                         SOL_OwH1w_dy = SOL_OwH1w_dy/mod_vecW
                         SOL_OwH1w_dz = SOL_OwH1w_dz/mod_vecW
    
                         DIST_SOL_OjH1i = mod_vecW # Distance Between H1i -- Oj
    
    #-------------------------------------------------------------------------------------------------------------
    # 5. Calculate the vector distance between i th H2 and j th O in Water : H2i --> Oj (Vect_5)
    #-------------------------------------------------------------------------------------------------------------
    
                         SOL_OwH2w_dx = SOL_OXY_POS[j][0] - SOL_HYD2_POS[i][0]
                         SOL_OwH2w_dy = SOL_OXY_POS[j][1] - SOL_HYD2_POS[i][1]
                         SOL_OwH2w_dz = SOL_OXY_POS[j][2] - SOL_HYD2_POS[i][2]
    
                         SOL_OwH2w_dx  = SOL_OwH2w_dx - (round(SOL_OwH2w_dx/box_x))*box_x  
                         SOL_OwH2w_dy  = SOL_OwH2w_dy - (round(SOL_OwH2w_dy/box_y))*box_y 
                         SOL_OwH2w_dz  = SOL_OwH2w_dz - (round(SOL_OwH2w_dz/box_z))*box_z   
                         mod_vecW = np.sqrt( SOL_OwH2w_dx**2 + SOL_OwH2w_dy**2 + SOL_OwH2w_dz**2 )
    
    # Checking with cut off distance, create unit vectors   
                         SOL_OwH2w_dx = SOL_OwH2w_dx/mod_vecW
                         SOL_OwH2w_dy = SOL_OwH2w_dy/mod_vecW
                         SOL_OwH2w_dz = SOL_OwH2w_dz/mod_vecW
    
                         DIST_SOL_OjH2i = mod_vecW  # Distance Between H2i --- Oj
    
    ## ----------- FINDING DOT PRODUCT -----------------------------------------
    
    
    #    find all H1i ---- Oi ----- Oj ANGLE 
    
                 
                         ANGLE_OOH1 = np.degrees(np.arccos(SOL_OwOw_dx * SOL_H1wOw_dx + 
                                                           SOL_OwOw_dy * SOL_H1wOw_dy +
                                                           SOL_OwOw_dz * SOL_H1wOw_dz))
                         
    #    find all H2i ---- Oi ----- Oj ANGLE 
                         
                         ANGLE_OOH2 = np.degrees(np.arccos(SOL_OwOw_dx * SOL_H2wOw_dx +
                                                           SOL_OwOw_dy * SOL_H2wOw_dy +
                                                           SOL_OwOw_dz * SOL_H2wOw_dz))
    
    
                         if (ANGLE_OOH1 < THETA_CUT ) and (DIST_SOL_OjH1i < R_CUT_OwHw ):
                             flag_H1w = j
    
                         elif (ANGLE_OOH2 < THETA_CUT) and (DIST_SOL_OjH2i < R_CUT_OwHw):
                             flag_H2w = j                         
          
    
    #######################################################################################################################################
                     ###############################          WATER , ETHANOL       ###############################
    #######################################################################################################################################
    
    #-------------------------------------------------------------------------------------------------------------
    # 6. Calculate the vector distance between i th OXYGEN in WATER and k th OXYGEN in ETHANOL : Vect Owi --> Oek
    #-------------------------------------------------------------------------------------------------------------
    
            for k in range (ETH_OXY_NUM):
              
                     ETH_SOL_OwOe_dx = ETH_OXY_POS[k][0] - SOL_OXY_POS[i][0]
                     ETH_SOL_OwOe_dy = ETH_OXY_POS[k][1] - SOL_OXY_POS[i][1]
                     ETH_SOL_OwOe_dz = ETH_OXY_POS[k][2] - SOL_OXY_POS[i][2]
    
                     ETH_SOL_OwOe_dx  = ETH_SOL_OwOe_dx - (round(ETH_SOL_OwOe_dx/box_x))*box_x  
                     ETH_SOL_OwOe_dy  = ETH_SOL_OwOe_dy - (round(ETH_SOL_OwOe_dy/box_y))*box_y 
                     ETH_SOL_OwOe_dz  = ETH_SOL_OwOe_dz - (round(ETH_SOL_OwOe_dz/box_z))*box_z  
                      
                     mod_vecE = np.sqrt( ETH_SOL_OwOe_dx**2 + ETH_SOL_OwOe_dy**2 + ETH_SOL_OwOe_dz**2 )
    
    # Checking with cut off distance btw two molecule, create unit vectors   
                     if mod_vecW < R_CUT_OwOe:
                         ETH_SOL_OwOe_dx = ETH_SOL_OwOe_dx/mod_vecE 
                         ETH_SOL_OwOe_dy = ETH_SOL_OwOe_dy/mod_vecE
                         ETH_SOL_OwOe_dz = ETH_SOL_OwOe_dz/mod_vecE    
    
    #----------------------------------------------------------------------------------------------------------------
    # 7. Calculate the vector distance between k th OXYGEN in ETHANOL and i th HYDROGEN1 in WATER : Vect Oek --> H1wi
    #----------------------------------------------------------------------------------------------------------------
                         ETH_SOL_OeH1w_dx = SOL_HYD1_POS[i][0] - ETH_OXY_POS[k][0]
                         ETH_SOL_OeH1w_dy = SOL_HYD1_POS[i][1] - ETH_OXY_POS[k][1]
                         ETH_SOL_OeH1w_dz = SOL_HYD1_POS[i][2] - ETH_OXY_POS[k][2]
              
                         ETH_SOL_OeH1w_dx  = ETH_SOL_OeH1w_dx - (round(ETH_SOL_OeH1w_dx/box_x))*box_x  
                         ETH_SOL_OeH1w_dy  = ETH_SOL_OeH1w_dy - (round(ETH_SOL_OeH1w_dy/box_y))*box_y 
                         ETH_SOL_OeH1w_dz  = ETH_SOL_OeH1w_dz - (round(ETH_SOL_OeH1w_dz/box_z))*box_z   
                         mod_vecE = np.sqrt( ETH_SOL_OeH1w_dx**2 + ETH_SOL_OeH1w_dy**2 + ETH_SOL_OeH1w_dz**2 )
    # Checking with cut off distance, create unit vectors   
    
                         ETH_SOL_OeH1w_dx = ETH_SOL_OeH1w_dx/mod_vecE
                         ETH_SOL_OeH1w_dy = ETH_SOL_OeH1w_dy/mod_vecE
                         ETH_SOL_OeH1w_dz = ETH_SOL_OeH1w_dz/mod_vecE
    
    
                         DIST_SOL_OkH1i = mod_vecE # Distance Between H1i --- Oek
    
    
    #----------------------------------------------------------------------------------------------------------------
    # 8. Calculate the vector distance between k th OXYGEN in ETHANOL and i th HYDROGEN2 in WATER : Vect Oek --> H2wi
    #----------------------------------------------------------------------------------------------------------------
    
    
                         ETH_SOL_OeH2w_dx = SOL_HYD2_POS[i][0] - ETH_OXY_POS[k][0]
                         ETH_SOL_OeH2w_dy = SOL_HYD2_POS[i][1] - ETH_OXY_POS[k][1]
                         ETH_SOL_OeH2w_dz = SOL_HYD2_POS[i][2] - ETH_OXY_POS[k][2]
              
                         ETH_SOL_OeH2w_dx  = ETH_SOL_OeH2w_dx - (round(ETH_SOL_OeH2w_dx/box_x))*box_x  
                         ETH_SOL_OeH2w_dy  = ETH_SOL_OeH2w_dy - (round(ETH_SOL_OeH2w_dy/box_y))*box_y 
                         ETH_SOL_OeH2w_dz  = ETH_SOL_OeH2w_dz - (round(ETH_SOL_OeH2w_dz/box_z))*box_z   
                         mod_vecE = np.sqrt( ETH_SOL_OeH2w_dx**2 + ETH_SOL_OeH2w_dy**2 + ETH_SOL_OeH2w_dz**2 )
    # Checking with cut off distance, create unit vectors   
    
                         ETH_SOL_OeH2w_dx = ETH_SOL_OeH2w_dx/mod_vecE
                         ETH_SOL_OeH2w_dy = ETH_SOL_OeH2w_dy/mod_vecE
                         ETH_SOL_OeH2w_dz = ETH_SOL_OeH2w_dz/mod_vecE
    
    
                         DIST_SOL_OkH2i = mod_vecE # Distance Between H2i --- Oek
    
    ## ----------- FINDING DOT PRODUCT -----------------------------------------
    
    
    #    find all Oe ---- Ow ----- H1wi ANGLE 
    
                 
                         ANGLE_OeOwH1 = np.degrees(np.arccos(ETH_SOL_OwOe_dx * SOL_H1wOw_dx + 
                                                             ETH_SOL_OwOe_dy * SOL_H1wOw_dy +
                                                             ETH_SOL_OwOe_dz * SOL_H1wOw_dz))
                         
    #    find all Oe ---- Oi ----- H2wi ANGLE 
                         
                         ANGLE_OeOwH2 = np.degrees(np.arccos(ETH_SOL_OwOe_dx * SOL_H2wOw_dx +
                                                             ETH_SOL_OwOe_dy * SOL_H2wOw_dy +
                                                             ETH_SOL_OwOe_dz * SOL_H2wOw_dz))
                         
    
                         if (ANGLE_OeOwH1 < THETA_CUT ) and (DIST_SOL_OkH1i < R_CUT_OwOe ):
                             flag_H1w = k + SOL_OXY_NUM
    
                         elif (ANGLE_OeOwH2 < THETA_CUT) and (DIST_SOL_OkH2i < R_CUT_OwOe):
                             flag_H2w = k + SOL_OXY_NUM                     
             
            f.write(f' {i+1:<5d}               {flag_H1w:<5d}               {flag_H2w:<5d}\n')
           # f.write(f'"Et" {i+SOL_OXY_NUM}             {flag_H1e}            {flag_H2e}\n')
            f.write('\n')
