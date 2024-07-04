
homo_WW = 0
hetero_WE = 0
homo_EE = 0
free_W =0
bound_W = 0
dang_W =0
nwat = 1600
norm_fact = 1600.0 * 10001

with open ("tst.dat",'r') as r:
    lines = r.readlines()
    for line in lines:
        if line.startswith("--"):
            continue
        parts = line.split()
        if len(parts) == 3 :
            h1w = int(parts[1])
            h2w = int(parts[2])
### CHECK FOR FREE WATER ###
            if h1w == 0 and h2w ==0:
                free_W = free_W+1                      
### CHECK FOR DANGLING WATER ### 
            elif h1w == 0 and h2w != 0:
                dang_W = dang_W+1
            elif h1w != 0 and h2w == 0:
                dang_W = dang_W+1            
### CHECK FOR BOUND WATER ### 
            elif h1w != 0 and h2w != 0:
                bound_W = bound_W+1
            ### CHECK FOR HOMO WATER ### 
                if h1w <= nwat and h2w <= nwat:
                    homo_WW = homo_WW+1
            ### CHECK FOR HETERO WATER ### 
                elif h1w <= nwat and h2w > nwat:
                    hetero_WE = hetero_WE+1
                elif h1w > nwat and h2w <= nwat:
                    hetero_WE = hetero_WE+1    
            ### CHECK FOR HETERO WATER ### 
                elif h1w > nwat and h2w > nwat:
                    homo_EE = homo_EE+1    

         
free_W = free_W / norm_fact
dang_W = dang_W / norm_fact
bound_W = bound_W / norm_fact
homo_WW = homo_WW / norm_fact
hetero_WE = hetero_WE / norm_fact
homo_EE = homo_EE / norm_fact


print ("Number of Free bonds in water = ", free_W)
print ("Number of Dangling bonds in water = ", dang_W)
print ("Number of Bound bonds in water = ", bound_W)
print ("Number of Homo bonds with Water = ", homo_WW)
print ("Number of Hetero bonds with Ethanol = ", hetero_WE)
print ("Number of Homo bonds with Ethanol = ", homo_EE)






