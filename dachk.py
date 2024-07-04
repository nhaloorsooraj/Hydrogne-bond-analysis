
# Initialize dictionary to store connection
connection = {}
DA = 0
DAA =0
DDA =0
DDAA =0
# Read the file and store connection
with open("tst.dat", 'r') as r:
    lines = r.readlines()
    for line in lines:
        if line.startswith("--"):
            continue
        parts = line.split()
        if len(parts) == 3:
            M = int(parts[0])
            D1 = int(parts[1])
            D2 = int(parts[2])
            connection[M] = (D1, D2)

# Check for DA condition and verify connection
with open("tst.dat", 'r') as r:
    lines = r.readlines()
    for line in lines:
        if line.startswith("--"):
            continue
        parts = line.split()
        if len(parts) == 3:
            M = int(parts[0])
            D1 = int(parts[1])
            D2 = int(parts[2])

            if D1 ==0 and D2 == 0:
                 continue

            # CHECK FOR DA #
            if (D1 == 0 and D2 != 0) or (D1 != 0 and D2 == 0):
                # Check if M is connected to any of D1 or D2 elsewhere in the file
                if  (D1 in connection and M in connection[D1]) ^ (D2 in connection and M in connection[D2]): # used xor operation " ^ "  to check if M is connected to either D1 or D2 but not both
                    DA = DA+1
            # CHECK FOR DAA #
                elif (D1 in connection and M in connection[D1]) and (D2 in connection and M in connection[D2]):
                        DAA = DAA+1
            # CHECK FOR DDA #
            elif D1 != 0  and D2 != 0 :
                if  (D1 in connection and M in connection[D1]) ^ (D2 in connection and M in connection[D2]): # used xor operation " ^ "  to check if M is connected to either D1 or D2 but not both
                    DDA = DDA+1
            # CHECK FOR DDAA #  
                elif (D1 in connection and M in connection[D1]) and (D2 in connection and M in connection[D2]):
                        DDAA = DDAA+1


print("DA = ",DA )
print("DAA = ",DAA )
print("DDA = ",DDA )
print("DDAA = ",DDAA )

