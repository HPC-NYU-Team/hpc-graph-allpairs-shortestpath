import os


for x in range(100,10100,100):
    command = "python3 createGraph.py "
    command = command+str(x)+" "+"5"
    os.system(command)
