import os
for i in range(0+1,22+1):
    command = 'python genle.py P0101299008{:02d}'.format(i)
    os.system(command)