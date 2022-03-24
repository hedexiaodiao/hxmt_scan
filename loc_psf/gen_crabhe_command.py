import os
for i in range(0+1,22+1):
    command = 'python genhe.py P0101299008{:02d} note.xml'.format(i)
    os.system(command)