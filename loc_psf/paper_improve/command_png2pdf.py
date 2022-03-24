import glob
import os

png_list = glob.glob('*.png')
jpg_list = glob.glob('*.jpg')
eps_list = glob.glob('*.eps')

tot_list = png_list + jpg_list
tot_list = tot_list + eps_list
print(tot_list)
for i in range(len(tot_list)):
    bf_file = tot_list[i]
    command = 'convert %s -resize 400x %s.pdf'%(bf_file,bf_file[:-4])
    print(command)
    os.system(command)
    print(bf_file,'convert down')