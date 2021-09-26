#!/hxmt/home/lick/soft/anaconda2/bin/python
'''
Find concerned file and filepath relates to Keyword!
usage:
    examlpe=Findfile.findfile()
    list=example()
return a list
'''
import os.path
import sys
import re
class findfile(object):
    def __call__(self,keywords=None,filepath=None):
        filelist = []
        #print 'haha' 
        if keywords==None: 
            print("Find concerned file and filepath relates to Keyword! Please:")
            filekeywd= raw_input('Input keyword to search-->')
            filepath = raw_input('Input file-path to find-->')
        else:
            filekeywd = keywords
            filepath = filepath
        #print filepath 
        for dirpath, dirnames, filenames in os.walk(filepath):
            #print dirpath
            for filename in filenames:
                #print filekeywd,filename
                #if ( filekeywd.lower() in filename.lower()):
                if re.search(filekeywd,filename,re.I) is not None:
                    #print filekeywd,filename
                    filelist.append(os.path.join(dirpath,filename))
        filelist = sorted(filelist)
	#print filelist
        return filelist

