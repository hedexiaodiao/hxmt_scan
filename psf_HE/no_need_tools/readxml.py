#!//home/hxmt/lick/soft/anaconda2/bin/python
'''
Copyright 2016, HXMT Ground System
All right reserved
File Name: hxmtpsf.py for sample

.......

Author: lick
Version: 0.1
Date:  2016-12-11
Email: lick@ihep.ac.cn
......
'''
from xml.dom import minidom 
import sys
class loadDom:
    '''this is load xml class'''
    def __init__(self,config):
        self.file = config
    def readMiniDom(self,child): 
        result=[]
        try:
            doc = minidom.parse( self.file ) 
            # get root element: <employees/> 
            root = doc.documentElement 
        except Exception,e:
            print ("the input xml file error")
            sys.exit(1)
        # get all children elements: <employee/> <employee/> 
        babys = root.getElementsByTagName( child ) 
        for baby in babys: 
            print ( " ------------------------------------------- " ) 
        # element name : employee 
            print (baby.nodeName,baby.nodeValue,baby.nodeType) 
        # element xml content : <employee><name>windows</name><age>20</age></employee> 
        # basically equal to toprettyxml function 
        #    print (baby.toxml()) 
        #    tempNode = baby.getElementsByTagName( "templatepath" )[0] 
        #    print (tempNode.childNodes) 
        #    print (tempNode.nodeName +  ":"  + tempNode.childNodes[0].nodeValue) 
        #    inNode = baby.getElementsByTagName( "inpath" )[0] 
        #    print (inNode.childNodes) 
        #    print (inNode.nodeName +  ":"  + inNode.childNodes[0].nodeValue) 
        #    print ( " ------------------------------------------- " ) 
            for i,n in enumerate(baby.childNodes): 
                try:
                    if n.hasChildNodes() == True:
                        for j,m in enumerate(n.childNodes):
                            #print m.nodeType 
                            print m.nodeName 
                            print m.data 
                            print (i,n,j,m)
                            print ("**** *****************")
                            result.append(m.data.encode().replace("\n",""))
                except Exception,e:
                    #print(i,"have no child nodes")
                    continue
        return result
    def getTagText(self, child):  
        try:
            doc = minidom.parse( self.file ) 
            # get root element: <employees/> 
            root = doc.documentElement 
        except Exception,e:
            print ("the input xml file error")
            sys.exit(1)
        try:
            nodes = root.getElementsByTagName(child)  
        except Exception,e:
            pt = "the name '"+child+"' is not in the xml."
            print (pt)
            sys.exit(1)
        rc = ""  
        for node in nodes:
            for n in node.childNodes:  
                if n.nodeType in ( n.TEXT_NODE, n.CDATA_SECTION_NODE):  
                    rc = rc + n.data  
        return rc 
if __name__ == "__main__":
    test = loadDom("config.xml")
    r2 = test.readMiniDom("infilename")
    r = test.readMiniDom("PATH")
    r3 = test.getTagText("infilenamelist")
    print r2
    print r
    print r3

