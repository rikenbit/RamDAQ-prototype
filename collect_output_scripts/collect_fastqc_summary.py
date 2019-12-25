import sys
import os
print "python :" + sys.version
import ConfigParser
import pandas as pd
import fnmatch

argvs = sys.argv
argc = len(argvs)

if (argc < 2):
    print 'Usage: # python %s input_dir out_name' % argvs[0]
    quit()

inp_dir = argvs[1]
out_fname = argvs[2]

def makeDir(dname):
    if os.path.exists(dname) is False:
        os.mkdir(dname)
        print '%s (dir) created.' % dname
    else:
        print '%s (dir) is already exists.' % dname

def dataout(data, filename):
    data.to_csv(filename + '.txt', sep="\t")

def fild_all_files(directory):
    for root, dirs, files in os.walk(directory, onerror=error_cb):
        for file in files:
            if fnmatch.fnmatch(file, 'fastqc_data.txt'):
                yield os.path.join(root, file)

def error_cb(e):
    print 'could not open "{0}" [{1}]: {2}'.format(e.filename, e.errno, e.strerror)
    raise oserror(e)

def main():
    filelist = []
    for gettxt in fild_all_files(inp_dir):
        print gettxt
        filelist.append(gettxt)
    
    print 'file num of dir : ' + str(len(filelist))
    
    df = pd.DataFrame() 
    for path in filelist:
        with open(path) as fp:
            fname=""
            totalseq=0
            perGC=0
            for i, line in enumerate(fp):
                if i == 3:
                    fname = line.split()[1].split(".")[0]
                elif i == 6:
                    totalseq = line.split()[2]
                elif i == 9:
                    perGC = line.split()[1]
                    break
            data = pd.DataFrame([[fname,totalseq,perGC]],columns=['name','totalseq','perGC'])
            #print data
        df = pd.concat([df,data],ignore_index=True)    
    print df.shape
    dataout(df, out_fname)

if __name__ == '__main__':
        main()

print u"....All Done. End of Script"
