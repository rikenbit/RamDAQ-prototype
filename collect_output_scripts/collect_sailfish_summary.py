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
    for root, dirs, files in os.walk(directory):
        for file in files:
            if fnmatch.fnmatch(file, 'meta_info.json'):
                yield os.path.join(root, file)

def main():
    filelist = []
    for gettxt in fild_all_files(inp_dir):
        print gettxt
        filelist.append(gettxt)
    
    print 'file num of dir : ' + str(len(filelist))
    
    df = pd.DataFrame() 
    for path in filelist:
        with open(path) as fp:
            fname=os.path.split(os.path.dirname(path).replace('/aux',''))[1].replace('sailfish_','').replace('_trim','')
            num_processed=0
            percent_mapped=0
    
            for i, line in enumerate(fp):

                if i == 8:
                    num_processed = line.split(':')[1].split(',')[0]
                elif i == 10:
                    percent_mapped = line.split(':')[1].split(',')[0]
                    break
            
            data = pd.DataFrame([[fname,num_processed,percent_mapped]], \
                                  columns=['name','num_processed','percent_mapped'])
        df = pd.concat([df,data],ignore_index=True)       
    print df.shape
    dataout(df, out_fname)

if __name__ == '__main__':
        main()

print u"....All Done. End of Script"
