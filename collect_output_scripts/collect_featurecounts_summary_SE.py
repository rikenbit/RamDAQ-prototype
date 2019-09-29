import sys
import os
print "python :" + sys.version
import ConfigParser
import pandas as pd
import fnmatch

argvs = sys.argv
argc = len(argvs)

if (argc < 3):
    print 'Usage: # python %s input_dir out_name pipeline_class' % argvs[0]
    quit()

inp_dir = argvs[1]
out_fname = argvs[2]
pipeline_class = argvs[3]

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
            if fnmatch.fnmatch(file, '*.log'):
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
            
            addline = 1 if pipeline_class == "stranded" else 0
            print pipeline_class
            pull = "T" if len(open(path).readlines()) > 50 else "F"
            print len(open(path).readlines())
            loglen = len(open(path).readlines())
            fname=os.path.basename(path).split(".")[0].replace('fcounts_','').replace('_trim','')
            print fname
            total_reads=0
            successfully_assigned=0
            assigned_rate=0

            for i, line in enumerate(fp):
                if pull == "T":
                    if i == loglen-50+39-addline:
                        print line.split()
                        total_reads = line.split()[4]
                    elif i == loglen-50+40-addline:
                        successfully_assigned = line.split()[5]
                        assigned_rate = line.split()[6].lstrip('(').split('%')[0]
                        break
                else:
                    if i == 38 + addline:
                        total_reads = line.split()[4]
                    elif i == 39 + addline:
                        successfully_assigned = line.split()[5]
                        #assigned_rate = line.split()[6].lstrip('(').split('%')[0]
                        break
        assigned_rate = (float(successfully_assigned) / float(total_reads)) * 100    
        data = pd.DataFrame([[fname,total_reads,successfully_assigned,str(assigned_rate)]], \
                        columns=['name','total_reads','successfully_assigned','assigned_rate'])
        df = pd.concat([df,data],ignore_index=True)
    dataout(df, out_fname)

if __name__ == '__main__':
        main()

print u"....All Done. End of Script"
