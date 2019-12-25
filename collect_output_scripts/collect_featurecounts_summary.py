import sys
import os
print "python :" + sys.version
import ConfigParser
import pandas as pd
import fnmatch

argvs = sys.argv
argc = len(argvs)

if (argc < 3):
    print 'Usage: # python %s input_dir out_name ' % argvs[0]
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
            if fnmatch.fnmatch(file, '*.log'):
                yield os.path.join(root, file)

def error_cb(e):
    print 'Could not open "{0}" [{1}]: {2}'.format(e.filename, e.errno, e.strerror)
    raise OSError(e)


def main():
    filelist = []
    for gettxt in fild_all_files(inp_dir):
        print gettxt
        filelist.append(gettxt)
    
    print 'file num of dir : ' + str(len(filelist))
    
    df = pd.DataFrame() 
    for path in filelist:
        with open(path) as fp:

            fname=os.path.basename(path).split(".")[0].replace('fcounts_','').replace('_trim','')
            print fname
            total_reads=0
            successfully_assigned=0
            assigned_rate=0

            for i, line in enumerate(fp):
                if line.find("Total alignments") >= 0:
                    #print line.split()
                    total_reads = line.split()[4]
                elif line.find("Successfully assigned alignments") >= 0:
                    #print line.split()
                    successfully_assigned = line.split()[5]
                    break

        assigned_rate = (float(successfully_assigned) / float(total_reads)) * 100    
        data = pd.DataFrame([[fname,total_reads,successfully_assigned,str(assigned_rate)]], \
                        columns=['name','total_reads','successfully_assigned','assigned_rate'])
        df = pd.concat([df,data],ignore_index=True)
    dataout(df, out_fname)

if __name__ == '__main__':
        main()

print u"....All Done. End of Script"
