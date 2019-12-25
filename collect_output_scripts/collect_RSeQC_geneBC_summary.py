import sys
import os
print "python :" + sys.version
import ConfigParser
import pandas as pd
import fnmatch
from numpy import std,mean
import operator

argvs = sys.argv
argc = len(argvs)

if (argc < 3):
    print 'Usage: # python %s input_dir strand_option out_name' % argvs[0]
    quit()

inp_dir = argvs[1]
strand_option = argvs[2]
out_fname = argvs[3]

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
            if fnmatch.fnmatch(file, '*.' + strand_option + '.geneBodyCoverage.txt'):
                yield os.path.join(root, file)

def error_cb(e):
    print 'Could not open "{0}" [{1}]: {2}'.format(e.filename, e.errno, e.strerror)
    raise OSError(e)

def pearson_moment_coefficient(lst):
    '''measure skewness'''
    mid_value = lst[int(len(lst)/2)]
    sigma = std(lst, ddof=1)
    tmp = []
    for i in lst:
        tmp.append(((i - mid_value)/sigma)**3)
    return mean(tmp)

def main():
    filelist = []
    for gettxt in fild_all_files(inp_dir):
        print gettxt
        filelist.append(gettxt)
    
    print 'file num of dir : ' + str(len(filelist))
    
    column_name = range(1, 101)
    column_name.insert(0, 'name')
    df = pd.DataFrame(columns=column_name) 
    index = 0
    for path in filelist:
        with open(path) as fp:
            fname=os.path.basename(path).split(".")[0]
            dataset=[]
            for i, line in enumerate(fp):
                line = line.strip()
                #print line
                if line.startswith("Percentile"):
                    continue
                f = line.split()
                name =  fname
                dat = [float(i) for i in f[1:]]
                skewness = pearson_moment_coefficient(dat)
                dataset.append((name, [(i -min(dat))/(max(dat) - min(dat)) for i in dat], skewness))
                dataset.sort(key = operator.itemgetter(2), reverse=True)
                data = dataset[0][1]
                data.insert(0, dataset[0][0])
                df.loc[index] = data
        index = index + 1

    print df.shape
    dataout(df, out_fname)

if __name__ == '__main__':
        main()

print u"....All Done. End of Script"
