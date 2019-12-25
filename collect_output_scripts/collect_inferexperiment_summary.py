import sys
import os
print "python :" + sys.version
import ConfigParser
import pandas as pd
import fnmatch

argvs = sys.argv
argc = len(argvs)

if (argc < 3):
    print 'Usage: # python %s input_dir out_name seqtype' % argvs[0]
    quit()

inp_dir = argvs[1]
out_fname = argvs[2]
seqtype = argvs[3]

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
            if fnmatch.fnmatch(file, '*.inferexp.txt'):
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
            fname=os.path.basename(path).split(".")[0]

            if seqtype == 'SE':
 
                for i, line in enumerate(fp):
                    if line.find("failed to determine") >= 0:
                        undetermined_fraction = line.split()[6]
                    elif line.find("++,--") >= 0:
                        sense_fraction = line.split()[6]
                    elif line.find("+-,-+") >= 0:
                        antisense_fraction = line.split()[6]
                        break

            elif seqtype == 'PE':

                for i, line in enumerate(fp):
                    if line.find("failed to determine") >= 0:
                        undetermined_fraction = line.split()[6]
                    elif line.find("1++,1--,2+-,2-+") >= 0:
                        antisense_fraction = line.split()[6]
                    elif line.find("1+-,1-+,2++,2--") >= 0:
                        sense_fraction = line.split()[6]
                        break

            data = pd.DataFrame([[fname, undetermined_fraction, sense_fraction, antisense_fraction]], columns=['name','undetermined_fraction','sense_fraction', 'antisense_fraction'])
    
        df = pd.concat([df,data],ignore_index=True)       
    print df.shape
    dataout(df, out_fname)

if __name__ == '__main__':
        main()

print u"....All Done. End of Script"
