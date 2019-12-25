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
            if fnmatch.fnmatch(file, '*hisat2.command.err'):
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
            fname=os.path.basename(path).split(".")[0].replace('.hisat2.command.err','')
            
            alined_rate=0
            if seqtype == 'SE':
                for i, line in enumerate(fp):
                    if line.find("reads; of these:") >= 0:
                        total_reads = line.split()[0]
                    elif line.find("aligned 0 times") >= 0:
                        unmapped_read = line.split()[0]
                    elif line.find("overall alignment rate") >= 0:
                        overall_alignment_rate = line.split()[0].replace('%','')
                        break
                percent_unmapped = float(unmapped_read)/float(total_reads) * 100
            elif seqtype == 'PE':
                for i, line in enumerate(fp):
                    if line.find("reads; of these:") >= 0:
                        total_reads = float(line.split()[0]) * 2
                    elif line.find("pairs aligned concordantly 0 times; of these:") >= 0:
                        unmapped_read = line.split()[0]
                    elif line.find("overall alignment rate") >= 0:
                        overall_alignment_rate = line.split()[0].replace('%','')
                        break
                percent_unmapped = float(unmapped_read)/total_reads * 100

            data = pd.DataFrame([[fname, unmapped_read, total_reads, percent_unmapped, overall_alignment_rate]], columns=['Sample_ID','unmapped_read','total_reads', 'percent_unmapped', 'overall_alignment_rate'])
    
        df = pd.concat([df,data],ignore_index=True)       
    print df.shape
    dataout(df, out_fname)

if __name__ == '__main__':
        main()

print u"....All Done. End of Script"
