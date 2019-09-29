import sys
import os
print "python :" + sys.version
import ConfigParser
import pandas as pd
import fnmatch

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
    for root, dirs, files in os.walk(directory):
        for file in files:
            if fnmatch.fnmatch(file, '*.' + strand_option + '_readdist.txt'):
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
            fname=os.path.basename(path).split(".")[0]
            totalread=0
            totaltag=0
            assignedtag=0
            cds=0
            utr5=0
            utr3=0
            intron=0
            tssup1=0
            tssup5=0
            tssup10=0
            tssdown1=0
            tssdown5=0
            tssdown10=0
    
            for i, line in enumerate(fp):
                if i == 0:
                    totalread = line.split()[2]
                elif i == 1:
                    totaltag = line.split()[2]
                elif i == 2:
                    assignedtag = line.split()[3]
                elif i == 5:
                    cds = line.split()[2]
                elif i == 6:
                    utr5 = line.split()[2]
                elif i == 7:
                    utr3 = line.split()[2]
                elif i == 8:
                    intron = line.split()[2]
                elif i == 9:
                    tssup1 = line.split()[2]
                elif i == 10:
                    tssup5 = line.split()[2]
                elif i == 11:
                    tssup10 = line.split()[2]
                elif i == 12:
                    tesdown1 = line.split()[2]
                elif i == 13:
                    tesdown5 = line.split()[2]
                elif i == 14:
                    tesdown10 = line.split()[2]
                    break
            data = pd.DataFrame([[fname,totalread,totaltag,assignedtag,cds,utr5,utr3,intron,tssup1,tssup5,tssup10,tesdown1,tesdown5,tesdown10]], columns=['name','totalread','totaltag','assignedtag','cds','utr5','utr3','intron','tssup1','tssup5','tssup10','tesdown1','tesdown5','tesdown10'])
        df = pd.concat([df,data],ignore_index=True)       
    print df.shape
    dataout(df, out_fname)

if __name__ == '__main__':
        main()

print u"....All Done. End of Script"
