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

join_col = 'Geneid'
out_col = '_fcount'

def makeDir(dname):
    if os.path.exists(dname) is False:
        os.mkdir(dname)
        print '%s (dir) created.' % dname
    else:
        print '%s (dir) is already exists.' % dname

def makeDataframe(fname):
    dataf = pd.read_table(fname, delimiter='\t', skiprows=1)
    pref = os.path.basename(fname).split(".")[0].replace("fcounts_","").replace("_trim","") + out_col
    dataf.rename(columns=(lambda x: pref if '_trim' in x else x), inplace=True)
    return dataf

def dataout(data, filename):
    data.to_csv(filename + '.txt', sep="\t")

def fild_all_files(directory):
    print directory
    for root, dirs, files in os.walk(directory, onerror=error_cb):
        for file in files:
            if fnmatch.fnmatch(file, '*_trim.txt'):
                yield os.path.join(root, file)

def error_cb(e):
    print 'Could not open "{0}" [{1}]: {2}'.format(e.filename, e.errno, e.strerror)
    raise OSError(e)

def merge_countdata(df_list, outcol):
    m_df = df_list[0]
    df_list.pop(0)
    for df in df_list:
        cols = [x for x in df.columns if join_col in x or x.endswith(str(outcol))]
        m_df = pd.merge(m_df, df[cols], on=join_col, how='outer', suffixes=['',''])

    #if inp_toint == "TRUE":
    #    m_df.ix[:,1:] = m_df.ix[:,1:].astype(int)   
    #print m_df.columns
    #m_df = m_df.reindex_axis(sorted(m_df.columns), axis=1)
    m_df.rename(columns=lambda x: x.replace(str(outcol), ''), inplace=True)
    return m_df

def main():
    filelist = []
    for gettxt in fild_all_files(inp_dir):
        print gettxt
        filelist.append(gettxt)
    
    print 'file num of dir : ' + str(len(filelist))
    
    df_list = [makeDataframe(fname) for fname in filelist]
    print(str(len(df_list)) + ' files convert to dataframe.')
    #print df_list[0].head(10)
    #print df_list[2].head(10)
    
    out_count = merge_countdata(df_list, out_col)
    #print out_count.head(10)
    print out_count.shape
    dataout(out_count, out_fname)

if __name__ == '__main__':
        main()

print u"....All Done. End of Script"
