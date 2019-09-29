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
out_fname = argvs[2] + '_raw'
out_fname_norm = argvs[2] + '_TPM'

join_col = 'Name'
out_col = '_NumReads'
out_col_norm = '_TPM'

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
            if fnmatch.fnmatch(file, 'quant.sf'):
                yield os.path.join(root, file)

def makeDataframe(fname):
    dataf = pd.read_table(fname, delimiter='\t')
    samplename = os.path.split(os.path.dirname(fname))[1].replace('sailfish_','').replace('_trim','')
    pref = samplename + out_col
    pref_norm = samplename + out_col_norm

    dataf.rename(columns=(lambda x: pref if 'NumReads' in x else x), inplace=True)
    dataf.rename(columns=(lambda x: pref_norm if 'TPM' in x else x), inplace=True)
    return dataf

def merge_countdata(dset, outcol):
    df_list = list(dset)
    m_df = df_list[0]
    df_list.pop(0)
    for df in df_list:
        cols = [x for x in df.columns if join_col in x or x.endswith(str(outcol))]
        m_df = pd.merge(m_df, df[cols], on=join_col, how='outer', suffixes=['',''])

    trimcols = [x for x in m_df.columns if join_col in x or x.endswith(str(outcol))]
    m_df = m_df[trimcols]
    m_df.rename(columns=lambda x: x.replace(str(outcol), ''), inplace=True)
    return m_df

def main():
    filelist = []
    for gettxt in fild_all_files(inp_dir):
        #print gettxt
        filelist.append(gettxt)

    print 'file num of dir : ' + str(len(filelist))

    dset = [makeDataframe(fname) for fname in filelist]
    print(str(len(dset)) + ' files convert to dataframe.')
    print dset[0].head(10)

    out_count = merge_countdata(dset, out_col)
    print out_count.head(10)
    print out_count.shape

    out_count_norm = merge_countdata(dset, out_col_norm)
    print out_count_norm.head(10)
    print out_count_norm.shape

    print list(set(out_count.columns) - set(out_count_norm.columns))

    dataout(out_count, out_fname)
    dataout(out_count_norm, out_fname_norm)


if __name__ == '__main__':
        main()

print u"....All Done. End of Script"
