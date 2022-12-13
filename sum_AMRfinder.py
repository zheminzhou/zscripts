# summarize amrfinder results. (amrfinder --plus -O <Organism> -o <output> -n <nucleotide file>)
import sys, os, numpy as np, pandas as pd, click, subprocess, multiprocessing


def read_amr(fname) :
    try :
        data = pd.read_csv(fname, sep='\t', header=0, na_filter=False)
    except :
        return None, None
    res = {}
    for d in data.values :
        if d[6].find('extended-spectrum') >=0 :
            d[11] += '/ESBL'
        g = d[8] if d[9] == d[8] else d[8] + ':' + d[9]
        for c in d[11].split('/') :
            key = '{0}/{1}'.format(g, c) if c != 'NA' else g
            if key not in res :
                res[key] = {d[5]:1}
            else :
                res[key][d[5]] = 1
    return fname, { key:','.join(sorted(r.keys())) for key, r in res.items() }


@click.command()
@click.option('-f', '--filelist', help='name for AMRfinder outputs')
def main(filelist) :
    pool = multiprocessing.Pool(10)
    with open(filelist, 'rt') as fin :
        fnames = [line.strip() for line in fin]
    results = {}
    fields = {}
    for fn, r in pool.imap_unordered(read_amr, fnames) :
        if fn :
            results[fn] = r
            for k, v in r.items() :
                fields[k] = 1
    fields = sorted(fields)
    print('#FNAME\t{0}'.format('\t'.join(fields)))
    for fn, res in sorted(results.items()) :
        print('{0}\t{1}'.format(fn, '\t'.join([ res.get(fld, '-') for fld in fields ])))


if __name__ == '__main__' :
    main()

