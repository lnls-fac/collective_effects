#!/usr/bin/env python3

import argparse
import subprocess as sp
import psutil as psu

def main(fname, modes=(0,1), nr_th=1):

    # Find the line to change in the reference file
    lines, ind = [], -1
    with open(fname,'r') as f:
        lines = f.readlines()
    for i,line in enumerate(lines):
        if line.strip().lower().startswith('modes='):
            ind = i
    if ind<0:
        raise Exception('base file corrupted.')

    # Write the input files and submit calculations:
    procs = []
    for i in range(nr_th):
        filename = '{0:s}_{1:02d}.txt'.format(fname.partition('.')[0],i)
        str_modes = ' '.join(['{0:d}'.format(x) for x in modes[i::nr_th]])
        lines[ind] = 'Modes= {0:s}\n'.format(str_modes)
        with open(filename,'w') as f:
            f.writelines(lines)

        procs.append(psu.Popen(['echo2d', filename], stdout=sp.PIPE))

    while procs:
        print()
        status = False
        for i, p in enumerate(procs):
            line = p.stdout.readline().decode()
            print('Thread {0:2d} : {1:s}'.format(i,line),end='')
            if p.poll() is not None: procs.remove(p)
        print(30*'-')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('fileName',action='store', default='input_in.txt',
        help='Name of the input file.',)
    parser.add_argument('-m','--modes',action='store',nargs='+',
        help='Modes to calculate the wake.', type=str)
    parser.add_argument('-t','--threads',action='store',
        help='Number of threads to use.', type=int)

    args = parser.parse_args()

    modes = []
    for m in args.modes:
        ms = m.split(':')
        if len(ms) == 1:
            modes.append(int(ms[0]))
        elif len(ms) == 2:
            modes.extend(  list( range(int(ms[0]), int(ms[1])) )  )
        elif len(ms) == 3:
            modes.extend(  list( range(int(ms[0]), int(ms[1]), int(ms[2])) )  )
        else:
            raise Exception('Pattern not recognized')
    modes = sorted(modes)
    modes = [x for i,x in enumerate(modes) if ( i==0 or (i>0 and x>modes[i-1]) )]

    print('Calulating Wake using echo2d for modes:')
    print(' '.join(['{0:2d}'.format(x) for x in modes]))
    print('Using {0:2d} processors'.format(args.threads))

    main(args.fileName,modes=modes,nr_th=args.threads)
