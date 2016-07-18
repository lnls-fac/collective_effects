#!/usr/bin/env python3

import optparse
import os
import pycolleff.process_wakes as ems

def main(pth2sv='analysis',silent = False):
    simul_data = ems.EMSimulData()
    ems.load_raw_data(simul_data,silent=False)
    ems.calc_impedance(simul_data,silent=False)
    ems.save_processed_data(simul_data,silent=False,pth2sv=pth2sv)
    return simul_data

if __name__ == '__main__':

    # configuration of the parser for the arguments
    parser = optparse.OptionParser()
    parser.add_option('-p','--noplot',dest='plot',action='store_true',
                      help="Show results", default=False)
    parser.add_option('-s','--silent',dest='silent',action='store_true',
              help="Print progress", default=False)
    parser.add_option('-c','--calc',dest='calc',action='store_true',
                      help="Calculate results", default=False)
    parser.add_option('--pth2sv',dest='pth2sv',type='str',
                      help="Path to save the data. Relative to the current folder",
                      default = 'analysis')
    (opts, _) = parser.parse_args()

    plot = not opts.plot
    silent = opts.silent
    pth2sv = opts.pth2sv
    path = os.path.abspath('.')
    fname  = os.path.sep.join([path,pth2sv,'SimulData.pickle'])

    if opts.calc:
        simul_data = main(pth2sv = pth2sv, silent=silent)
        salva = True
    else:
        simul_data = ems.load_processed_data(fname)
        salva = False

    ot = dict(save_figs=salva,pth2sv=pth2sv,show=False)
    ems.plot_wakes(simul_data,**ot)
    ems.plot_impedances(simul_data,**ot)
    ems.plot_losskick_factors(simul_data,**ot)
    ems.create_make_fig_file(path=pth2sv)
    if plot: ems.show_now()
