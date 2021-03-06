#!/usr/bin/env python

from Builder import get_workspace
import argparse

parser = argparse.ArgumentParser(description='Build binned workspaces.')
parser.add_argument('argument', type=str, choices = ['bins','chans','samps','nps','events'],
                    help='The parameter to be scaled')
parser.add_argument('-l',dest='somerange', type=int, nargs='+',
                    help='A list of options to scale over')
parser.add_argument('-r',dest='somerange', type=int, nargs='+',
                    help='The range to scale over')

arg = parser.parse_args().argument
somerange = parser.parse_args().somerange
if len(somerange) <= 2:
    if len(somerange) > 1:
        somerange = range(somerange[0],somerange[1])
    else:
        somerange = range(somerange[0])        
print somerange

d = {'events':1000, 
     'chans':1,
     'samps':1,
     'nps':1,
     'bins':10}

for i in somerange:
    print "writing {} {}".format(i, arg)
    d[arg] = i
    workspace = get_workspace(nchannels = d['chans'], nsamples = d['samps'], events = d['events'], nbins = d['bins'], nnps = d['nps'])
    workspace.SetName('BinnedWorkspace')
    workspace.writeToFile("output/workspace{}channels{}samples{}events{}bins{}nps.root".format(d['chans'], d['samps'], d['events'], d['bins'], d['nps']))
