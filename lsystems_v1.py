#!/usr/bin/env python
######################################################################
#
# lsystems_v1.py
# Matt Zucker
#
######################################################################
#
# Based on documentation in https://en.wikipedia.org/wiki/L-system and
# http://paulbourke.net/fractals/lsys/

import sys
import argparse
from datetime import datetime
from collections import namedtuple
import numpy as np

from plot_segments import plot_segments

LSystem = namedtuple('LSystem', 'start, rules, turn_angle_deg, draw_chars')

# A few L-Systems found on pages linked above

KNOWN_LSYSTEMS = {

    'sierpinski_triangle': LSystem(
        start = 'F-G-G',
        rules = dict(F='F-G+F+G-F', G='GG'),
        turn_angle_deg = 120,
        draw_chars = None
    ),
    
    'sierpinski_arrowhead': LSystem(
        start = 'A',
        rules = dict(A='B-A-B', B='A+B+A'),
        turn_angle_deg = 60,
        draw_chars = None
    ),

    'dragon_curve': LSystem(
        start = 'FX',
        rules = dict(X='X+YF+', Y='-FX-Y'),
        turn_angle_deg = 90,
        draw_chars = None
    ),

    'barnsley_fern': LSystem(
        start = 'X',
        rules = dict(X='F+[[X]-X]-F[-FX]+X', F='FF'),
        turn_angle_deg = 25,
        draw_chars = None
    ),

    'sticks': LSystem(
        start = 'X',
        rules = dict(X='F[+X]F[-X]+X', F='FF'),
        turn_angle_deg = 20,
        draw_chars = 'F'
    ),

    'hilbert': LSystem(
        start = 'L',
        rules = dict(L='+RF-LFL-FR+', R='-LF+RFR+FL-'),
        turn_angle_deg = 90,
        draw_chars = 'F'
    ),

    'pentaplexity': LSystem(
        start = 'F++F++F++F++F',
        rules = dict(F='F++F++F+++++F-F++F'),
        turn_angle_deg = 36,
        draw_chars = None
    )
}

######################################################################
# make a big ol' string from an L-System starting from its start state
# using repeated string replacement.

def lsys_build_string(lsys, max_depth):

    lstring = lsys.start

    rules = lsys.rules

    for i in range(max_depth):

        output = ''

        for symbol in lstring:
            if symbol in rules:
                output += rules[symbol]
            else:
                output += symbol

        lstring = output

    return lstring

######################################################################
# "draw" a single symbol from an L-System string
#
# stack is read-write
# segments is for appending
# cur_state is read and returned

def _lsys_execute_symbol(lsys, symbol, cur_state, stack, segments):

    cur_pos, cur_angle_deg = cur_state

    if symbol.isalpha():
        
        if lsys.draw_chars is None or symbol in lsys.draw_chars:
            
            cur_theta = cur_angle_deg * np.pi / 180
            offset = np.array([np.cos(cur_theta), np.sin(cur_theta)])
            new_pos = cur_pos + offset
            segments.append([cur_pos, new_pos])
            cur_pos = new_pos
            
    elif symbol == '+':
        
        cur_angle_deg += lsys.turn_angle_deg
        
    elif symbol == '-':
        
        cur_angle_deg -= lsys.turn_angle_deg
        
    elif symbol == '[':
        
        stack.append( ( cur_pos, cur_angle_deg ) )
        
    elif symbol == ']':
        
        return stack.pop()
        
    else:
        
        raise RuntimeError('invalid symbol:' + symbol)

    return cur_pos, cur_angle_deg
    

######################################################################
# take a string and turn it into a set of line segments, returned as
# an n-by-2-by-2 array where each segment is represented as
#
#  [(x0, y0), (x1, y1)]

def lsys_segments_from_string(lsys, lstring):

    cur_pos = np.array([0., 0.])
    cur_angle_deg = 0

    cur_state = ( cur_pos, cur_angle_deg )

    # stack of pos, angle pairs
    stack = []

    segments = []

    for symbol in lstring:
        cur_state = _lsys_execute_symbol(lsys, symbol, cur_state,
                                         stack, segments)
        
    return np.array(segments)

######################################################################
# build up segment list recursively - don't call this function
# directly, use lsys_segments_recursive instead

def _lsys_segments_r(lsys, s,
                     remaining_steps,
                     cur_state,
                     state_stack,
                     segments):

    # for each symbol in input
    for symbol in s:

        # see if we can run a replacement rule for this symbol
        if remaining_steps > 0 and symbol in lsys.rules:

            # get the replacement
            replacement = lsys.rules[symbol]
            
            # recursively call this function with fewer remaining steps
            cur_state = _lsys_segments_r(lsys, replacement, 
                                         remaining_steps-1,
                                         cur_state,
                                         state_stack, segments)

        else: # execute symbol directly

            cur_state = _lsys_execute_symbol(lsys, symbol, cur_state,
                                             state_stack, segments)

    return cur_state

######################################################################
# build up segment list recursively by calling helper function above

def lsys_segments_recursive(lsys, max_depth):

    cur_pos = np.array([0., 0.])
    cur_angle_deg = 0

    cur_state = (cur_pos, cur_angle_deg)
    state_stack = []

    segments = []

    s = lsys.start

    _lsys_segments_r(lsys, s,
                     max_depth,
                     cur_state,
                     state_stack,
                     segments)

    return np.array(segments)

######################################################################
# parse command-line options for this program

def parse_options():

    parser = argparse.ArgumentParser(
        description='simple Python L-system renderer')

    parser.add_argument('lname', metavar='LSYSTEM', nargs=1,
                        help='name of desired L-system',
                        type=str,
                        choices=KNOWN_LSYSTEMS)

    parser.add_argument('max_depth', metavar='MAXDEPTH', nargs=1,
                        help='maximum depth to evaluate', type=int)

    parser.add_argument('-x', dest='max_segments', metavar='MAXSEGMENTS',
                        type=int, default=100000,
                        help='maximum number of segments to plot')
    
    parser.add_argument('-t', dest='text_only', action='store_true',
                        help='use text output instead of PNG')

    parser.add_argument('-s', dest='use_recursion', action='store_false',
                        default=False,
                        help='use string building method (default)')

    parser.add_argument('-r', dest='use_recursion', action='store_true',
                        default=False,
                        help='use recursive method')
                        
    opts = parser.parse_args()

    opts.lname = opts.lname[0]
    opts.max_depth = opts.max_depth[0]

    opts.lsys = KNOWN_LSYSTEMS[opts.lname]

    if opts.use_recursion:
        print('using recursive method')
    else:
        print('using string building method')

    return opts

######################################################################
# main function

def main():

    opts = parse_options()

    # time segment generation
    start = datetime.now()

    if opts.use_recursion:
        
        segments = lsys_segments_recursive(opts.lsys, opts.max_depth)

    else:

        lstring = lsys_build_string(opts.lsys, opts.max_depth)

        segments = lsys_segments_from_string(opts.lsys, lstring)

    # print elapsed time
    elapsed = (datetime.now() - start).total_seconds()

    print('generated {} segments in {:.6f} s ({:.3f} us/segment)'.format(
        len(segments), elapsed, 1e6 * elapsed/len(segments)))

    if opts.max_segments >= 0 and len(segments) > opts.max_segments:
        print('...maximum of {} segments exceeded, skipping output!'.format(
            opts.max_segments))
        return

    if opts.text_only:
        np.savetxt('segments.txt', segments.reshape(-1, 4))
    else:
        plot_segments(segments)

if __name__ == '__main__':
    main()

    

