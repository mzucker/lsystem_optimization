######################################################################
#
# lsystems.py
#
# Matt Zucker
# Written for ENGR 52, Spring 2020
#
######################################################################
#
# Based on documentation in https://en.wikipedia.org/wiki/L-system and
# http://paulbourke.net/fractals/lsys/
#
# This version uses iteration and explicit string representation to
# build L-Systems. See lsystems_v2.py for an example using recursion
# and implicit representations.

import sys
import argparse
from datetime import datetime
from collections import namedtuple
import numpy as np

from plot_segments import plot_segments

LSystem = namedtuple('LSystem', 'start, rules, turn_angle_deg, draw_chars')

# dictionary mapping names to a few L-Systems we found on Wikipedia
KNOWN_LSYSTEMS = {
    
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
    )
    
}

# s: string
# rules: dictionary maps variables -> replacements
def lsystem_build_string(lsys, max_depth):

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

def lsystem_execute_symbol(lsys, symbol, cur_state, stack, segments):

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
    

# take a string and turn it into a set of line segments
# turn_angle is provided in degrees
# segments returned as an n-by-2-by-2 array where each segment is represented as
# [(x0, y0), (x1, y1)]
def lsystem_segments_from_string(lsys, lstring):

    cur_pos = np.array([0., 0.])
    cur_angle_deg = 0

    cur_state = ( cur_pos, cur_angle_deg )

    # stack of pos, angle pairs
    stack = []

    segments = []

    for symbol in lstring:
        cur_state = lsystem_execute_symbol(lsys, symbol, cur_state, stack, segments)
        
    return np.array(segments)

######################################################################

def lsystem_segments_recursive_helper(lsys, s,
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
            cur_state = lsystem_segments_recursive_helper(lsys, replacement, 
                                                          remaining_steps-1,
                                                          cur_state,
                                                          state_stack, segments)

        else: # execute symbol directly

            cur_state = lsystem_execute_symbol(lsys, symbol, cur_state,
                                               state_stack, segments)

    return cur_state

# just calls the helper function above
def lsystem_segments_recursive(lsys, max_depth):

    cur_pos = np.array([0., 0.])
    cur_angle_deg = 0

    cur_state = (cur_pos, cur_angle_deg)
    state_stack = []

    segments = []

    s = lsys.start

    lsystem_segments_recursive_helper(lsys, s,
                                      max_depth,
                                      cur_state,
                                      state_stack,
                                      segments)

    return np.array(segments)


######################################################################

def parse_options():

    parser = argparse.ArgumentParser(description='simple Python L-system renderer')

    parser.add_argument('-x', dest='max_segments', metavar='MAXSEGMENTS',
                        type=int, default=100000,
                        help='maximum number of segments to plot')

    parser.add_argument('lname', metavar='LSYSTEM', nargs=1,
                        help='name of desired L-system',
                        type=str,
                        choices=KNOWN_LSYSTEMS)

    parser.add_argument('max_depth', metavar='MAXDEPTH', nargs=1,
                        help='maximum depth to evaluate', type=int)

    parser.add_argument('-t', dest='text_only', action='store_true',
                        help='use text output instead of PNG')

    parser.add_argument('-s', dest='use_recursion', action='store_false',
                        default=False,
                        help='use string building method')

    parser.add_argument('-r', dest='use_recursion', action='store_true',
                        default=False,
                        help='use recursive method')
                        
    opts = parser.parse_args()

    opts.lsys = KNOWN_LSYSTEMS[opts.lname[0]]
    opts.max_depth = opts.max_depth[0]

    if opts.use_recursion:
        print('using recursive method')
    else:
        print('using string building method')

    return opts


# main function
def main():

    opts = parse_options()

    # time segment generation
    start = datetime.now()

    if opts.use_recursion:
        
        segments = lsystem_segments_recursive(opts.lsys, opts.max_depth)

    else:

        lstring = lsystem_build_string(opts.lsys, opts.max_depth)

        segments = lsystem_segments_from_string(opts.lsys, lstring)

    # print elapsed time
    elapsed = (datetime.now() - start).total_seconds()

    print('generated {} segments in {:.6f} seconds'.format(
        len(segments), elapsed))

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

    

