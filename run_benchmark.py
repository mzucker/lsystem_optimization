#!/usr/bin/env python

import subprocess
import re

def main():
    
    test_cases = [
        ('sierpinski_arrowhead', 17, 12),
        ('sierpinski_triangle', 16, 11),
        ('dragon_curve', 26, 18),
        ('barnsley_fern', 13, 9),
        ('sticks', 16, 11),
        ('hilbert', 13, 9),
        ('pentaplexity', 9, 6)
    ]

    test_programs = [
        'lsystems_v0.py',
        'lsystems_v1.py',
        'lsystems_v2',
        'lsystems_v3',
        'lsystems_v4',
        'lsystems_v5',
        'lsystems_v6'
    ]

    test_programs = ['lsystems_v5', 'lsystems_v6']

    expected_python_slowdown = 100.

    expr = re.compile(r'generated (.*) segments in (.*) s \(')

    ostr = open('benchmark_results.txt', 'a')
    
    for prog in test_programs:

        first = True

        for name, max_depth_c, max_depth_py in test_cases:

            print(name, max_depth_c, max_depth_py)

            if prog.endswith('.py'):
                max_depth = max_depth_py
                args = [ 'python', prog ]
            else:
                max_depth = max_depth_c
                args = [ './' + prog ]

            args += [ '-x', '0', name, str(max_depth) ]

            if first:
                print('throwing away first run...')
                subprocess.run(args, capture_output=True, text=True)
                first = False

            print(args)
                
            result = subprocess.run(args, capture_output=True, text=True)

            output = result.stdout

            print(output)

            match = re.search(expr, output)

            assert match

            nseg = int(match.group(1))
            time = float(match.group(2))

            period = time / nseg
            ostr.write('{} {} {}\n'.format(prog, name, repr(period)))
            ostr.flush()

if __name__ == '__main__':
    main()
        
