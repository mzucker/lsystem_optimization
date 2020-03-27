#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import os

def main():

    istr = open('benchmark_results.txt', 'r')

    results = dict()

    for line in istr:

        line = line.rstrip()
        
        if not line:
            continue
        
        prog, test_name, period = line.split()
        period = float(period)

        if prog not in results:
            results[prog] = dict()

        assert test_name not in results[prog]
        
        results[prog][test_name] = period

    test_names = set( next(iter(results.values())).keys() )

    progs = list(sorted(results.keys()))

    nprogs = len(progs)

    versions = np.arange(nprogs)
    mean_times = np.zeros(nprogs)
    errs = np.zeros((2, nprogs))

    for pidx, prog in enumerate(progs):
        presults = results[prog]
        assert set(presults.keys()) == test_names
        result_data = np.array(list(presults.values())) * 1e6
        mean_times[pidx] = np.power(result_data.prod(), 1.0/len(result_data))
        errs[0, pidx] = mean_times[pidx] - result_data.min()
        errs[1, pidx] = result_data.max() - mean_times[pidx]

    assert np.all(errs >= 0)

    is_c = np.array([not p.endswith('.py') for p in progs])

    plt.errorbar(versions,
                 mean_times, errs, fmt='-', capsize=2,
                 color=[.7, .7, .7], ecolor='k',
                 linewidth=0.5, elinewidth=0.5, zorder=0,
                 barsabove=True)


    for i in range(2):
        mask = (is_c == bool(i))
        label = 'Python versions' if i == 0 else 'C versions'
        plt.plot(versions[mask], mean_times[mask], '.', zorder=1, label=label)

    plt.legend(loc='upper right')
    plt.xlabel('version')
    plt.ylabel('μs / segment')
    plt.title('L-System benchmark results')
        
    plt.yscale('log')
    #plt.show()

    plt.savefig('benchmark_results.png', dpi=300)

    os.system('mogrify -trim -bordercolor White -border 10x10 benchmark_results.png')
    


if __name__ == '__main__':
    main()
