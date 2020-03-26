import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

def plot_segments(segments, image_filename='segment_plot.png'):

    assert len(segments.shape) == 3 and segments.shape[1:] == (2, 2)

    lc = LineCollection(segments, color='b')

    ax = plt.axes()

    ax.add_collection(lc)
    ax.autoscale()

    plt.axis('equal')
    plt.axis('off')

    plt.savefig('segment_plot.png')
    print('wrote segment_plot.png')


if __name__ == '__main__':

    segments = np.genfromtxt('segments.txt').reshape(-1, 2, 2)
    
    plot_segments(segments)

