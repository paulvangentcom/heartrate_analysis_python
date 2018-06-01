import numpy as np
import sys
sys.path.append('../')
import heartbeat as hb
import matplotlib.pyplot as plt

def mark_clipping(data):
    '''function that marks start and end of clipping part
    
    possibly with start/end of slope?
    '''
    pass


if __name__ == '__main__':
    data = hb.get_data('clipping_part.csv')
    plt.plot(data)
    plt.show()