import numpy as np
from pulse.square_wave import build_square_wave


def comparator(signal, num_samplesperbit, v_ref):
    if num_samplesperbit % 2 == 0:
        step = int(num_samplesperbit / 2)
    else:
        step = int((num_samplesperbit + 1) / 2 - 1)
    sample = []
    for i in np.arange(step, signal.shape[0], num_samplesperbit, dtype='int'):
        if signal[i] >= v_ref:
            sample.append(1)
        else:
            sample.append(0)

    return build_square_wave(num_samplesperbit, input_seq=sample)
