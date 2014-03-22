import numpy as np


def get_dims(FACTOR=0.45, WIDTH=360):
    WIDTH =350.0  # the number latex spits out
    FACTOR = 0.45  # the fraction of the width you'd like the figure to occupy
    fig_width_pt  = WIDTH * FACTOR

    inches_per_pt = 1.0 / 72.27
    golden_ratio  = (np.sqrt(5) - 1.0) / 2.0  # because it looks good

    fig_width_in  = fig_width_pt * inches_per_pt  # figure width in inches
    fig_height_in = fig_width_in * golden_ratio   # figure height in inches
    return [fig_width_in, fig_height_in] # fig dims as a list

