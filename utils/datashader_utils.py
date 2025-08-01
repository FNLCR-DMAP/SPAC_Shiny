import pandas as pd
import datashader as ds
import datashader.transfer_functions as tf
from colorcet import fire

def scatter_heatmap(x, y, color=None, width=800, height=600):
    df = pd.DataFrame({'x': x, 'y': y})
    cvs = ds.Canvas(plot_width=width, plot_height=height)
    if color is not None:
        df['color'] = color
        agg = cvs.points(df, 'x', 'y', ds.mean('color'))
        img = tf.shade(agg, cmap=fire, how='eq_hist')
    else:
        agg = cvs.points(df, 'x', 'y', ds.count())
        img = tf.shade(agg, cmap=fire)
    return img.to_pil()