"""
Convert .h5 data from SLEAP to .csv for R
Following the example at https://sleap.ai/notebooks/Analysis_examples.html
load data using h5py package The shape of the locations matrix (the tracks dataset) 
(#frames, nodes, x/y, #ind) after it is transposed (with the .T). 
"""

# ------------------------------------------
# These are r codes to use python in Rstudio
library(reticulate)
use_python("C:/Users/nobua/.conda/envs/tmanual/python.exe", required = T)
repl_python()
# ------------------------------------------


import glob
import h5py
import numpy as np
from numpy.linalg import norm
from scipy.interpolate import interp1d
import pandas as pd


in_dir = "data/Cf-control-sleap/h5"
out_dir = "data/Cf-control-sleap/csv"
coord_name = ["x", "y"]
files = glob.glob(in_dir + "/*.h5")


def fill_missing(Y, kind="linear"):
    initial_shape = Y.shape
    Y = Y.reshape((initial_shape[0], -1))
    # Interpolate along each slice.
    for i in range(Y.shape[-1]):
        y = Y[:, i]
        # Build interpolant.
        x = np.flatnonzero(~np.isnan(y))
        f = interp1d(x, y[x], kind=kind, fill_value=np.nan, bounds_error=False)
        # Fill missing
        xq = np.flatnonzero(np.isnan(y))
        y[xq] = f(xq)
        # Fill leading or trailing NaNs with the nearest non-NaN values
        mask = np.isnan(y)
        y[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), y[~mask])
        Y[:, i] = y
        if sum(np.isnan(y)) > 0:
          print("error"+str(i))
          print("error"+(i))
    # Restore to initial shape.
    Y = Y.reshape(initial_shape)
    return Y


for f_name in files:
  with h5py.File(f_name, "r") as f:
    dset_names = list(f.keys())
    locations = f["tracks"][:].T
    node_names = [n.decode() for n in f["node_names"][:]]
  
  print(f_name)
  locations = fill_missing(locations)
  
  # identify female-male identity
  head_0_x = locations[:, 6, 0, 0]
  head_0_y = locations[:, 6, 1, 0]
  head_1_x = locations[:, 6, 0, 1]
  head_1_y = locations[:, 6, 1, 1]
  abd_0_x = locations[:, 10, 0, 0]
  abd_0_y = locations[:, 10, 1, 0]
  abd_1_x = locations[:, 10, 0, 1]
  abd_1_y = locations[:, 10, 1, 1]
  
  dis_h0_t1 = np.sqrt( np.square(head_0_x - abd_1_x) + np.square(head_0_y - abd_1_y) )
  dis_h1_t0 = np.sqrt( np.square(head_1_x - abd_0_x) + np.square(head_1_y - abd_0_y) )
  
  if sum(dis_h0_t1 < dis_h1_t0) < sum(dis_h1_t0 < dis_h0_t1):
    ind_name = ["f", "m"]
  else:
    ind_name = ["m", "f"]


  for i_ind in range(2):
    label_name = []
    coordinates = []
    treatment  = ind_name[i_ind]
    
    # headtip, pronotumfront, abdomentip
    for i_body in [6,7,10]: 
      for i_coord in range(locations.shape[2]):
        label = node_names[i_body] + "_" + coord_name[i_coord]
        loc_ndar = locations[:, i_body, i_coord, i_ind]
        loc_list = loc_ndar.tolist()

        coordinates.append(loc_list)
        label_name.append(label)

    df = pd.DataFrame(coordinates)
    df = df.T
    df = df.set_axis(label_name, axis='columns')
    out_name = f_name.replace(".h5", "_" + treatment + ".csv")
    out_name = out_name.replace(in_dir, out_dir)

    df.to_csv(out_name)

