from sklearn.linear_model import Lasso
from skglm import WeightedLasso

import numpy as np
import math
import random


def defn_net(d, p, n, grp_sparsity = 0.5, sparsity = None, grp = None):
  if grp == None:
    grp = np.array(range(1, p + 1))
  edge = np.zeros((d * p * p))
  weight = np.array([1, 1, 1])
  signum = [1, -1]
  grpCt = len(set(grp))
  sparsity = max(min((n / (d * grpCt * p)), (0.05)), 0.01)
  num_non_zero = math.ceil(d * p * p * sparsity)

  # This piece of code is not checked
  if np.array_equal(grp, np.array(range(1, p + 1))) == False:
    grpCt_idx = list(range(1, grpCt + 1))
    grp_sel = random.sample(grpCt_idx, math.floor(grpCt * grp_sparsity))
    j_grp = np.where(np.isin(grp, grp_sel))[0]
    num_pergrp = math.floor(num_non_zero / len(j_grp))
    for j in j_grp:
      ind_seq = list(range((j-1) * p * d + 1, d * p * j + 1))
      ind_sel = random.sample(ind_seq, num_pergrp)
      edge_sel = np.array([random.choice(signum) for _ in range(num_pergrp)]) * 0.5
      edge[ind_sel] = edge_sel

  ind_seq = list(range(1, p * p * d + 1))
  ind_sel = random.sample(ind_seq, num_non_zero)
  edge_sel = np.array([random.choice(signum) for _ in range(num_non_zero)]) * 0.5
  edge[ind_sel] = edge_sel

  bottom_diag = np.diag(np.ones(p * (d - 1)))
  bottom_zero = np.zeros((p * (d - 1), p))
  bottom = np.column_stack((bottom_diag, bottom_zero))
  edge = edge.reshape(d, p, p)

  while True:
    A = np.concatenate(edge, axis=1)
    At = np.vstack((A, bottom))
    At_eigen = np.linalg.eigvals(At)
    if max(abs(At_eigen)) < 0.95: break
    edge = edge * 0.95

  return edge

def sim_data(n, edge, T, error_sd = 0.1, cutt = 30):
  d = edge.shape[0]
  p = edge.shape[1]
  x = np.zeros((n, p, T + (cutt * d)))
  for i in range(n):
    x[i, :, 0:d] = np.random.normal(0, 1, size=(p * d)).reshape((p, d))
    for j in range(d, x[0, :, :].shape[1]):
      x[i, :, j] = np.random.normal(0, error_sd, p)
      for l in range(d):
        x[i, :, j] = x[i, :, j] + edge[l, :, : ] @ x[i, :, j - l - 1]
  return x[:, :, cutt * d:]
