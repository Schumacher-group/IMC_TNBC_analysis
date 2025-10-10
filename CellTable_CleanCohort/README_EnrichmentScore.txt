Enrichment score per ROI for:

pre treatment samples (Responders vs Non-Responders)
Non-Responder samples (pre vs post treatment samples)

pkl files: a serialized pickle file.

To open:

import pickle

with open('neighbours_matrix.pkl', 'rb') as f:
    neigh_all = pickle.load(f)

ROIs are all together, i.e. in pre-treatment samples, the file contains both Responders and Non-Responders. They should be split when analysing.