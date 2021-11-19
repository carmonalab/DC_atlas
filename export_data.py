#conda install scanpy

import scanpy as sc
import anndata
import pandas as pd
import numpy as np
from scipy import io
import glob
import re
import os

ddir = "/Users/mass/Dropbox/CSI/Datasets/Gerhard_JEM"
#print(glob.glob("%s*.h5ad" % ddir))

files = glob.glob("/Users/mass/Dropbox/CSI/Datasets/Gerhard_JEM/*.h5ad")


for f in files:

	root = re.sub('.*adata_', '', f)
	root = re.sub('.h5ad', '', root)

	path = os.path.join(ddir, "exported", root)
	if not os.path.exists(path):
		os.mkdir(path)

		# Read scanpy .h5ad object
		adata = sc.read_h5ad(f)

		# Write count matrix 
		io.mmwrite("%s/matrix.mtx" % path, adata.X)

		# Write meta annotation dataframe to .csv
		adata.obs.to_csv("%s/metadata.csv" % path)

		# Write gene names to .tsv
		adata.var_names.to_frame().to_csv("%s/genes.tsv" % path, sep= '\t')	

		# Write barcodes to .tsv 
		adata.obs.index.to_frame().to_csv("%s/barcodes.tsv" % path, sep= '\t')

print("Done!")