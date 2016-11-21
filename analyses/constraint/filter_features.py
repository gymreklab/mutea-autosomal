#!/usr/bin/env python

import sys
import pandas as pd

try:
    ffile = sys.argv[1]
except:
    sys.stderr.write("Usage: ./filter_features.py <featurefile>\n")
    sys.exit(1)

# Load data
features = pd.read_csv(ffile, sep="\t")
features["filter"] = False
features.ix[(features["length"]>=80), "filter"] = True
features.ix[np.isnan(features["reptiming"]), "reptiming"] = 0 # set to average

# Motifs
feature_counts = features.groupby("motif", as_index=False).agg({"filter": len})
feature_counts.columns = ["motif","count"]
motifs = list(feature_counts[feature_counts["count"]>=100]["motif"])
features.ix[features["motif"].apply(lambda x: x not in motifs), "filter"] = True

# Output
features.to_csv(sys.stdout, sep="\t", index=False)
