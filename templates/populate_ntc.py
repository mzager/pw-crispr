#!/usr/bin/env python3

import os
import pandas as pd

# Counts table
counts_fp = "input.counts.txt"
assert os.path.exists(counts_fp)
counts = pd.read_csv(counts_fp, sep='\t')

# Library definition file
library_fp = "${library}"
assert os.path.exists(library_fp)
library = pd.read_csv(library_fp, sep=',')

# Iterate over each of the guides in the library which share a sequence
for guide_seq, shared_guides in library.groupby('sgRNA'):

    # If there is more than one
    if shared_guides.shape[0] > 1:

        # Get the ID of the guide which does have valid counts
        guide_counts = counts.set_index('sgRNA').reindex(index=shared_guides.guide.values).dropna()

        # There should be no more than one set of counts
        msg = f"Found multiple rows of counts from mageck for the identical guide {guide_seq}"
        assert guide_counts.shape[0] <= 1, shared_guides

        # If there are no counts
        if guide_counts.shape[0] == 0:

            # Skip it
            continue

        # Drop the guide and gene columns from the counts, convert to a dict
        guide_counts = guide_counts.drop(columns=['Gene']).iloc[0].to_dict()
        
        # Add rows to the counts table for each of the guides which were omitted by mageck count
        counts = pd.concat(
            [
                counts,
                pd.DataFrame([
                    {
                        'sgRNA': r.guide,
                        'Gene': r.Gene,
                        **guide_counts
                    }
                    for _, r in shared_guides.iterrows()
                    if r.guide not in counts.sgRNA.values
                ])
            ]
        )

# Make sure that all counts are numeric
counts = counts.astype(
    {
        cname: 'int'
        for cname in counts.columns.values
        # Skipping the columns with gene and guide names
        if cname not in ['Gene', 'sgRNA']
    }
)

# Sort by gene and guide
counts = counts.sort_values(by=['Gene', 'sgRNA'])

# Write out the counts table
counts.to_csv('counts.txt', sep='\t', index=None)
