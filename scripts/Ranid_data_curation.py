#!/usr/bin/env python

# make sure you are in the directory that contains your file
# start the python interpreter from your terminal by simply running 'python'
# then you can run all of the following code!

import pandas as pd
import sys

# import file
dat = pd.read_table('RanidaeTissues-c95a4f5d3b9148fab16495517f44a505.txt')

# select taxon of interest
mask = (dat['specificepithet'] == 'catesbeianus') | (dat['specificepithet'] == 'catesbeiana')

# filter data for taxon of interest
filtered_dat = dat[mask]

# EXAMPLE FOR LATITUDE
## Get rid of missing latitude data
filtered_dat_cleaned = filtered_dat.dropna(subset=['decimallatitude'])
## Sort values
filtered_dat_sorted = filtered_dat_cleaned.sort_values(by='decimallatitude', ascending=True)
## Print data to see range
filtered_dat_sorted.decimallatitude

# EXAMPLE FOR LATITUDE
## Get rid of missing latitude data
filtered_dat_cleaned = filtered_dat.dropna(subset=['maximumelevationinmeters'])
## Sort values
filtered_dat_sorted = filtered_dat_cleaned.sort_values(by='maximumelevationinmeters', ascending=True)
## Print data to see range
filtered_dat_sorted.maximumelevationinmeters
