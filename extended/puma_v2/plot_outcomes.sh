#!/bin/sh
plot_outcomes.py --matched_cats=mrc,vlssr,sumss,nvss \
	--pref_cats=nvss,sumss,mrc,vlssr \
	--input_bayes=/home/jline/Documents/cataloguing/RTS/mrc_extras/puma_v2/bayes_m-v-s-n_v2.txt \
	--cat_freqs=408,74,843,1400 \
	--prob_thresh=0.8,0.95 \
	--epsilon_thresh=0.1 --chi_thresh=10 \
	--resolution=00:03:00 \
	--split=00:01:15 --query=2221-023 ##--accept --num_matches=8
