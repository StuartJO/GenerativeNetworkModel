function szvals = scaled_sigmoid_edge_vals(vals)

kk = 1./(1+exp(-zscore(vals)));
szvals = (kk-min(kk))/(max(kk)-min(kk));