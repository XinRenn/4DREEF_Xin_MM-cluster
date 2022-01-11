# Xin_MM-cluster
to cluster multimodels results

1.use monthly data to calculate the distance between models
	(script: cluster-var.ncl  output: cluster_var.nc)

2.use heritation clustering method and distance file(cluster_var.nc) to cluster and draw the tree figure.
	(script: cluster-PCC-avglink.ncl)

3.draw clustered contour figures
	(script: contour_var_cluster.ncl)
