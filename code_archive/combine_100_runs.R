# to combine the 100 runs together
for (parameters in c('networks_FALSE_kmeans_silhouette_euclidean_',
                     'networks_TRUE_kmeans_silhouette_euclidean_',
                     'random_FALSE_kmeans_silhouette_euclidean_',
                     'random_TRUE_kmeans_silhouette_euclidean_',
                     'random_FALSE_k3_',
                     'random_TRUE_k3_')){
  files <- list.files(pattern = parameters)
  mnem <- combine.mnem (parameters,files)
  print(paste(parameters,'done'))
}