## Functions for classifying cells using a kNN classifier

#' Run PCA
#'
#' @param x Input data matrix, should be genes x samples
#' @param pcs.use Number of PCs to calculate
#'
#' @return List of PCA embeddings, loadings, and variance explained
#' @import Matrix
#' @import irlba
#'
fast_pca <- function(x, pcs.use) {
  x <- Matrix::t(x)

  cm <- Matrix::colMeans(x)
  pca.obj <- irlba::irlba(x, nv = pcs.use, center = cm)

  var.exp <- pca.obj$d
  pc.emb <- t(pca.obj$u); colnames(pc.emb) <- rownames(x);
  pc.load <- pca.obj$v; rownames(pc.load) <- colnames(x);
  names(var.exp) <- rownames(pc.emb) <- colnames(pc.load) <- paste("PC", 1:pcs.use)

  return(list(pc.load = pc.load, pc.emb = pc.emb))
}


#' Projects new data onto existing PCs
#'
#' @param newx Dataset to project, should be genes x samples
#' @param pc.load Existing pca loadings
#'
#' @return Projected PCA embeddings
#'
#' @import Matrix
#'
project_pca <- function(newx, pc.load) {
  if(!all(rownames(newx) == rownames(pc.load))) stop("New data must have same rownames")
  newx <- Matrix::t(newx)

  cm <- Matrix::colMeans(newx)
  pc.emb <- t((newx - cm) %*% pc.load)
  return(pc.emb)
}


#' Compute distance weighted voting scores for each reference class for each cell
#'
#' @param knn.res List of neighbors and distances, list output of FNN::get.knnx
#' @param ident Reference identity
#'
#' @return Reference classes by cells matrix of voting scores
#'
weighted_neighbor_voting <- function(knn.res, ident) {
  unique.ident <- unique(ident)

  sapply(1:nrow(knn.res$nn.index), function(i) {
    x <- knn.res$nn.index[i,]
    d <- 1/knn.res$nn.dist[i,]
    ident.tbl <- tapply(d, factor(ident[x]), sum)

    ident.counts <- vector(mode = "integer", length = length(unique.ident))
    names(ident.counts) <- unique.ident
    ident.counts[names(ident.tbl)] <- ident.tbl

    ident.counts <- ident.counts/sum(ident.counts)
    return(ident.counts)
  })
}


#' Pairwise correlation of the columns of two matrices
#'
#' @param data.1 First input matrix
#' @param data.2 Second input matrix
#' @param metric Metric to use (either pearson or cosine)
#'
#' @return Matrix of column to column correlations
#'
correlate_cols <- function(data.1, data.2, metric = "pearson") {
  cor.mat <- sapply(1:ncol(data.1), function(i) {
    sapply(1:ncol(data.2), function(j) {
      if (metric == "pearson") {
        cor(data.1[,i],data.2[,j])
      } else if (metric == "cosine") {
        lsa::cosine(data.1[,i],data.2[,j])
      }
    })
  })
  colnames(cor.mat) <- colnames(data.1)
  rownames(cor.mat) <- colnames(data.2)

  return(cor.mat)
}

correlate_cols <- compiler::cmpfun(correlate_cols)



#' Runs PCA on the training dataset
#'
#' @param norm.counts Normalized gene expression matrix, should be genes x samples
#' @param genes.use Genes to use for PCA
#' @param ident Cell type identity
#' @param pcs.use Number of PCs to calculate
#'
#' @return List of gene expression, cell identity, PCA embeddings, PCA loadings
#'
#' @import Matrix
#' @import irlba
#' @export
#'
TrainKNN <- function(norm.counts, genes.use, ident, pcs.use = 40) {
  norm.counts <- norm.counts[genes.use,]
  pca <- fast_pca(norm.counts, pcs.use = pcs.use)
  list(norm.counts = norm.counts, ident = ident,
       pc.emb = pca$pc.emb, pc.load = pca$pc.load)
}


#' Maps query clusters to reference clusters using correlation
#'
#' @param query.data Query gene expression matrix
#' @param query.clusters Factor of query clusters
#' @param train.data Reference gene expression matrix
#' @param train.clusters Factor of reference clusters
#' @param genes.use Genes to use for classification
#' @param metric Metric to use (either pearson or cosine)
#'
#' @return List of reference and query cluster means and cluster correlations
#'
#' @import compiler
#' @export
#'
MapClustersCor <- function(query.data, query.clusters, train.data, train.clusters,
                           genes.use = NULL, metric = "pearson") {
  stopifnot(metric %in% c("pearson", "cosine"))
  if (is.null(genes.use)) genes.use <- intersect(rownames(train.data), rownames(query.data))
  else genes.use <- genes.use[genes.use %in% rownames(query.data) & genes.use %in% rownames(train.data)]

  query.cluster.data <- t(apply(query.data[genes.use,], 1, function(x) tapply(x, query.clusters, mean)))
  ref.cluster.data <- t(apply(train.data[genes.use,], 1, function(x) tapply(x, train.clusters, mean)))
  cluster.correlations <- correlate_cols(query.cluster.data, ref.cluster.data, metric = metric)

  return(list(query.cluster.data = query.cluster.data,
              ref.cluster.data = ref.cluster.data,
              cluster.correlations = cluster.correlations))
}

MapClustersCor <- compiler::cmpfun(MapClustersCor)




#' Map cells to a reference dataset using a kNN classifier
#'
#' @param query.norm.counts Query gene expression matrix
#' @param train.knn Output of TrainKNN
#' @param genes.use Genes to use for classification
#' @param use.pca Whether or not to project to PCs first
#' @param k Number of neighbors
#'
#' @return Matrix of cell type probabilities
#'
#' @import FNN
#' @import compiler
#' @export
#'
MapKNN <- function(query.norm.counts, train.knn, genes.use = NULL, use.pca = T, k = 30) {
  if (is.null(genes.use)) genes.use <- rownames(query.norm.counts)

  if (use.pca) {
    genes.use <- intersect(genes.use, rownames(train.knn$pc.load))
    query.pc.emb <- project_pca(query.norm.counts[genes.use,], train.knn$pc.load)
    knn.res <- FNN::get.knnx(t(train.knn$pc.emb), t(query.pc.emb), k = k)
  } else {
    genes.use <- intersect(genes.use, rownames(train.knn$norm.counts))
    knn.res <- FNN::get.knnx(t(as.matrix(train.knn$norm.counts[genes.use,])),
                             t(as.matrix(query.norm.counts[genes.use,])), k = k)
  }

  ident <- as.character(train.knn$ident)
  cell.mapped.ident <- weighted_neighbor_voting(knn.res, ident)
  colnames(cell.mapped.ident) <- colnames(query.norm.counts)

  return(cell.mapped.ident)
}

MapKNN <- compiler::cmpfun(MapKNN)


#' Summarize cell classification by group
#'
#' @param cell.mapped.ident Mapped cell type probabilities
#' @param query.ident Either a list or a factor of cell type identities for the query dataset
#' @param ident.thresh Filter out cells that don't map to a reference identity
#' @param assign.cells Assigns each cell to a single reference group
#' @param min.ident.frac Filter out reference identities
#'
#' @return A reference groups x query groups matrix of probabilities
#' @export
#'
SummarizeKNN <- function(cell.mapped.ident, query.ident, ident.thresh = 0, assign.cells = T,
                         min.ident.frac = 0.05) {
  cell.mapped.ident <- cell.mapped.ident[,apply(cell.mapped.ident, 2, max) > ident.thresh]

  if(assign.cells) {
    cell.mapped.ident <- apply(cell.mapped.ident, 2, function(x) {
      y <- rep(0, length(x)); names(y) <- names(x);
      y[which.max(x)] <- 1
      y
    })
  }

  if (is.list(query.ident)) {
    query.ident <- lapply(query.ident, function(cells) intersect(cells, colnames(cell.mapped.ident)))
    query.mean.ident <- apply(cell.mapped.ident, 1, function(x) {
      sapply(query.ident, function(cells) mean(x[cells]))
    })
  } else {
    query.ident <- factor(query.ident)
    query.ident <- query.ident[colnames(cell.mapped.ident)]
    query.mean.ident <- apply(cell.mapped.ident, 1, function(x) tapply(x, query.ident, mean))
  }

  query.mean.ident <- t(query.mean.ident)
  return(query.mean.ident[apply(query.mean.ident, 1, max) >= min.ident.frac,])
}

SummarizeKNN <- compiler::cmpfun(SummarizeKNN)


#' Cross validation to find the optimal k
#'
#' @param train.knn List output of TrainKNN with gene expression matrix, PCs, and reference ident
#' @param k.range Range of k values to test
#' @param n.folds Number of folds to split the data
#' @param seed Random seed for splitting folds
#' @param n.cores Number of multithreading cores
#'
#' @return List of k value accuracies and the optimal k
#'
#' @import caret
#' @import snow
#' @export
#'
CrossValidateKNN <- function(train.knn, k.range, n.folds = 4, seed = NULL, n.cores = 4) {
  ident <- as.character(train.knn$ident); names(ident) <- names(train.knn$ident);
  cell.names <- names(ident)

  if (!is.null(seed)) set.seed(seed)
  test.folds <- caret::createFolds(ident, k = n.folds, list = T, returnTrain = F)

  test.folds <- lapply(test.folds, function(idx) cell.names[idx])
  train.folds <- lapply(test.folds, function(test.cells) cell.names[!cell.names %in% test.cells])

  knn.folds <- lapply(1:n.folds, function(i) {
    test.cells <- test.folds[[i]]
    train.cells <- train.folds[[i]]

    pc.train <- fast_pca(train.knn$norm.counts[rownames(train.knn$pc.load), train.cells],
                         pcs.use = nrow(train.knn$pc.emb))
    pc.train.emb <- pc.train$pc.emb
    pc.test.emb <- project_pca(train.knn$norm.counts[rownames(train.knn$pc.load), test.cells],
                               pc.train$pc.load)

    list(pc.train = pc.train.emb, pc.test = pc.test.emb)
  })

  cl <- snow::makeCluster(n.cores, type = "SOCK")
  snow::clusterExport(cl, c("knn.folds", "ident"), envir = environment())
  k.acc <- snow::parSapply(cl, k.range, function(k) {
    fold.acc <- sapply(knn.folds, function(fld) {
      train.ident <- ident[colnames(fld$pc.train)]
      test.ident <- ident[colnames(fld$pc.test)]

      knn.res <- FNN::get.knnx(t(fld$pc.train), t(fld$pc.test), k = k)
      classes <- apply(weighted_neighbor_voting(knn.res, train.ident), 2, function(x) names(x)[which.max(x)])
      sum(classes == test.ident)/length(classes)
    })
    mean(fold.acc)
  })
  snow::stopCluster(cl)

  best.idx <- which.max(k.acc)
  return(list(acc = k.acc, k = k.range[[best.idx]]))
}

CrossValidateKNN <- compiler::cmpfun(CrossValidateKNN)


#' Map mouse genes to human genes
#'
#' @param mouse.genes Vector of mouse genes
#'
#' @return vector of human genes with mouse genes as names
#'
#' @import compiler
#' @export
#'
MouseHumanMapping <- function(mouse.genes) {
  data(mouse.human.table)
  mouse.human.table <- mouse.human.table[!duplicated(mouse.human.table$Gene.name),]

  mouse2human <- mouse.human.table$Human.gene.name
  names(mouse2human) <- mouse.human.table$Gene.name

  mouse.genes <- mouse.genes[mouse.genes %in% names(mouse2human)]
  return(mouse2human[mouse.genes])
}

MouseHumanMapping <- compiler::cmpfun(MouseHumanMapping)
