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
  pc.load <- pca.obj$v; rownames(pc.load) <- genes.use;
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
  pca <- fast_pca(norm.counts[genes.use,], pcs.use = pcs.use)
  list(norm.counts = norm.counts[genes.use,], ident = ident,
       pc.emb = pca$pc.emb, pc.load = pca$pc.load)
}


#' Map cells to a reference dataset using a kNN classifier
#'
#' @param query.norm.counts Query gene expression matrix
#' @param train.knn Output of TrainKNN
#' @param genes.use Genes to use for classification
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
  unique.ident <- unique(ident)

  cell.mapped.ident <- apply(knn.res$nn.index, 1, function(x) {
    ident.counts <- vector(mode = "integer", length = length(unique.ident))
    names(ident.counts) <- unique.ident

    ident.tbl <- table(ident[x])
    ident.counts[names(ident.tbl)] <- ident.tbl
    ident.counts/k
  })
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
