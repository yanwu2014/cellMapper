## Functions for classifying cells using a kNN classifier

#' Run PCA
#'
#' @param norm.data Input data matrix, should be genes x samples
#' @param genes.use Subset of genes to use
#' @param pcs.use Number of PCs to calculate
#'
#' @return List of PCA embeddings, loadings, and variance explained
#' @import Matrix
#' @import irlba
#' @export
#'
FastPCA <- function(norm.data, genes.use, pcs.use = 40) {
  genes.use <- intersect(genes.use, rownames(norm.data))
  norm.data <- t(norm.data[genes.use,])

  cm <- Matrix::colMeans(norm.data)
  pca.obj <- irlba::irlba(norm.data, nv = pcs.use, center = cm)

  var.exp <- pca.obj$d
  pc.emb <- t(pca.obj$u); colnames(pc.emb) <- rownames(norm.data);
  pc.load <- pca.obj$v; rownames(pc.load) <- genes.use;
  names(var.exp) <- rownames(pc.emb) <- colnames(pc.load) <- paste("PC", 1:pcs.use)

  return(list(var.exp = var.exp, pc.load = pc.load, pc.emb = pc.emb))
}


#' Projects new data onto existing PCs
#'
#' @param new.data Dataset to project, should be genes x samples
#' @param pc.load Existing pca loadings
#' @param genes.use Genes to use for projection
#'
#' @return Projected PCA embeddings
#'
#' @import Matrix
#' @export
#'
ProjectPCA <- function(new.data, pc.load, genes.use) {
  new.data <- t(new.data[genes.use,])
  pc.load <- pc.load[genes.use,]

  cm <- Matrix::colMeans(new.data)
  pc.emb <- t((new.data - cm) %*% pc.load)
  return(pc.emb)
}



#' Map cells to a reference dataset using a kNN classifier
#'
#' @param query.data Query matrix
#' @param ref.data Reference matrix
#' @param ref.ident Reference cell types
#' @param k Number of neighbors
#'
#' @return Matrix of cell type probabilities
#'
#' @import FNN
#' @import compiler
#' @export
#'
MapCellsKNN <- function(query.data, ref.data, ref.ident, k = 30) {
  knn.res <- FNN::get.knnx(t(ref.data), t(query.data), k = k)

  ref.ident <- as.character(ref.ident)
  unique.ident <- unique(ref.ident)

  cell.ident.frac <- apply(knn.res$nn.index, 1, function(x) {
    ident.counts <- vector(mode = "integer", length = length(unique.ident))
    names(ident.counts) <- unique.ident

    ident.tbl <- table(ref.ident[x])
    ident.counts[names(ident.tbl)] <- ident.tbl
    ident.counts/k
  })
  colnames(cell.ident.frac) <- colnames(query.data)

  return(cell.ident.frac)
}

MapCellsKNN <- compiler::cmpfun(MapCellsKNN)


#' Summarize cell classification by group
#'
#' @param cell.ident.frac A reference groups x cells matrix of reference group probabilities
#' @param query.ident Either a list or a factor of cell type identities for the query dataset
#' @param ident.cell.thresh Filter out cells that don't map to a reference group with at least this probability
#' @param hard.cell.thresh Assigns each cell to a single reference group
#' @param min.ref.frac Filter out reference groups that are not mapped to with at least this probability
#'
#' @return A reference groups x query groups matrix of probabilities
#' @export
#'
SummarizeCellMapping <- function(cell.ident.frac, query.ident, ident.cell.thresh = 0, hard.cell.thresh = T,
                                 min.ref.frac = 0.05) {
  cell.ident.frac <- cell.ident.frac[,apply(cell.ident.frac, 2, max) > ident.cell.thresh]

  if(hard.cell.thresh) {
    cell.ident.frac <- apply(cell.ident.frac, 2, function(x) {
      y <- rep(0, length(x)); names(y) <- names(x);
      y[which.max(x)] <- 1
      y
    })
  }

  if (typeof(query.ident) == "list") {
    query.ident <- lapply(query.ident, function(cells) intersect(cells, colnames(cell.ident.frac)))
    query.ident.frac <- apply(cell.ident.frac, 1, function(x) {
      sapply(query.ident, function(cells) mean(x[cells]))
    })
  } else {
    query.ident <- factor(query.ident)
    query.ident <- query.ident[colnames(cell.ident.frac)]
    query.ident.frac <- apply(cell.ident.frac, 1, function(x) tapply(x, query.ident, mean))
  }

  query.ident.frac <- t(query.ident.frac)
  return(query.ident.frac[apply(query.ident.frac, 1, max) >= min.ref.frac,])
}

SummarizeCellMapping <- compiler::cmpfun(SummarizeCellMapping)


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
