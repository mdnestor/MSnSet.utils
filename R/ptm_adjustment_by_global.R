
#' Adjusting PTM profiles by global
#'
#'
#' Adjusting PTM relative abundance profiles by the global profiles of the
#' parent protein. The motivation behind this adjustment is that there are two
#' contributions to the observed PTM abundance changes - change in the
#' PTM stoichiometry and the change in the parent protein abundance.
#'
#'
#' @param m_ptm MSnSet with PTM data
#' @param m_global MSnSet with global data
#'
#'
#' @return MSnSet object
#' @importFrom Biobase fData featureNames exprs
#' @export adjust_for_parent_protein_profiles
#'
#'


adjust_for_parent_protein_profiles <- function(m_ptm, m_global){

   # sample alignment
   common_samples <- intersect(sampleNames(m_ptm), sampleNames(m_global))
   stopifnot(length(common_samples) > 2)
   m_ptm <- m_ptm[,common_samples]
   m_global <- m_global[,common_samples]

   # gene matching
   stopifnot("Gene" %in% fvarLabels(m_ptm))

   gene_overlap <- intersect(fData(m_ptm)$Gene, featureNames(m_global))
   m_ptm <- m_ptm[fData(m_ptm)$Gene %in% gene_overlap,]
   x_global <- exprs(m_global)[fData(m_ptm)$Gene,]

   x_ptm <- exprs(m_ptm)

   pf <- function(i){
      x <- as.numeric(x_global[i,])
      y <- as.numeric(x_ptm[i,])
      idx <- !is.na(x) & !is.na(y)
      if(sum(idx) > 1){
         mdl <- lm(y ~ x)
         ynot <- predict(mdl, newdata = data.frame(x))
         return(y - as.numeric(ynot))
      }else{
         return(rep(NA,ncol(x_ptm)))
      }
   }

   res <- lapply(1:nrow(m_ptm), pf)

   x_ptm_adj <- matrix(unlist(res), byrow=TRUE, nrow=nrow(x_ptm))
   exprs(m_ptm) <- x_ptm_adj

   idx_NAs <- apply(is.na(exprs(m_ptm)), 1, all)

   return(m_ptm[!idx_NAs,])

}


