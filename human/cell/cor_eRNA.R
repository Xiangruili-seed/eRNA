library(BiocFileCache)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(psych)

.validate_dimensions <- function(dims){
    if (length(dims) != 2) stop("Only two dimensions can be plotted")
}


.search_dimensions <- function(dims, cell_embeddings, reduction) {
    # Validate dimensions
    i <- dims %in% seq_len(ncol(cell_embeddings))
    if (!all(i)) {
        missing_dims <- dims[which(!i)]
        stop("Dimension(s) ", missing_dims, " not present in", reduction, "\n")
    }

    cell_embeddings <- cell_embeddings[, dims]
    cell_embeddings
}

wkde2d <- function(x, y, w, h, adjust = 1, n = 100,
                   lims = c(range(x), range(y))) {

    # Validate values and dimensions
    nx <- length(x)
    if (!all(all(nx == length(y)), all(nx == length(w)))) {
        stop("data vectors must be the same length")
    }
    if (any(!is.finite(x)) || any(!is.finite(y))) {
        stop("missing or infinite values in the data are not allowed")
    }
    if (any(!is.finite(lims))) {
        stop("only finite values are allowed in 'lims'")
    }

    h <- c(
        ks::hpi(x),
        ks::hpi(y)
    )
    h <- h * adjust

    # Get grid
    gx <- seq.int(lims[1L], lims[2L], length.out = n)
    gy <- seq.int(lims[3L], lims[4L], length.out = n)

    # weight
    ax <- outer(gx, x, "-") / h[1L]
    ay <- outer(gy, y, "-") / h[2L]

    w <- Matrix::Matrix(rep(w, n), nrow = n, ncol = nx, byrow = TRUE)

    z <- Matrix::tcrossprod(dnorm(ax) * w, dnorm(ay) * w) /
        (sum(w) * h[1L] * h[2L])

    dens <- list(x = gx, y = gy, z = z)
    dens
}


get_dens <- function(data, dens, method) {
    if (method == "ks") {
        ix <- findInterval(data[, 1], dens$eval.points[[1]])
        iy <- findInterval(data[, 2], dens$eval.points[[2]])
        ii <- cbind(ix, iy)
        z <- dens$estimate[ii]
    } else if (method == "wkde") {
        ix <- findInterval(data[, 1], dens$x)
        iy <- findInterval(data[, 2], dens$y)
        ii <- cbind(ix, iy)
        z <- dens$z[ii]
    }
    z
}


calculate_density <- function(w, x, method, adjust = 1, map = TRUE) {
    if (method == "ks") {
        dens <- kde(x[, c(1, 2)],
                    w = w / sum(w) * length(w)
        )
    } else if (method == "wkde") {
        dens <- wkde2d(
            x = x[, 1],
            y = x[, 2],
            w = w / sum(w) * length(w),
            adjust = adjust
        )
    }

    if (map) {
        get_dens(x, dens, method)
    } else {
        dens
    }
}

.extract_feature_data <- function(exp_data, features) {
    # Extract data for input features
    i <- colnames(exp_data) %in% features

    # Test existence of feature in gene expression data
    j <- !features %in% colnames(exp_data)
    if (any(j)) {
        stop(
            "'", paste(features[j], collapse = ", "),
            "' feature(s) not present in meta.data or expression data"
        )
    }
    vars <- exp_data[, i, drop = FALSE]
    vars <- vars[, features, drop = FALSE]
    vars
}
plot_density_ <- function(z, feature, cell_embeddings, dim_names, shape, size,
                          legend_title,
                          pal = c(
                              "viridis", "magma", "cividis",
                              "inferno", "plasma"
                          ), ...) {
    p <- ggplot(data.frame(cell_embeddings, feature = z)) +
        aes_string(dim_names[1], dim_names[2], color = "feature") +
        geom_point(shape = shape, size = size) +
        xlab(gsub("_", " ", dim_names[1])) +
        ylab(gsub("_", " ", dim_names[2])) +
        ggtitle(feature) +
        labs(color = guide_legend(legend_title)) +
        theme(
            text = element_text(size = 14),
            panel.background = element_blank(),
            axis.text.x = element_text(color = "black"),
            axis.text.y = element_text(color = "black"),
            axis.line = element_line(size = 0.25),
            strip.background = element_rect(color = "black", fill = "#ffe5cc")
        )

    pal <- match.arg(pal)
    p + scale_color_viridis_c(option = pal, ...)
    
}

#' @importFrom patchwork wrap_plots
.plot_final_density <- function(vars, cell_embeddings, features, joint, method,
                                adjust, shape, size, pal, combine, ...) {
    dim_names <- colnames(cell_embeddings)
    if (ncol(vars) > 1) {
        res <- apply(vars, 2, calculate_density, 
                     cell_embeddings, method, adjust)
        p <- mapply(plot_density_, as.list(as.data.frame(res)), colnames(res),
                    MoreArgs = list(cell_embeddings, dim_names, shape, size,
                                    "Density", pal = pal, ...), 
                    SIMPLIFY = FALSE)

        if(joint){

            z <- apply(res, 1, prod)
            joint_label <- paste0(paste(features, "+", sep = ""), 
                                  collapse = " ")
            pz <- plot_density_(z, joint_label, cell_embeddings,
                                dim_names,
                                shape,
                                size,
                                "Joint density",
                                pal = pal,
                                ...
            )


            if (combine) {
                p <- wrap_plots(p) + pz
            } else {
                p <- c(p, list(pz))
                names(p) <- c(features, joint_label)
            }
        }else{
            if (combine) {
                p <- wrap_plots(p)
            } else {
                names(p) <- features
            }
        }

    } else {
        z <- calculate_density(vars[, 1], cell_embeddings, method, adjust)
        p <- plot_density_(z,
                           features,
                           cell_embeddings,
                           dim_names,
                           shape,
                           size,
                           "Density",
                           pal = pal,
                           ...
        )
    }
    p
}


setGeneric("plot_density", function(object, features, slot = NULL,
                                    joint = FALSE, reduction = NULL,
                                    dims = c(1, 2),
                                    method = c("ks", "wkde"),
                                    adjust = 1, size = 1, shape = 16,
                                    combine = TRUE, pal = "viridis",
                                    ...)
    standardGeneric("plot_density"))

#' @importFrom Seurat GetAssayData Reductions Embeddings FetchData
#' @export
#' @describeIn plot_density Plot gene-weighted 2D kernel density
setMethod("plot_density", signature("Seurat"),
          function(object,
                   features, slot = NULL, joint = FALSE, reduction = NULL,
                   dims = c(1, 2), method = c("ks", "wkde"), adjust = 1,
                   size = 1, shape = 16, combine = TRUE, pal = "viridis",
                   ...) {

              # Validate dimensions -----
              .validate_dimensions(dims)
              # Validate existence of reduction
              if (!is.null(reduction)) {
                  if (!reduction %in% Reductions(object)) {
                      stop("No reduction named '", reduction,
                           "' found in object")
                  }
              }

              # Set up default reduction -----
              reduction_list <- Reductions(object)
              if (!length(reduction_list)){
                  stop("No reduction has been computed!")
              }
              if (is.null(reduction)) {
                  reduction <- reduction_list[length(reduction_list)]
              }
              cell_embeddings <- as.data.frame(Embeddings(object[[reduction]]))

              # Search for dimensions -----
              cell_embeddings <- .search_dimensions(dims, cell_embeddings,
                                                    reduction)
              # Set up default assay -----
              if (is.null(slot)) slot <- "data"

              # Match kde method
              method <- match.arg(method)

              # Extract metadata
              metadata <- object[[]]

              # Determine type of feature and extract
              if (all(features %in% colnames(metadata))) {
                  vars <- metadata[, features, drop = FALSE]
              }else{
                  exp_data <- FetchData(object, vars = features, slot = slot)
                  vars <- .extract_feature_data(exp_data, features)
              }
              z <- calculate_density(vars[, 1], cell_embeddings, method, adjust)
          })



args<-commandArgs(T)

dir<-args[1]
load(paste(dir,'/.Rdata',sep = ""))


data<-read.table(paste(dir,"s_gene_id.txt",sep = ""),header=F,stringsAsFactors = F)


fdr<-c()
co<-c()
for (i in 1:nrow(data)) {
  plot1<-plot_density(pbmc, as.character(data[i,1]),joint = F,method = "wkde",pal = "viridis")
  plot2<-plot_density(pbmc, as.character(data[i,2]),joint = F,method = "wkde",pal = "viridis")
  fdr[i]<-corr.test(plot1, plot2, method = "pearson", adjust = "fdr")$p
  co[i]<-cor(plot1,plot2,method="pearson")

}


result<-cbind(data,co,fdr)
p.adjust<-p.adjust(result[,4],method="fdr")
result1<-cbind(result,p.adjust)
write.table(result1,paste(dir,"/top_0.05_cor.txt",sep = ""),quote=F,col.names=F,row.names=F)


result2<-result1[which(result1[,3]>0.6),]
result3<-result2[which(result2[,5]<0.1),]

pdf(file = paste(dir,"/Nebulosa_cor.pdf",sep = ""))

for (i in 1:nrow(result3)) {
  plot1<-plot_density(pbmc, as.character(result3[i,1]),joint = F,method = "wkde",pal = "viridis")
  plot2<-plot_density(pbmc, as.character(result3[i,2]),joint = F,method = "wkde",pal = "viridis")
  date=as.data.frame(cbind(plot1,plot2))
  p<-ggplot(date,aes(x=plot1,y=plot2))+xlab(as.character(result2[i,1]))+ylab(as.character(result2[i,2]))+geom_point(size=1,shape=15)+ ggtitle(paste('correlation: ',result2[i,3],' & fdr: ',result2[i,4],sep = ""))
  print(p)

}
dev.off()

write.table(result3,paste(dir,"/top_0.05_cor_cut.txt",sep = ""),quote=F,col.names=F,row.names=F)
