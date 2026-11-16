#!/usr/bin/env Rscript

# Rank reference populations and other CaG clusters for every CaG cluster at
# one ADMIXTURE K.
# Uses .Q files (individual proportions); .P files are not appropriate here.

d <- list(
  admixture_dir="/lustre07/scratch/justinp/IBD_clusters_project/IBD_segment_pipeline/ADMIXTURE/results/06_admixture",
  cag_metadata="/home/justinp/scratch/IBD_clusters_project/IBD_segment_pipeline/IBD-segment/results_3cM_Fst_0001/clustering/leiden/leiden.final_membership.POP_for_ADMIXTURE.tsv",
  reference_metadata="/lustre06/project/6072431/shared/REF_for_PCA/HGDP_1000G/metadata/HGDP_1KG.PassedQC.subpop.metadata.tsv",
  output_dir="/lustre07/scratch/justinp/IBD_clusters_project/IBD_segment_pipeline/ADMIXTURE/post_hoc_analyses/ranked_reference_pop_and_clusters",
  k=7L, top_n=10L, plot_top_n=5L, min_reference_n=2L)

help <- function() cat(
  "Usage: Rscript rank_reference_populations_and_clusters.R [options]\n",
  "Default: K=7. The selected K directory must contain admixture_input.fam\n",
  "and admixture_input.K.Q. Options: --admixture-dir, --cag-metadata,\n",
  "--reference-metadata, --output-dir, --k, --top-n, --plot-top-n,\n",
  "--min-reference-n.\n", sep="")

args <- commandArgs(trailingOnly=TRUE)
if (any(args %in% c("-h","--help"))) { help(); quit(status=0) }
if (length(args) %% 2) stop("Every option requires a value")
if (length(args)) {
  nm <- gsub("-","_",sub("^--","",args[seq(1,length(args),2)]),fixed=TRUE)
  if (any(!nm %in% names(d))) stop("Unknown option: ", nm[!nm %in% names(d)][1])
  val <- args[seq(2,length(args),2)]
  for (i in seq_along(nm)) d[[nm[i]]] <- val[i]
}
for (x in c("k","top_n","plot_top_n","min_reference_n")) d[[x]] <- as.integer(d[[x]])
if (anyNA(unlist(d[c("k","top_n","plot_top_n","min_reference_n")])))
  stop("Invalid numeric option")

read_meta <- function(path, type) {
  x <- read.delim(path, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE,
                  quote="", comment.char="")
  ids <- c("sample_id","IID","iid","sample","Sample","ID","id")
  labs <- if (type=="CaG") c("population","cluster","final_cluster","community","POP","pop")
          else c("population","subpopulation","subpop","POP","pop")
  ic <- ids[ids %in% names(x)][1]; lc <- labs[labs %in% names(x)][1]
  if (is.na(ic) || is.na(lc)) {
    x <- read.delim(path, header=FALSE, stringsAsFactors=FALSE, quote="", comment.char="")
    if (ncol(x)<2) stop(type," metadata needs at least two columns")
    message(type," headers not recognized; using columns 1 and 2")
    y <- data.frame(sample_id=x[[1]], label=x[[2]], stringsAsFactors=FALSE)
  } else y <- data.frame(sample_id=x[[ic]], label=x[[lc]], stringsAsFactors=FALSE)
  y$sample_id <- trimws(as.character(y$sample_id)); y$label <- trimws(as.character(y$label))
  y <- y[!is.na(y$sample_id)&nzchar(y$sample_id)&!is.na(y$label)&nzchar(y$label),]
  if (anyDuplicated(y$sample_id)) stop("Duplicate IDs in ",type," metadata")
  y$source <- type; y
}

norm_rows <- function(x) {
  if (any(!is.finite(x)) || any(x<0) || any(rowSums(x)<=0)) stop("Invalid Q proportions")
  x/rowSums(x)
}
centroids <- function(q, lab) {
  g <- split(seq_len(nrow(q)),lab)
  m <- t(vapply(g,function(i) colMeans(q[i,,drop=FALSE]),numeric(ncol(q))))
  m <- norm_rows(m)
  disp <- vapply(names(g),function(z) mean(sqrt(rowSums(sweep(q[g[[z]],,drop=FALSE],2,m[z,],"-")^2))),numeric(1))
  list(q=m,n=lengths(g),disp=disp)
}

# Return Cluster_1, Cluster_2, ..., Cluster_10 instead of lexicographic order.
natural_cluster_order <- function(labels) {
  suffix <- suppressWarnings(as.integer(sub("^.*_([0-9]+)$","\\1",labels)))
  order(is.na(suffix),suffix,labels)
}

order_cluster_centroids <- function(x) {
  idx <- natural_cluster_order(rownames(x$q))
  x$q <- x$q[idx,,drop=FALSE]
  x$n <- x$n[idx]
  x$disp <- x$disp[idx]
  x
}
jsd <- function(p,q) {
  m <- (p+q)/2
  kl <- function(a,b) sum(ifelse(a>0,a*log(a/b),0))
  sqrt((.5*kl(p,m)+.5*kl(q,m))/log(2))
}
compare <- function(a,b,k) {
  out <- vector("list",nrow(a$q)*nrow(b$q)); z <- 1L
  for(i in seq_len(nrow(a$q))) for(j in seq_len(nrow(b$q))) {
    p<-a$q[i,]; q<-b$q[j,]; jd<-jsd(p,q)
    out[[z]] <- data.frame(K=k,dataset_cluster=rownames(a$q)[i],dataset_n=unname(a$n[i]),
      dataset_dispersion=unname(a$disp[i]),reference_population=rownames(b$q)[j],
      reference_n=unname(b$n[j]),reference_dispersion=unname(b$disp[j]),
      js_similarity=1-jd,js_distance=jd,
      hellinger_distance=sqrt(sum((sqrt(p)-sqrt(q))^2))/sqrt(2),
      euclidean_distance=sqrt(sum((p-q)^2)),
      cosine_similarity=sum(p*q)/(sqrt(sum(p^2))*sqrt(sum(q^2))))
    z<-z+1L
  }
  x<-rbind_fill(out)
  cluster_position <- match(x$dataset_cluster,rownames(a$q))
  x<-x[order(cluster_position,-x$js_similarity,x$reference_population),]
  x$rank<-ave(x$js_similarity,x$dataset_cluster,FUN=seq_along)
  x[c("K","dataset_cluster","dataset_n","dataset_dispersion","rank","reference_population",
      "reference_n","reference_dispersion","js_similarity","js_distance","hellinger_distance",
      "euclidean_distance","cosine_similarity")]
}
w <- function(x,p) write.table(x,p,sep="\t",quote=FALSE,row.names=FALSE,na="NA")

# Stack tables that contain different Q columns (e.g., K2 has Q1-Q2 while
# K9 has Q1-Q9). Columns absent at a given K are filled with NA.
rbind_fill <- function(tables) {
  all_names <- unique(unlist(lapply(tables,names),use.names=FALSE))
  aligned <- lapply(tables,function(x) {
    missing <- setdiff(all_names,names(x))
    for(nm in missing) x[[nm]] <- NA_real_
    x[all_names]
  })
  result <- aligned[[1L]]
  if(length(aligned)>1L) {
    for(i in 2:length(aligned)) result <- rbind(result,aligned[[i]])
  }
  result
}

needed <- c(d$cag_metadata,d$reference_metadata)
if (any(!file.exists(needed))) stop("Missing file(s):\n",paste(needed[!file.exists(needed)],collapse="\n"))
cm <- read_meta(d$cag_metadata,"CaG"); rm <- read_meta(d$reference_metadata,"Reference")
if(length(intersect(cm$sample_id,rm$sample_id))) stop("IDs overlap between metadata files")
all_meta <- rbind_fill(list(cm,rm))
dir.create(d$output_dir,recursive=TRUE,showWarnings=FALSE)
k <- d$k
{
  kdir <- file.path(d$admixture_dir,paste0("K",k))
  qf <- file.path(kdir,paste0("admixture_input.",k,".Q"))
  fam_file <- file.path(kdir,"admixture_input.fam")
  missing_k_files <- c(qf,fam_file)[!file.exists(c(qf,fam_file))]
  if(length(missing_k_files)) {
    stop("Missing input for K",k,": ",paste(missing_k_files,collapse=", "))
  }
  message("Processing K",k)
  fam <- read.table(fam_file,header=FALSE,stringsAsFactors=FALSE,quote="",comment.char="")
  if(ncol(fam)<2) stop("Invalid FAM file: ",fam_file)
  ids <- as.character(fam[[2]])
  if(anyDuplicated(ids)) stop("Duplicate IDs in ",fam_file)
  meta <- all_meta[match(ids,all_meta$sample_id),]
  known <- !is.na(meta$sample_id)
  if(any(!known)) warning("K",k,": ",sum(!known)," FAM samples lack metadata and will be excluded")
  q <- as.matrix(read.table(qf,header=FALSE)); storage.mode(q)<-"double"
  if(nrow(q)!=length(ids)||ncol(q)!=k) stop("Unexpected dimensions for ",qf)
  colnames(q)<-paste0("Q",seq_len(k)); q<-norm_rows(q)
  ca<-centroids(q[known&meta$source=="CaG",,drop=FALSE],meta$label[known&meta$source=="CaG"])
  ca<-order_cluster_centroids(ca)
  re<-centroids(q[known&meta$source=="Reference",,drop=FALSE],meta$label[known&meta$source=="Reference"])
  keep<-re$n>=d$min_reference_n; re$q<-re$q[keep,,drop=FALSE]; re$n<-re$n[keep]; re$disp<-re$disp[keep]
  # Use one unified comparison pool: all reference populations plus all CaG
  # clusters. The focal cluster itself is removed below.
  comparison<-list(
    q=rbind(re$q,ca$q),
    n=c(re$n,ca$n),
    disp=c(re$disp,ca$disp))
  comparison_type<-c(setNames(rep("reference_population",nrow(re$q)),rownames(re$q)),
                     setNames(rep("CaG_cluster",nrow(ca$q)),rownames(ca$q)))
  ranked<-compare(ca,comparison,k)
  names(ranked)[names(ranked)=="reference_population"]<-"comparison_group"
  names(ranked)[names(ranked)=="reference_n"]<-"comparison_n"
  names(ranked)[names(ranked)=="reference_dispersion"]<-"comparison_dispersion"
  ranked$comparison_type<-unname(comparison_type[ranked$comparison_group])
  ranked<-ranked[!(ranked$comparison_type=="CaG_cluster" &
                     ranked$dataset_cluster==ranked$comparison_group),]
  cluster_position<-match(ranked$dataset_cluster,rownames(ca$q))
  ranked<-ranked[order(cluster_position,-ranked$js_similarity,ranked$comparison_group),]
  ranked$rank<-ave(ranked$js_similarity,ranked$dataset_cluster,FUN=seq_along)
  ct<-rbind_fill(list(
    data.frame(K=k,type="dataset_cluster",group=rownames(ca$q),n=unname(ca$n),dispersion=unname(ca$disp),ca$q,check.names=FALSE),
    data.frame(K=k,type="reference_population",group=rownames(re$q),n=unname(re$n),dispersion=unname(re$disp),re$q,check.names=FALSE)))
  prefix<-paste0("K",k,".")
  w(ranked,file.path(d$output_dir,paste0(prefix,"population_and_cluster_similarity.all_ranked.tsv")))
  w(ranked[ranked$rank<=d$top_n,],
    file.path(d$output_dir,paste0(prefix,"population_and_cluster_similarity.top",d$top_n,".tsv")))
  w(ct,file.path(d$output_dir,paste0(prefix,"admixture_centroids.tsv")))
  w(ranked[ranked$rank==1,],file.path(d$output_dir,paste0(prefix,"best_population_or_cluster_match.tsv")))
}

if(!requireNamespace("ggplot2",quietly=TRUE))
  stop("The ggplot2 package is required to create plots. Install it or load an R module that provides it.")
library(ggplot2)

# Complete heatmap: every focal CaG cluster versus the unified comparison pool.
heat <- ranked
heat$dataset_cluster <- factor(heat$dataset_cluster,levels=unique(heat$dataset_cluster))
reference_heat <- heat[heat$comparison_type=="reference_population",]
ref_order <- aggregate(js_similarity~comparison_group,reference_heat,mean)
ref_order <- ref_order$comparison_group[order(ref_order$js_similarity)]
comparison_levels <- c(ref_order,rev(rownames(ca$q)))
heat$comparison_group <- factor(heat$comparison_group,levels=comparison_levels)
p_heat <- ggplot(heat,aes(dataset_cluster,comparison_group,fill=js_similarity))+
  geom_tile()+
  scale_fill_viridis_c(option="magma",direction=1,limits=c(0,1),name="JS similarity")+
  labs(title=paste0("Population and cluster similarity at K=",k),
       subtitle="Reference populations and other CaG clusters",x="Focal CaG cluster",y="Comparison group")+
  theme_minimal(base_size=11)+
  theme(axis.text.x=element_text(angle=45,hjust=1),panel.grid=element_blank())
ggsave(file.path(d$output_dir,paste0(prefix,"similarity_heatmap.png")),p_heat,
       width=max(10,0.45*length(unique(heat$dataset_cluster))),
       height=min(24,max(8,0.25*length(unique(heat$comparison_group)))),dpi=300,limitsize=FALSE)

# Top-match overview. References and other CaG clusters compete in one ranking.
top_plot <- ranked[ranked$rank<=d$plot_top_n,]
top_plot <- top_plot[order(top_plot$dataset_cluster,-top_plot$rank),]
top_plot$dataset_cluster <- factor(top_plot$dataset_cluster,levels=rownames(ca$q))
top_plot$key <- factor(paste(top_plot$dataset_cluster,top_plot$comparison_group,sep="|||"),
                       levels=unique(paste(top_plot$dataset_cluster,top_plot$comparison_group,sep="|||")))
p_top <- ggplot(top_plot,aes(js_similarity,key))+
  geom_segment(aes(x=0,xend=js_similarity,yend=key),colour="#B3C3C8",linewidth=0.6)+
  geom_point(size=2.4,colour="#2F4858")+
  facet_wrap(~dataset_cluster,scales="free_y")+
  scale_y_discrete(labels=function(x) sub("^.*\\|\\|\\|","",x))+
  scale_x_continuous(limits=c(0,1))+
  labs(title=paste0("Top ",d$plot_top_n," population or cluster matches at K=",k),
       x="Jensen-Shannon similarity",y=NULL)+
  theme_minimal(base_size=10)+theme(panel.grid.major.y=element_blank(),legend.position="none")
ggsave(file.path(d$output_dir,paste0(prefix,"top_population_and_cluster_matches.png")),p_top,
       width=14,height=max(8,2.2*ceiling(length(unique(top_plot$dataset_cluster))/4)),
       dpi=300,limitsize=FALSE)

# One profile plot per CaG cluster versus its leading groups of either type.
for(cluster in rownames(ca$q)) {
  refs <- ranked$comparison_group[ranked$dataset_cluster==cluster & ranked$rank<=d$plot_top_n]
  groups <- c(cluster,refs)
  mats <- rbind(ca$q[cluster,,drop=FALSE],comparison$q[refs,,drop=FALSE])
  profile <- data.frame(group=rep(groups,each=k),component=rep(colnames(mats),times=length(groups)),
                        proportion=as.vector(t(mats)),stringsAsFactors=FALSE)
  profile$group <- factor(profile$group,levels=groups)
  reference_colours <- setNames(scales::hue_pal()(length(refs)),refs)
  plot_colours <- c(setNames("black",cluster),reference_colours)
  p_profile <- ggplot(profile,aes(component,proportion,group=group,colour=group))+
    geom_line(aes(linewidth=group==cluster))+geom_point(aes(size=group==cluster))+
    scale_colour_manual(values=plot_colours)+
    scale_linewidth_manual(values=c(`TRUE`=1.3,`FALSE`=0.7),guide="none")+
    scale_size_manual(values=c(`TRUE`=2.8,`FALSE`=1.8),guide="none")+
    scale_y_continuous(limits=c(0,1),expand=expansion(mult=c(0,0.04)))+
    labs(title=paste0(cluster," versus its top comparison groups (K=",k,")"),
         subtitle="Thicker line: CaG cluster centroid",x="ADMIXTURE component",y="Mean proportion",colour="Group")+
    theme_minimal(base_size=11)+theme(panel.grid.minor=element_blank(),legend.position="bottom")
  safe_cluster <- gsub("[^A-Za-z0-9_.-]+","_",cluster)
  ggsave(file.path(d$output_dir,paste0(prefix,safe_cluster,".top_comparison_profiles.png")),
         p_profile,width=9,height=6,dpi=300)
}


message("Done: ",d$output_dir)
