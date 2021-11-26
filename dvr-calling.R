# dvr-calling.R
#
# Copyright (c) 2021 Yuttapong Thawornwattana
#
# Infer DVR deletion patterns of Mtb from read depth file 
# from mapping short reads to 14723_8_51.6


suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(gridExtra))


# filter options for DVR calling; see main() below
depth.q50.cutoff <- 5
depth.q25.cutoff <- 3
depth.q05.cutoff <- 1
prop.cutoff <- 0.5
depth.hi <- 50
depth.lo <- 15


# quantile functions
q01 <- function(x) quantile(x, 0.01) %>% unname
q05 <- function(x) quantile(x, 0.05) %>% unname
q10 <- function(x) quantile(x, 0.10) %>% unname
q25 <- function(x) quantile(x, 0.25) %>% unname
q75 <- function(x) quantile(x, 0.75) %>% unname
q90 <- function(x) quantile(x, 0.90) %>% unname
q95 <- function(x) quantile(x, 0.95) %>% unname


# define DR & DVR regions
get_dvr_def <- function() {
  ll <- read_lines(file.path("~/proj/mtb-l1/dvr", "14723_8_51.6_split.txt"))
  n.ch <- sapply(ll, nchar) %>% unname
  
  dr.dvr <- tibble(end=cumsum(n.ch)) %>%
    mutate(start=lag(end, default=0) + 1) %>%
    select(start, end)
  
  dr.dvr %>% 
    # remove dr, only keep dvr
    filter(row_number() %% 2 != 0) %>% 
    
    # remove starting and ending flanking regions
    mutate(length=end - start + 1) %>% filter(length < 50) %>% select(-length) %>%
    
    # add dvr id
    mutate(id=row_number(), dvr_id=rev(id))
}


get_dvr_position <- function(dvr.info) {
  foreach(i=seq_len(nrow(dvr.info)), .combine=rbind) %do% {
    dvr.i <- dvr.info %>% slice(i)
    tibble(
      pos=seq(dvr.i$start, dvr.i$end),
      dvr_id=dvr.i$dvr_id
    )
  }
}


# return names of spo positions (using letters A-O) for 68 DVR positions
get_spo_matching_dvr <- function() {
  out <- rep("", 68)
  out[2:4]   <- rep(LETTERS[1], 3)
  out[12:14] <- rep(LETTERS[2], 3)
  out[c(15, 18, 19)] <- rep(LETTERS[3], 3)
  out[20:43] <- rep(LETTERS[4:11], each=3)
  out[c(44, 46, 47)] <- rep(LETTERS[12], 3)
  out[51:53] <- rep(LETTERS[13], 3)
  out[62:64] <- rep(LETTERS[14], 3)
  out[65] <- LETTERS[15]
  tibble(dvr_id=1:68, spo_name=out)
}


# main function that extract DVRs from read depths
extract_dvr <- function(fin, is.plot=FALSE, fig.dir="") {
  # length of sq 14723_8_51.6 used as ref in read mapping
  n.ref.site <- 7316
  pos.all <- tibble(pos=seq_len(n.ref.site))
  
  # positions of dvr
  dvr.info <- get_dvr_def()
  dvr.pos <- get_dvr_position(dvr.info)
  
  spo.pos <- get_spo_matching_dvr()
  dvr.info <- left_join(dvr.info, spo.pos, by="dvr_id")
  
  d.col.type <- cols(
    ref=col_character(), 
    pos=col_double(), 
    depth=col_double()
  )

  if (is.plot) {
    plot_theme <- list(
      theme_bw(),
      theme(strip.background=element_blank(), panel.grid=element_blank(),
              legend.position="top"),
      scale_x_continuous(expand=c(0, 0), limits=c(1, n.ref.site), n.breaks=10),
      scale_y_continuous(expand=c(0, 0)),
      xlab("Position"), ylab("Depth"))
  }
  
  foreach(f=fin, .combine=rbind, .errorhandling="remove") %do% {
    # get read depth data
    d <- read_tsv(f, col_names=c("ref", "pos", "depth"), d.col.type) %>% 
      select(-ref)
    sample.id <- basename(f) %>% sub("_depth.zip", "", .) %>% 
      sub("_trim", "", .) %>% sub("_markdup", "", .)

    if (is.plot) {
      p1 <- ggplot(d) + 
        geom_rect(data=dvr.info, aes(xmin=start, xmax=end, ymin=0, ymax=Inf), 
                  colour="grey90", fill="grey90") +
        geom_text(data=dvr.info, size=2, check_overlap=TRUE,
                  aes(label=dvr_id, x=start + (end - start)/2, y=1)) +
        geom_text(data=dvr.info, size=2,
                  aes(label=spo_name, x=start + (end - start)/2, y=max(d$depth))) +
        geom_hline(yintercept=5, linetype="dashed") +
        geom_point(aes(pos, depth), size=0.2) +
        ggtitle(sample.id) +
        plot_theme
    }
    
    
    # step 1: fill in 0 at positions without depth
    d.full <- d %>% 
      right_join(pos.all, by="pos") %>%
      mutate(depth=ifelse(is.na(depth), 0, depth))
    
    
    # step 2: add dvr id to each position
    d.dvr <- d.full %>% 
      left_join(dvr.pos, by="pos") %>%
      filter(!is.na(dvr_id)) %>%
      group_by(dvr_id)

    if (is.plot) {
      p2 <- ggplot(d.dvr) + 
        geom_rect(data=dvr.info, aes(xmin=start, xmax=end, ymin=0, ymax=Inf), 
                  colour="grey90", fill="grey90") +
        geom_text(data=dvr.info, size=2, check_overlap=TRUE,
                  aes(label=dvr_id, x=start + (end - start)/2, y=1)) +
        geom_text(data=dvr.info, size=2, 
                  aes(label=spo_name, x=start + (end - start)/2, y=max(d.dvr$depth))) +
        geom_hline(yintercept=5, linetype="dashed") +
        geom_boxplot(aes(pos, depth, group=dvr_id), outlier.size=0.5) +
        plot_theme
      
      g <- arrangeGrob(p1, p2, ncol=1)
      ggsave(file.path(fig.dir, paste0(sample.id, ".pdf")), 
        width=12, height=5, plot=g)

      # log scales
      g <- suppressWarnings(suppressMessages(
        arrangeGrob(
          p1 + scale_y_log10(), 
          p2 + scale_y_log10(), 
          ncol=1
        )
      ))
      ggsave(file.path(fig.dir, paste0(sample.id, "-log.pdf")), 
        width=12, height=5, plot=g)
    }
    
    
    # step 3: summarise depth at each dvr region using quantiles
    # count of positions with non-zero depth
    d.dvr.n <- d.dvr %>% ungroup %>% 
      mutate(dvr_id=factor(dvr_id)) %>%
      filter(depth > 0) %>%
      count(dvr_id, name="n_pos", .drop=FALSE) %>%
      mutate(dvr_id=as.integer(as.vector(dvr_id)))
    
    # denominator (legnth of each dvr)
    d.dvr.n.all <- d.dvr %>% ungroup %>% 
      mutate(dvr_id=factor(dvr_id)) %>%
      count(dvr_id, name="n_all", .drop=FALSE) %>%
      mutate(dvr_id=as.integer(as.vector(dvr_id)))
    
    d.dvr.sum <- d.dvr %>% 
      summarise_at(vars(depth), lst(q05, q10, q25, median, q75, q90, q95)) %>%
      left_join(d.dvr.n, by="dvr_id") %>%
      left_join(d.dvr.n.all, by="dvr_id") 
    
    # add sample-level summaries, pooling across dvr regions
    d.dvr.sum.med.all <- d.dvr %>% ungroup %>% 
      summarise_at(vars(depth), list(median_all=median, mean_all=mean,
                                     q05_all=q05, q10_all=q10, q25_all=q25))
    d.dvr.sum.med.all.pos <- d.dvr %>% ungroup %>% filter(depth > 0) %>%
      summarise_at(vars(depth), list(median_all_pos=median, mean_all_pos=mean,
                                     q05_all_pos=q05, q10_all_pos=q10, 
                                     q25_all_pos=q25))
    
    # step 4: report dvr as deleted if median depth drops below the threshold
    # n_prop = proportion of positions with depth > 0 in each dvr
    d.dvr.sum %>% 
      mutate(id=sample.id, n_prop=n_pos / n_all) %>%
      select(id, everything()) %>%
      add_column(d.dvr.sum.med.all, d.dvr.sum.med.all.pos)
  }
}


# convert list of deleted DVRs to binary DVR pattern
# input: x = dvr_del from summarise_dvr()
dvr_del_to_bin <- function(dvr.del) {
  dvr.del <- dvr.del %>% strsplit(",") %>% unlist %>% as.integer
  
  dvr.bin <- rep(1L, 68)
  dvr.bin[dvr.del] <- 0

  dvr.bin
}
dvr_del_to_bin_vec <- Vectorize(dvr_del_to_bin)


spo_bin_to_oct <- function(x) {
  # check input
  if (length(x) == 1 && nchar(x) == 43) {
    x <- str_split(x, "") %>% unlist
  }
  if (!is.character(x)) x <- as.character(x)
  stopifnot(length(x) == 43)
  stopifnot(any(!is.na(match(x, c("0", "1")))))
  
  x %>% 
    head(-1) %>% 
    matrix(ncol=3, byrow=TRUE) %>%
    apply(1, function(x) paste(x, collapse="") %>% strtoi(base=2L)) %>% 
    paste(collapse="") %>% 
    paste0(tail(x, 1))  # add 43th value
}


# convert binary DVR pattern to spoligotype
# dvr.bin is the output from dvr_del_to_bin()
dvr_bin_to_spo <- function(dvr.bin, convert=FALSE) {
  spo.pos <- c(2:4, 12:15, 18:44, 46:47, 51:53, 62:65)

  if (is.character(dvr.bin)) {
    dvr.bin <- strsplit(dvr.bin, "") %>% unlist %>% as.integer
  }
  
  spo <- dvr.bin[spo.pos]
  stopifnot(identical(length(spo), 43L))

  if (convert) {
    spo <- spo_bin_to_oct(spo)
  } else {
    spo <- paste(spo, collapse="")
  }

  spo
}
dvr_bin_to_spo_vec <- Vectorize(dvr_bin_to_spo)


# convert runs of consecutive integers into ranges
# adapted from https://stackoverflow.com/questions/16911773/collapse-runs-of-consecutive-numbers-to-ranges
dvr_int_range <- function(x) {
  # process input
  if (!all(sapply(x, is.integer))) {
    x <- str_extract_all(x, "[0-9]+") %>% unlist %>% as.integer %>% sort
  }
  
  x.diff <- c(1, diff(x))
  diff.list <- split(x, cumsum(x.diff != 1))
  
  lapply(diff.list, function(x) {
    if (length(x) == 1) as.character(x) else paste0(x[1], "-", x[length(x)])
  }) %>% 
    unlist(use.names=FALSE) %>%
    paste0(collapse=",")
}
dvr_int_range <- Vectorize(dvr_int_range)


# summarise dvr outputs
summarise_dvr <- function(out) {
  out %>% 
  mutate(summary=paste0(dvr_id, " (q: ", q25, ",", median, ",", q75, ")")) %>%
  group_by(id) %>%
  summarise(
    dvr_n_del=n(),
    dvr_n=68 - dvr_n_del,
    dvr_del=paste(dvr_id, collapse=","),
    dvr_del_stat=paste(summary, collapse=", "),
    dvr=dvr_del %>% dvr_del_to_bin %>% paste(collapse=""),
    spo=dvr_del %>% dvr_del_to_bin %>% dvr_bin_to_spo(convert=TRUE)
  ) %>%
    mutate(dvr_del_range=dvr_int_range(dvr_del)) %>%
    select(id, spo, dvr, dvr_n_del, dvr_del_range, everything())
}


# column types for reading dvr data
dvr.cols <- cols(
  id = col_character(),
  spo = col_character(),
  dvr = col_character(),
  dvr_n = col_character(),
  dvr_n_del = col_character(),
  dvr_del = col_character(),
  dvr_del_stat = col_character()
)


main <- function() {
  # parse input arguments
  args <- commandArgs(trailingOnly=TRUE)
  if (length(args) >= 2) {
    din <- args[1]
    dout <- args[2]

    is.plot <- ifelse(any(grepl("--plot", args)), TRUE, FALSE)

  } else {
    stop("usage: Rscript dvr-calling.R input_dir output_dir")
  }

  stopifnot(dir.exists(din))
  if (!dir.exists(dout)) dir.create(dout)

  dfig <- file.path(dout, "fig")
  if (is.plot) {
    if (!dir.exists(dfig)) dir.create(dfig)  
  }

  fin <- list.files(din, "*_depth.zip", full.names=TRUE)


  fout <- file.path(dout, "dvr-full.csv")
  ffilt <- file.path(dout, "dvr-filt.csv")
  fsum <- file.path(dout, "dvr.csv")

  out <- extract_dvr(fin, is.plot=is.plot, fig.dir=dfig)
  out %>% write_csv(fout)

  out.filt <- out %>% 
    mutate(cutoff=pmin(depth.q50.cutoff, median_all / 2), 
           cutoff_hi=median_all / 4) %>%
    filter((median < cutoff & n_prop < prop.cutoff) | 
             median_all > depth.lo & median < cutoff |
             q25 < depth.q25.cutoff |
             q05 < depth.q05.cutoff |
             median_all > depth.hi & median < cutoff_hi) %>% 
    filter(n_prop < 1 | median <= 2)


  # format output
  out.filt %>% write_csv(ffilt)

  out.fmt <- summarise_dvr(out.filt)
  out.fmt %>% write_csv(fsum)
}

main()

