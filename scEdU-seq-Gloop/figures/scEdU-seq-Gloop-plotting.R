#!/usr/bin/env Rscript

# dependency

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))

mixdist   <- function(Estimate, diffusion){
  
  if(!diffusion){
    if(length(Estimate) == 4){
      
      X <- data.table(dist = seq(0, 4e5, 5000)
      )[,.(dist = dist,
           E = dexp(dist, 1 / Estimate[1]) %>% {. / sum(.) * Estimate[3]},
           U = dunif(dist, 0, 4e5) %>% {. / sum(.) * Estimate[4]} )] %>%
        melt(id.vars = "dist")
    } else
      if(length(Estimate) == 7){
        
        X <-  data.table(dist = seq(0, 4e5, 5000)
        )[,.(dist = dist,
             E = dexp(dist, 1 / Estimate[1]) %>% {. / sum(.) * Estimate[5]} ,
             M = dnorm(dist, Estimate[2], sqrt(Estimate[3])) %>% {. / sum(.) * Estimate[6] },
             U = dunif(dist, 0, 4e5) %>% {. / sum(.) * Estimate[7]})] %>%
          melt(id.vars = "dist")
      } else 
        if(length(Estimate) == 10){
          X <-   data.table(dist = seq(0, 4e5, 5000)
          )[,.(dist = dist,
               E = dexp(dist, 1 / Estimate[1]) %>% {. / sum(.) * Estimate[7]} ,
               M0 = dnorm(dist, Estimate[2], sqrt(Estimate[3])) %>% {. / sum(.) * Estimate[8] },
               M = dnorm(dist, Estimate[4], sqrt(Estimate[5])) %>% {. / sum(.) * Estimate[9] },
               U = dunif(dist, 0, 4e5) %>% {. / sum(.) * Estimate[10]})] %>%
            melt(id.vars = "dist")}
  } 
  
  if(diffusion){
    if(length(Estimate) == 4){
      
      X <- data.table(dist = seq(0, 4e5, 5000)
      )[,.(dist = dist,
           E = dexp(dist, 1 / Estimate[1]) %>% {. / sum(.) * Estimate[3]},
           U = dunif(dist, 0, 4e5) %>% {. / sum(.) * Estimate[4]} )] %>%
        melt(id.vars = "dist")
      
      
    } else 
      if(length(Estimate) == 7){
        
        X <- data.table(dist = seq(0, 10e5, 5000)
        )[,.(dist = dist,
             E = dexp(dist, 1 / Estimate[1]) %>% {. / sum(.) * Estimate[5]} ,
             M0 = dnorm(dist, Estimate[2], sqrt(Estimate[3])) %>% {. / sum(.) * Estimate[6] },
             U = dunif(dist, 0, 10e5) %>% {. / sum(.) * Estimate[7]})] %>%
          melt(id.vars = "dist")
        
      } else 
        if(length(Estimate) == 10){
          X <-   data.table(dist = seq(0, 10e5, 5000)
          )[,.(dist = dist,
               E = dexp(dist, 1 / Estimate[1]) %>% {. / sum(.) * Estimate[7]} ,
               M0 = dnorm(dist, Estimate[2], sqrt(Estimate[3])) %>% {. / sum(.) * Estimate[8] },
               M = dnorm(dist, Estimate[4], sqrt(Estimate[5])) %>% {. / sum(.) * Estimate[9] },
               U = dunif(dist, 0, 10e5) %>% {. / sum(.) * Estimate[10]})] %>%
            melt(id.vars = "dist")}
  } 
  
  
  return(X[,variable := factor(variable, c("M", "E","M0", "U"))][])
  
  
  
}
frollconf <- function(X, Y, N){
  Ypad <- c(rep(NA, round(1.5 * N)), Y, rep(NA, round(1.5 * N)))
  
  data.table(ymeanpad = frollmean(Ypad, n = 2 * N + 1, align = "center", na.rm = T, hasNA = T)[],
             y95pad = frollapply(Ypad, n = 2 * N + 1, FUN = sd, align = "center", na.rm = T) * 2
  )[,.(x = X,
       ymean = ymeanpad[!is.na(Ypad)],
       y95 = y95pad[!is.na(Ypad)],
       yse = c(y95pad / sqrt(frollsum(!is.na(Ypad) * 1, n = 2 * N + 1, align = "center")))[!is.na(Ypad)])
  ]
  
}

dependencies <- setNames(c("_filterer",    "_overlapscore", "_filterer", "_pc",      "_pc",    "_pc", "_filterer", "_hmmfit"),
                         c("overlapscore", "sphaseorder", "pc",       "psmashc", "mmfit", "BICfit", "hmmfit", "singlefork"))

# arguments

arguments <- commandArgs(trailingOnly=TRUE)

print(arguments)

exp_pattern <- "mESC"

base_folder <- "~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/"

pc_sim <- fread("~/Desktop/post-doc/Experiments/exp140_CAF1_10m-2hr-30hr-15E-60WO-15E-F/scEdU-seq/pc_sim_super_all.tsv") %>% 
  dplyr::filter(U !="U")%>%
  mutate_if(is.character, as.numeric)%>% as.data.table()

str(pc_sim)
pc_sim[, c("u_R", "s_R")       := .(abs(u - mean(u)) / sd(u), abs(s - mean(s)) / sd(s)), .(U,S)
][, c("u_rank", "s_rank") := .(rank(u_R), rank(s_R)), .(U,S)]

U_model <- loess(U ~ u + s, pc_sim[u_rank < 10 | s_rank < 10])
S_model <- loess(S ~ u + s, pc_sim[u_rank < 10 | s_rank < 10])
#
#all_plots <- list()

#for(i in 1:26) all_plots[[paste0(LETTERS, LETTERS)[i]]] <- "empty"

##### PLOT files
#cat("\n1. PLOT files \n")
#try({
AA <- data.table(file = list.files(base_folder, pattern = exp_pattern, recursive = T, include.dirs = F)
)[, exp := gsub(".*/|_.*|-pl.*","", file)
][, type := sub("_|-", "", gsub(paste0(exp, "|.tsv.gz|.RDS"), "", sub(".*/", "", file))), .(file)
][, size := file.size(paste0(base_folder, file)) / 1e6
][str_detect(type, "memstats|plotted", T) & str_detect(exp, "plot", T)] %>%
  ggplot() + geom_raster(aes(x = type, y = exp, fill = (size))) +
  facet_grid(cols = vars(str_detect(type, "pl")), scales = "free")
#})

AA
##### PLOT memory usage
#cat("\n2. PLOT memory usage \n")
#try({all_plots$
BB <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*memstats"), recursive = T), 
                    function(x){fread(paste0(base_folder, x))[, exp := gsub(".*/|_.*", "", x)]})
)[, job := sub("_.*", "", JobID)
][][, jobtype := sub(".*\\.", "", JobName)
][][,.(mem_used = sum(as.numeric(sub("K", "", MaxRSS)), na.rm = T) / 1e6,
       time_req = Timelimit[1], 
       jobtype = jobtype[1],
       time_elapsed = Elapsed[1]), .(exp, job, AllocCPUS, mem_req = as.numeric(sub("Gn", "", ReqMem)))
][][, dependency := grep(dependencies[jobtype], 
                         list.files(base_folder, pattern = exp, full.names = T, recursive = T), value = T)[1], .(exp, job)
][][, mem_depend := file.size(dependency) / 1e9
][!is.na(mem_depend), depend_ratio := coef(lm(mem_used ~ mem_depend, data = .SD))[2], .(jobtype)
][] %>%
  ggplot() + 
  #geom_point(aes(x = as.POSIXct(time_req), y = as.POSIXct(time_elapsed), col = jobtype)) + 
  geom_point(aes(x = mem_depend, y = mem_used, col = jobtype)) + 
  geom_smooth(aes(x = mem_depend, y = mem_used, col = jobtype), method = "lm", se = F) + 
  geom_label(aes(x = 0.6, y = as.numeric(as.factor(jobtype)) * 2, label = paste(jobtype, round(depend_ratio, 2)), col = jobtype)) +
  #geom_point(aes(x = mem_req, y = mem_used, col = jobtype)) + 
  geom_abline(slope = 1, intercept = 0)

BB
#})


##### PLOT poisson test
CC <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*poissontest"), recursive = T), 
                    function(x){fread(paste0(base_folder, x))[, exp := gsub(".*/|_.*", "", x)]})
) %>%
  ggplot() +
  geom_point(aes(x = log(mean), y = poisson_score, col = exp),
             size = 0.1) +
  facet_wrap(vars(exp), scales = "free") +
  theme(legend.position = "none") + 
  geom_hline(yintercept = 0.1) +
  geom_vline(xintercept = c(-2.5, 2))

CC
##### PLOT posterior per experiment
DD <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*filterer"), recursive = T), 
                    function(x){fread(paste0(base_folder, x))[, exp := gsub(".*/|_.*", "", x)]})
)[,.(N = as.double(.N)), .(exp, posterior)
][,N := N / max(N), .(exp)] %>%
  ggplot() + geom_tile(aes(x = 1, y = posterior, height = 1, width = N)) +
  facet_wrap(vars(exp))
DD

##### PLOT shpase order QCs
EE <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*sphaseorder"), recursive = T), 
                    function(x){readRDS(paste0(base_folder, x))$bootstrap[, exp := gsub(".*/|_.*", "", x)]})
) %>% 
  ggplot() + geom_point(aes(x = rank, y = value, col = N), size = 0.2) + 
  facet_wrap(vars(exp)) +
  scale_color_viridis_c()
EE

FF <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*sphaseorder"), recursive = T), 
                    function(x){as.data.table(readRDS(paste0(base_folder, x))$final_order)[, exp := gsub(".*/|_.*", "", x)]})
) %>% ggplot(aes(y = X1.1, x = X2.1)) +
  geom_segment(aes(yend = X2, xend = Y2), size = 0.05) +
  geom_point(aes(col = rank), size = 1) +
  coord_fixed() +
  facet_wrap(vars(exp)) +
  scale_color_viridis_c(name = "Tour index")
FF
####### plotting cut tracks
GG <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*filterer"), recursive = T), 
                    function(x){
                      fread(paste0(base_folder, x)
                      )[chr == 2 & posterior > 800
                      ][, exp := gsub(".*/|_.*", "", x)
                      ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphase"), recursive = T, full.names = T)
                      )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                      ][, .(N = as.double(.N)), .(exp, rank, bp = round(bp, -5))
                      ][, N := (N) / quantile(N, 0.80), .(rank)
                      ][bp > 4.5e7 & bp < 6.5e7
                      ][N > 1, N := 1]
                    })
)[str_detect(exp, "133-WT-mESC|133-FANCJKO-mESC|133-DHX36KO-mESC"), rank := abs(rank - 1)  #"DMSO-2hr|V1-2hr"
]%>% 
  ggplot() + 
  geom_raster(aes(y = rank, x = bp / 1e6, fill = (N)) ) +
  scale_y_continuous(trans = "reverse") + ylab("S-phase Progression") + xlab("Chromosome 2 (Mbp)") +
  facet_wrap(vars(fct_rev(exp)), scales = "free", ncol = 4) +
  scale_fill_gradientn(colours = c("white", "grey30", "grey5")) +
  theme_bw() +
  coord_cartesian(expand = F)
GG

####### plotting genome pile up
HH <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*filterer"), recursive = T), 
                    function(x){
                      fread(paste0(base_folder, x)
                      )[chr %in% 1:22, .(N = as.double(.N), exp = gsub(".*/|_.*", "", x)), .(cell, chr, bp = round(bp / 1e6))
                      ][, N := N / sum(N), .(cell)
                      ][, .(N = sum(N)), .(exp, chr, bp)
                      ][, N := N / quantile(N, 0.97), .(chr)
                      ][N > 1, N := 1][]
                    })) %>% 
  ggplot() + 
  geom_col(aes(y = N, x = bp, fill = log10(N)) ) +
  facet_grid(rows = vars(exp), cols = vars(chr), scales = "free") +
  scale_fill_viridis_c() +
  coord_cartesian(expand = F)
HH

####### plotting model fits
II <- merge.data.table(
  rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
                function(x){fread(paste0(base_folder, x)
                )[, exp := gsub(".*/|_.*", "", x)
                ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                ]})
  )[((!diffusion) & model != "quadruple") | (diffusion)
  ][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
  ][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
  ][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
  ][parameter != 11, mixdist(Estimate, diffusion = diffusion[1]), .(exp, cell, rank, diffusion)],
  rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*psmashc"), recursive = T), 
                function(x){fread(paste0(base_folder, x))[, exp := gsub(".*/|_.*", "", x)]})),
  by = c("exp", "cell", "dist")
)[dist > 0
][,  N := as.double(N)
][, N := N / sum(N), .(exp, cell, diffusion, variable)
][str_detect(exp, "DMSO|24"), rank := abs(rank - 1)   
  #][rank > 0.4 & rank < 0.6
][,.SD[cell %in% sample(unique(cell), 3)], .(exp)
][dist > 0 & dist < 4e5] %>%
  ggplot() + 
  geom_col(aes(x = dist, y = value, fill = variable)) + 
  geom_line(aes(x = dist , y = N)) +
  facet_wrap(vars(exp, rank), scales = "free", nrow = 4) 
#})

II

##### PLOT speeds
JJ <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
                    function(x){fread(paste0(base_folder, x)
                    )[, exp := gsub(".*/|_.*", "", x)
                    ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                    )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                    ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 3e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4 
][str_detect(exp,"133-WT-mESC|133-FANCJKO-mESC|133-DHX36KO-mESC"), rank := abs(rank - 1)  
][rank < 0.8
  #][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
] %>% 
  {
    .[!is.na(rank)][order(rank), frollconf(rank, Estimate, 15), .(exp)]  %>% 
      ggplot() + 
      geom_point(data = ., aes(x = rank, y = Estimate )) + 
      geom_ribbon(aes(x = x, ymin = ymean - y95, ymax = ymean + y95), alpha = 0.3) +
      geom_ribbon(aes(x = x, ymin = ymean - yse, ymax = ymean + yse), alpha = 0.3) +
      geom_line(aes(x = x, y = ymean)) +
      facet_grid(cols = vars(exp)) + theme_bw() + ylab("Speed (kb/min)") + xlab("S-phase progression")
    
  }

JJ

KK <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
                    function(x){fread(paste0(base_folder, x)
                    )[, exp := gsub(".*/|_.*", "", x)
                    ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                    )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                    ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4 
][str_detect(exp,"133-WT-mESC|133-FANCJKO-mESC|133-DHX36KO-mESC"), rank := abs(rank - 1)  
][rank < 0.8
  
][!is.na(rank)][order(rank), frollconf(rank, Estimate, 15), .(exp)]  %>% 
  ggplot() + 
  #geom_ribbon(aes(x = x, ymin = ymean - y95, ymax = ymean + y95, col = exp), alpha = 0.3) +
  geom_ribbon(aes(x = x, ymin = ymean - yse, ymax = ymean + yse, fill = exp), alpha = 0.3) + 
  geom_line(aes(x = x, y = ymean, col = exp)) + 
  theme_bw() + ylab("Speed (kb/min)") + xlab("S-phase progression")

KK

pc_sim <- fread("~/Desktop/post-doc/Experiments/exp140_CAF1_10m-2hr-30hr-15E-60WO-15E-F/scEdU-seq/pc_sim_super_all.tsv") %>% 
  dplyr::filter(U !="U")%>%
  mutate_if(is.character, as.numeric)%>% as.data.table()

str(pc_sim)
pc_sim[, c("u_R", "s_R")       := .(abs(u - mean(u)) / sd(u), abs(s - mean(s)) / sd(s)), .(U,S)
][, c("u_rank", "s_rank") := .(rank(u_R), rank(s_R)), .(U,S)]

U_model <- loess(U ~ u + s, pc_sim[u_rank < 10 | s_rank < 10])
S_model <- loess(S ~ u + s, pc_sim[u_rank < 10 | s_rank < 10])


LL <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
                    function(x){fread(paste0(base_folder, x)
                    )[, exp := gsub(".*/|_.*", "", x)
                    ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                    )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                    ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4 
][str_detect(exp, "133-WT-mESC|133-FANCJKO-mESC|133-DHX36KO-mESC"), rank := abs(rank - 1)  
][rank < 0.8
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
] %>% # dplyr::filter(str_detect(exp, "133-WT-mESC|133-DHX36-FANCJKO-mESC")) %>%
  {
    .[!is.na(rank)][order(rank), frollconf(rank, u, 30), .(exp)]  %>% 
      ggplot() + 
      geom_point(data = ., aes(x = rank, y = u, col = fct_rev(exp)), size=1, alpha=0.5) + 
      #geom_ribbon(aes(x = x, ymin = ymean - y95, ymax = ymean + y95, fill=exp), alpha = 0.2) +
      geom_ribbon(aes(x = x, ymin = ymean - yse, ymax = ymean + yse, fill=fct_rev(exp)), alpha = 0.2) +
      geom_line(aes(x = x, y = ymean, col=fct_rev(exp)), size = 1) +
      facet_grid(cols = vars(fct_rev(exp))) + theme_bw() + ylab("Speed (kb/min)") + xlab("S-phase progression")+
      scale_fill_manual(values = c( "#C0BDBA", "#6D90CA", "#F9A11B", 
                                    "#E92724"))+
      scale_color_manual(values = c("#C0BDBA", "#6D90CA", "#F9A11B", 
                                    "#E92724" ))+
      theme(legend.position = "bottom")+ylim(0.5,2.5)+
      coord_cartesian()
    
  }

LL

MM<- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
                   function(x){fread(paste0(base_folder, x)
                   )[, exp := gsub(".*/|_.*", "", x)
                   ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                   )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                   ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4 
][str_detect(exp,"133-WT-mESC|133-FANCJKO-mESC|133-DHX36KO-mESC"), rank := abs(rank - 1)  
][rank < 0.8
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
] %>% #dplyr::filter(str_detect(exp, "133-WT-mESC|133-DHX36-FANCJKO-mESC")) %>%
  #{
  # .[!is.na(rank)][order(rank), frollconf(rank, u, 25), .(exp)]  %>% 
  ggplot() + 
  
  geom_violin( aes(x = fct_rev(exp), y = u, group =fct_rev(exp)), size=0.33, alpha=0.5, trim = F)+
  geom_jitter( aes(x = fct_rev(exp), y = u, col =fct_rev(exp)), size=0.5, alpha=0.66) + 
  #geom_ribbon(aes(x = x, ymin = ymean - y95, ymax = ymean + y95, fill=exp), alpha = 0.2) +
  #geom_ribbon(aes(x = x, ymin = ymean - yse, ymax = ymean + yse, fill=exp), alpha = 0.2) +
  #geom_line(aes(x = x, y = ymean, col=exp), size = 1) +
  #facet_grid(cols = vars(fct_rev(exp))) + 
  theme_bw() + ylab("Speed (kb/min)") + xlab("S-phase progression")+
  scale_fill_manual(values = c( "#C0BDBA", "#6D90CA", "#F9A11B", 
                                "#E92724"))+
  scale_color_manual(values = c("#C0BDBA", "#6D90CA", "#F9A11B", 
                                "#E92724" ))+
  theme(legend.position = "bottom")+ylim(0.5,2.5)+coord_cartesian()

#  }

ML <- (MM|LL|OO2)+plot_layout(widths = c(2,6,2))
ML

(GG/ML)+plot_layout(heights = c(6,3))

rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
              function(x){fread(paste0(base_folder, x)
              )[, exp := gsub(".*/|_.*", "", x)
              ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
              )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
              ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4 
][str_detect(exp,"133-WT-mESC|133-FANCJKO-mESC|133-DHX36KO-mESC"), rank := abs(rank - 1)  
][rank < 0.9
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
] %>% rstatix::t_test(u ~ exp)




LL



S_model <- loess(S ~ u + s, pc_sim[u_rank < 10 | s_rank < 10])


OO2 <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
                     function(x){fread(paste0(base_folder, x)
                     )[, exp := gsub(".*/|_.*", "", x)
                     ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                     )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                     ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 5e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4 
][str_detect(exp,"133-WT-mESC|133-FANCJKO-mESC|133-DHX36KO-mESC"), rank := abs(rank - 1)  
][rank < 0.9
][, s := predict(S_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
] %>% #dplyr::filter(str_detect(exp, "133-WT-mESC|133-DHX36-FANCJKO-mESC")) %>%
  ggplot()+
  geom_boxplot(aes(x=fct_rev(exp), y=s, fill=fct_rev(exp)), outlier.shape = NA, alpha=0.25)+
  geom_jitter(aes(x=fct_rev(exp), y=s, col=fct_rev(exp)), size =0.6, alpha=0.5)+
  #geom_ribbon(aes(x = x, ymin = ymean - y95, ymax = ymean + y95, fill=exp), alpha = 0.2) +
  #geom_ribbon(aes(x = x, ymin = ymean - yse, ymax = ymean + yse, fill=exp), alpha = 0.2) +
  #geom_line(aes(x = x, y = ymean, col=exp), size = 1) +
  #facet_grid(cols = vars(fct_rev(exp))) + 
  theme_bw() + ylab("Variance (kb/min)") + xlab('S-phase Progression')+
  scale_fill_manual(values = c( "#C0BDBA", "#6D90CA", "#F9A11B", 
                                "#E92724"))+
  scale_color_manual(values = c("#C0BDBA", "#6D90CA", "#F9A11B", 
                                "#E92724" ))+
  theme(legend.position = "none")+
  #xlim(0,0.8)+
  ylim(0.15,0.45)



d1 <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
                    function(x){fread(paste0(base_folder, x)
                    )[, exp := gsub(".*/|_.*", "", x)
                    ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                    )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                    ]})
)%>% dplyr::filter(!is.na(rank))%>% 
  #dplyr::filter(str_detect(exp, "DMSO|10m"))%>%
  distinct(cell, exp)%>% mutate(plate = ceiling(cell/384)) %>% group_by(plate, exp) %>% summarize(rank=n())

d2 <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
                    function(x){fread(paste0(base_folder, x)
                    )[, exp := gsub(".*/|_.*", "", x)
                    ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                    )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                    ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2.5 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4 
][#str_detect(exp, "DMSO|24"), rank := abs(rank - 1)  
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
]%>%  #dplyr::filter(str_detect(exp, "DMSO|10m"))%>% 
  mutate(plate = ceiling(cell/384)) %>% distinct(plate, cell, exp)%>% group_by(plate, exp) %>% summarize(speed=n())


OO <- left_join(d1, d2, by=c("exp", "plate"))%>% 
  dplyr::filter(str_detect(exp, "133-WT-mESC|133-DHX36-FANCJKO-mESC")) %>%
  mutate(ratio = speed/rank, ratioEdU = rank/1520)%>%
  ggplot()+
  geom_boxplot(aes(x=fct_rev(exp),y=ratio, fill= fct_rev(exp)),alpha=0.5)+
  geom_jitter(aes(x=fct_rev(exp),y=ratio, col= fct_rev(exp)),size =3)+
  theme_bw()+theme(legend.position = "none")+ ylim(0,1)+
  scale_fill_manual(values = c( "#C0BDBA",# "#6D90CA", "#F9A11B", 
                                "#E92724"))+
  scale_color_manual(values = c("#C0BDBA",# "#6D90CA", "#F9A11B", 
                                "#E92724" ))+
  ggtitle("Ratio of RFS+/EdU+")+xlab('')+ylab('Percentage cells')


OOO <- left_join(d1, d2, by=c("exp", "plate"))%>% 
  dplyr::filter(str_detect(exp, "133-WT-mESC|133-DHX36-FANCJKO-mESC")) %>%
  mutate(ratio = speed/rank, ratioEdU = rank/376)%>%
  ggplot()+
  geom_boxplot(aes(x=fct_rev(exp),y=ratioEdU, fill= fct_rev(exp)),alpha=0.5)+
  geom_jitter(aes(x=fct_rev(exp),y=ratioEdU, col= fct_rev(exp)))+
  theme_bw()+theme(legend.position = "none")+ ylim(0,1)+
  scale_fill_manual(values = c( "#C0BDBA",# "#6D90CA", "#F9A11B", 
                                "#E92724"))+
  scale_color_manual(values = c("#C0BDBA", #"#6D90CA", "#F9A11B", 
                                "#E92724" ))+
  ggtitle("QC PASS")+xlab('')+ylab('Percentage of cells')

OO|OOO

rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
              function(x){fread(paste0(base_folder, x)
              )[, exp := gsub(".*/|_.*", "", x)
              ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
              )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
              ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4 
][str_detect(exp,"133-WT-mESC|133-FANCJKO-mESC|133-DHX36KO-mESC"), rank := abs(rank - 1)  
][rank < 0.9
][, s := predict(S_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
][] %>% rstatix::t_test(s ~ exp)

#ggpol::geom
library(patchwork)

(OO/OO2)+plot_layout(heights = c(2,6))
rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
              function(x){fread(paste0(base_folder, x)
              )[, exp := gsub(".*/|_.*", "", x)
              ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
              )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
              ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 2e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4 
][str_detect(exp,"133-WT-mESC|133-FANCJKO-mESC|133-DHX36KO-mESC"), rank := abs(rank - 1)  
][rank < 1
][, s := predict(S_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
][] %>% dplyr::filter(str_detect(exp, "133-WT-mESC|133-DHX36-FANCJKO-mESC")) %>%
  mutate(rankbin = (round(rank*400))/400,
         plate = ceiling(cell/384))%>% 
  group_by(exp)%>% 
  arrange(rank)%>% 
  mutate(cumulative = cumsum(rep(1,n()))/n())%>%
  ggplot()+geom_line(aes(x=rank, y=cumulative, col=fct_rev(exp)), size=1.5)+theme_bw()+
  ylab('Cumulative RFS+ cells')+xlab('S-phase Progression')+
  #scale_fill_manual(values = c("#6D6E71", "#007D20"))+
  #scale_color_manual(values = c("#6D6E71", "#007D20"))+coord_cartesian()+
  theme(legend.position = 'bottom')+coord_equal()+
  scale_fill_manual(values = c( "#C0BDBA",# "#6D90CA", "#F9A11B", 
                                "#E92724"))+
  scale_color_manual(values = c("#C0BDBA", #"#6D90CA", "#F9A11B", 
                                "#E92724" ))

#PP

rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
              function(x){fread(paste0(base_folder, x)
              )[, exp := gsub(".*/|_.*", "", x)
              ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
              )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
              ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4 
][str_detect(exp,"133-WT-mESC|133-FANCJKO-mESC|133-DHX36KO-mESC"), rank := abs(rank - 1)  
][rank < 1
][, s := predict(S_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
][] %>% 
  mutate(rankbin = (round(rank*100))/100)%>% 
  group_by(exp)%>% 
  arrange(rank)%>% 
  mutate(cumulative = cumsum(rep(1,n()))/n())%>%
  ggplot()+geom_line(aes(x=rank, y=cumulative, col=exp), size=1.5)+theme_bw()+
  ylab('Cumulative RFS+ cells')+xlab('S-phase Progression')+
  #scale_fill_manual(values = c("#6D6E71", "#007D20"))+
  #scale_color_manual(values = c("#6D6E71", "#007D20"))+coord_cartesian()+
  theme(legend.position = 'none')+coord_equal()+xlim(0,0.6)+ylim(0,0.6)+
  scale_color_manual(values = c("green4",  "#1F97FF", "#FF7B45" , "grey60"))
#+
#  geom_smooth(aes(x=rankbin, y=speed_diff), span=0.5, level=0.75, col="#007D20", fill="#007D20", alpha=0.1)+
# xlim(c(-.05,0.8))+
#  ylim(c(0,0.9))+
# 3 theme_bw()#ylab("Absolute speed difference (kb/min)")+xlab('S-phase progression')

PP|((PPP|PPP)/(PPP|PPP))


bbb<- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
              function(x){fread(paste0(base_folder, x)
              )[, exp := gsub(".*/|_.*", "", x)
              ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
              )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
              ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 1.5| is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4 
][str_detect(exp,"133-WT-mESC|133-FANCJKO-mESC|133-DHX36KO-mESC"), rank := abs(rank - 1)  
][rank < 0.76
][, s := predict(S_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
][] %>% group_by(exp) %>%
  mutate(total =  n()) %>%
  mutate(rankbin = (round(rank*80))/80) %>%
  #dplyr::filter(str_detect(exp, "133-WT-mESC|133-DHX36-FANCJKO-mESC")) %>%
  group_by(rankbin, exp)%>%
  mutate(N = n()) %>%
  ungroup() %>%
  group_by(exp) %>%
  mutate(maxN = max(N),
         Nnorm = N/max(N))%>% 
  distinct(exp, rankbin, Nnorm) %>%
  pivot_wider(names_from = exp, values_from = Nnorm)%>%
  mutate(log2FC = log2(`JvB-133-DHX36-FANCJKO-mESC-15E-60WO-15E-F`/  `JvB-133-WT-mESC-15E-60WO-15E-F`)) %>%
  ggplot()+geom_point(aes(x=rankbin, y=log2FC), color = "grey10")+
  geom_smooth(aes(x=rankbin, y=log2FC), fill = "grey10", color = "grey10", span =0.25, se = 0.8)+
                theme_bw()+
  scale_fill_manual(values = "#E92724")+
  scale_color_manual(values = "#E92724")
  
# PPPP
  



PPPtest <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
              function(x){fread(paste0(base_folder, x)
              )[, exp := gsub(".*/|_.*", "", x)
              ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
              )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
              ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4 
][str_detect(exp,"133-WT-mESC|133-FANCJKO-mESC|133-DHX36KO-mESC"), rank := abs(rank - 1)  
][rank < 0.9
][, s := predict(S_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
][] %>% 
  mutate(rankbin = (round(rank*100))/100)%>% 
  group_by(exp)%>% 
  arrange(rank)%>% 
  mutate(cumulative = cumsum(rep(1,n()))/n())


PPPwt <- PPPtest %>% dplyr::filter(exp == "JvB-133-WT-mESC-15E-60WO-15E-F", rank >0.05,  rank <0.5) %>% arrange(rank)
PPP_DHX36 <- PPPtest %>% dplyr::filter(exp == "JvB-133-DHX36KO-mESC-15E-60WO-15E-F", rank >0.05, rank <0.5) %>% arrange(rank)
PPP_FANCJ <- PPPtest %>% dplyr::filter(exp == "JvB-133-FANCJKO-mESC-15E-60WO-15E-F", rank >0.05, rank <0.5) %>% arrange(rank)
PPP_dbl <- PPPtest %>% dplyr::filter(exp == "JvB-133-DHX36-FANCJKO-mESC-15E-60WO-15E-F", rank >0.05, rank <0.5) %>% arrange(rank)




ks.test(PPPwt$rank, PPP_DHX36$rank, alternative = "greater", exact=T)
ks.test(PPPwt$rank, PPP_FANCJ$rank, alternative = "greater", exact=T)
ks.test(PPPwt$rank, PPP_dbl$rank, alternative = "greater", exact=T)



forks_WT <- fread("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/JvB-133-WT-mESC-15E-60WO-15E-F/JvB-133-WT-mESC-15E-60WO-15E-F_singlefork.tsv.gz"
)[state == 2
][, width_est := width + width / N
][order(center)
][, rank := as.data.table(readRDS("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/JvB-133-WT-mESC-15E-60WO-15E-F/JvB-133-WT-mESC-15E-60WO-15E-F_sphaseorder.RDS"
)$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]][, genotype := 'WT'
][width_est < 1.5e6][]%>% mutate(start = center-(0.5*width_est), end = center+(0.5*width_est))%>% select(chr, start, end, rank)%>%
  mutate(rank = abs(rank-1))%>%
  as.data.table()






G4 <- fread('~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/G4_CnT_combined_peaks_DESeq_DKO_sig_lfc_base_cutoff.bed')%>%
  mutate(V1 = gsub("chr", "", V1)) %>% select(V1, V2, V3, V5) %>% setNames(c("chr", "start", "end", "score"))%>% as.data.table()




setkeyv(forks_WT, c('chr', 'start', 'end'))
setkeyv(G4, c('chr', 'start', 'end'))


foverlaps(forks_WT, G4)%>% mutate(enriched = "yes")%>%  drop_na()%>% 
  ggplot()+geom_histogram(aes(x=rank))+
  #scale_fill_viridis_c(option="A")+
  theme_bw()


G4non <- fread('~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/G4_CnT_combined_peaks_DESeq_nonsig.bed')%>%
  mutate(V1 = gsub("chr", "", V1)) %>% select(V1, V2, V3, V5) %>% setNames(c("chr", "start", "end", "score"))%>% as.data.table()

setkeyv(G4non, c('chr', 'start', 'end'))

bind_rows(foverlaps(forks_WT, G4)%>% mutate(enriched = "yes")%>%  drop_na(),
foverlaps(forks_WT, G4non)%>% mutate(enriched = "no")%>%  drop_na())%>%
  ggplot()+geom_histogram(aes(x=rank, fill =enriched), position = position_identity(),  bins=50)+
  facet_grid(cols = vars(enriched), scales = "fixed", space = "fixed")+
  #scale_fill_viridis_c(option="A")+
  theme_bw()

bind_rows(foverlaps(forks_WT, G4)%>% mutate(enriched = "background"),
          foverlaps(forks_WT, G4)%>% mutate(enriched = "yes")%>%  drop_na(),
          foverlaps(forks_WT, G4non)%>% mutate(enriched = "no")%>%  drop_na())%>%
  ggplot()+geom_histogram(aes(x=rank, y = ..density.., fill =enriched), position = position_identity(),  bins=60)+
  facet_grid(cols = vars(enriched), scales = "fixed", space = "fixed")+
  #scale_fill_viridis_c(option="A")+
  theme_bw()


bind_rows(foverlaps(forks_WT, G4)%>% mutate(enriched = "background"),
          foverlaps(forks_WT, G4)%>% mutate(enriched = "yes")%>%  drop_na(),
          foverlaps(forks_WT, G4non)%>% mutate(enriched = "no")%>% drop_na()) %>%
  group_by(enriched)%>%
  mutate(rankbin = round(rank*50),
         total  = n())%>%
  group_by(enriched, rankbin)%>%
  summarize(n=(n()/total)) %>%
  distinct(enriched, rankbin, n) %>%
  ggplot()+
  geom_col(aes(x=rankbin*2, y = n, fill =enriched), position = position_identity(), width =2)+
  facet_grid(cols = vars(enriched), scales = "fixed", space = "fixed")+
  #scale_fill_viridis_c(option="A")+
  theme_bw()
  
  
aa <- bind_rows(foverlaps(forks_WT, G4)%>% mutate(enriched = "background"),
          foverlaps(forks_WT, G4)%>% mutate(enriched = "yes")%>%  drop_na(),
          foverlaps(forks_WT, G4non)%>% mutate(enriched = "no")%>% drop_na()) %>%
  group_by(enriched)%>%
  mutate(rankbin = round(rank*100),
         total  = n())%>%
  group_by(enriched, rankbin, total)%>%
  summarize(n=(n()/total)) %>%
  distinct(enriched, rankbin, n) %>%  ungroup()%>% select(rankbin, enriched, n) %>%
  pivot_wider(names_from = enriched, values_from = n) %>%
  mutate(norm_nonDJ_G4 = no/background,
         normDJ_G4 = yes/background,
         log2FC = (normDJ_G4/norm_nonDJ_G4))%>%
  pivot_longer(cols = c(norm_nonDJ_G4, normDJ_G4#,log2FC
                        )) %>%
  #select(norm_G4, normDJ_G4)%>%
  ggplot()+
  geom_point(aes(x=rankbin, y=log2(value), col = name))+
  geom_smooth(aes(x=rankbin, y=log2(value), col = name, fill=name), alpha=0.33)+
  theme_bw()+
  scale_fill_manual(values = c( "#C0BDBA", #"#6D90CA", "#F9A11B", 
                                "#E92724"))+
  scale_color_manual(values = c("#C0BDBA", #"#6D90CA", "#F9A11B", 
                                "#E92724" ))+
  xlab("Replication Timing")+
  ylab("Normalized occurence of G4 type (log2)")+
  theme(legend.position = "bottom")

#"#C0BDBA", "#6D90CA", "#F9A11B", "#E92724"
aaa <- bind_rows(foverlaps(forks_WT, G4)%>% mutate(enriched = "background"),
          foverlaps(forks_WT, G4)%>% mutate(enriched = "yes")%>%  drop_na(),
          foverlaps(forks_WT, G4non)%>% mutate(enriched = "no")%>% drop_na()) %>%
  group_by(enriched)%>%
  mutate(rankbin = round(rank*100),
         total  = n())%>%
  group_by(enriched, rankbin, total)%>%
  summarize(n=(n()/total)) %>%
  distinct(enriched, rankbin, n) %>%  ungroup()%>% select(rankbin, enriched, n) %>%
  pivot_wider(names_from = enriched, values_from = n) %>%
  mutate(norm_G4 = no/background,
         normDJ_G4 = yes/background,
         log2FC = (normDJ_G4/norm_G4))%>%
  pivot_longer(cols = c(log2FC
  )) %>% dplyr::filter(rankbin <76) %>%
  #select(norm_G4, normDJ_G4)%>%
  ggplot()+
  geom_point(aes(x=rankbin, y=log2(value), col = name))+
  geom_smooth(aes(x=rankbin, y=log2(value), col = name, fill=name), alpha=0.33)+
  theme_bw()+
  scale_fill_manual(values = c(  "grey60" ))+
  scale_color_manual(values = c(  "grey60"  ))  +
  xlab("Replication Timing")+
  theme(legend.position = "bottom")+
  ylab("Enrichment DJ-G4 \n log2(DJ-G4/nonDJ-G4)")

(aa|aaa)

aaaa <- bind_rows(foverlaps(forks_WT, G4)%>% mutate(enriched = "background"),
          foverlaps(forks_WT, G4)%>% mutate(enriched = "yes")%>%  drop_na(),
          foverlaps(forks_WT, G4non)%>% mutate(enriched = "no")%>% drop_na()) %>%
  group_by(enriched)%>%
  mutate(rankbin = round(rank*100),
         total  = n())%>%
  group_by(enriched, rankbin, total)%>%
  summarize(n=(n()/total)) %>%
  distinct(enriched, rankbin, n) %>%  ungroup()%>% select(rankbin, enriched, n) %>%
  pivot_wider(names_from = enriched, values_from = n) %>%
  mutate(norm_G4 = no/background,
         normDJ_G4 = yes/background,
         log2FC = (normDJ_G4/norm_G4))%>%
  pivot_longer(cols = c(log2FC
  )) %>% #dplyr::filter(rankbin <81) %>%
  #select(norm_G4, normDJ_G4)%>%
  ggplot()+
  geom_raster(aes(x=rankbin, y=1, fill = log2(value)))+
  xlab("Replication Timing")+
  theme(legend.position = "right")+
  scale_fill_viridis_c(option="B")+
  #scale_fill_distiller(palette = "RdYlBu")+
  theme_bw()

(aaaa/aa)+plot_layout(heights = c(1,6))

bbbb<-rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
                    function(x){fread(paste0(base_folder, x)
                    )[, exp := gsub(".*/|_.*", "", x)
                    ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                    )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                    ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 1.5| is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4 
][str_detect(exp,"133-WT-mESC|133-FANCJKO-mESC|133-DHX36KO-mESC"), rank := abs(rank - 1)  
][rank < 0.76
][, s := predict(S_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
][] %>% group_by(exp) %>%
  mutate(total =  n()) %>%
  mutate(rankbin = (round(rank*100))/100) %>%
  #dplyr::filter(str_detect(exp, "133-WT-mESC|133-DHX36-FANCJKO-mESC")) %>%
  group_by(rankbin, exp)%>%
  mutate(N = n()) %>%
  ungroup() %>%
  group_by(exp) %>%
  mutate(maxN = max(N),
         Nnorm = N/max(N))%>% 
  distinct(exp, rankbin, Nnorm)%>%
  ggplot()+
  geom_col(aes(x=rankbin, y=Nnorm, group= fct_rev(exp), fill=fct_rev(exp)), width=0.01, alpha=1, position = position_identity())+
  facet_grid(rows=vars(exp))+
  theme_bw()+theme(legend.position= 'none')+
  scale_fill_manual(values = c( "#C0BDBA", "#6D90CA", "#F9A11B", 
                                "#E92724"))+
  scale_color_manual(values = c("#C0BDBA", "#6D90CA", "#F9A11B", 
                                "#E92724" ))+
  xlab('Replication Timing')+
  ylab('RFS+ normalized counts')



(aaaa/bbb/bbbb)+plot_layout(heights = c(0.5,1,4))

bbb1 <-rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
                    function(x){fread(paste0(base_folder, x)
                    )[, exp := gsub(".*/|_.*", "", x)
                    ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                    )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                    ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 1.5| is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4 
][str_detect(exp,"133-WT-mESC|133-FANCJKO-mESC|133-DHX36KO-mESC"), rank := abs(rank - 1)  
][rank < 0.76
][, s := predict(S_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
][] %>% group_by(exp) %>%
  mutate(total =  n()) %>%
  mutate(rankbin = (round(rank*80))/80) %>%
  #dplyr::filter(str_detect(exp, "133-WT-mESC|133-DHX36-FANCJKO-mESC")) %>%
  group_by(rankbin, exp)%>%
  mutate(N = n()) %>%
  ungroup() %>%
  group_by(exp) %>%
  mutate(maxN = max(N),
         Nnorm = N/max(N))%>% 
  distinct(exp, rankbin, Nnorm) %>%
  pivot_wider(names_from = exp, values_from = Nnorm)%>%
  mutate(log2FC = log2(`JvB-133-FANCJKO-mESC-15E-60WO-15E-F`/  `JvB-133-WT-mESC-15E-60WO-15E-F`)) %>%
  ggplot()+geom_point(aes(x=rankbin, y=log2FC), color = "#6D90CA")+
  geom_smooth(aes(x=rankbin, y=log2FC), fill = "#6D90CA", color = "#6D90CA", span =0.25, se = 0.8)+
  theme_bw()+ggtitle("FANCJKO")+ylim(-1.5,0.3)+
  scale_fill_manual(values = "#E92724")+
  scale_color_manual(values = "#E92724")

  bbb2<- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
                function(x){fread(paste0(base_folder, x)
                )[, exp := gsub(".*/|_.*", "", x)
                ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                ]})
  )[((!diffusion) & model != "quadruple") | (diffusion)
  ][log2(est_sd / stde_mean) > 1.5| is.na(stde_mean)
  ][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
  ][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
  ][parameter == 4 
  ][str_detect(exp,"133-WT-mESC|133-FANCJKO-mESC|133-DHX36KO-mESC"), rank := abs(rank - 1)  
  ][rank < 0.76
  ][, s := predict(S_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
  ][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
  ][] %>% group_by(exp) %>%
  mutate(total =  n()) %>%
  mutate(rankbin = (round(rank*80))/80) %>%
  #dplyr::filter(str_detect(exp, "133-WT-mESC|133-DHX36-FANCJKO-mESC")) %>%
  group_by(rankbin, exp)%>%
  mutate(N = n()) %>%
  ungroup() %>%
  group_by(exp) %>%
  mutate(maxN = max(N),
         Nnorm = N/max(N))%>% 
  distinct(exp, rankbin, Nnorm) %>%
  pivot_wider(names_from = exp, values_from = Nnorm)%>%
  mutate(log2FC = log2(`JvB-133-DHX36-FANCJKO-mESC-15E-60WO-15E-F`/  `JvB-133-WT-mESC-15E-60WO-15E-F`)) %>%
  ggplot()+geom_point(aes(x=rankbin, y=log2FC), color ="#E92724")+ #"#6D90CA", "#F9A11B", 
  #  "#E92724"
  geom_smooth(aes(x=rankbin, y=log2FC), fill = "#E92724", color = "#E92724", span =0.25, se = 0.8)+
  theme_bw()+ggtitle("DKO")+ylim(-1.5,0.3)+
  scale_fill_manual(values = "#E92724")+
  scale_color_manual(values = "#E92724")

  bbb3<-rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
                function(x){fread(paste0(base_folder, x)
                )[, exp := gsub(".*/|_.*", "", x)
                ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                ]})
  )[((!diffusion) & model != "quadruple") | (diffusion)
  ][log2(est_sd / stde_mean) > 1.5| is.na(stde_mean)
  ][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
  ][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
  ][parameter == 4 
  ][str_detect(exp,"133-WT-mESC|133-FANCJKO-mESC|133-DHX36KO-mESC"), rank := abs(rank - 1)  
  ][rank < 0.76
  ][, s := predict(S_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
  ][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
  ][] %>% group_by(exp) %>%
  mutate(total =  n()) %>%
  mutate(rankbin = (round(rank*80))/80) %>%
  #dplyr::filter(str_detect(exp, "133-WT-mESC|133-DHX36-FANCJKO-mESC")) %>%
  group_by(rankbin, exp)%>%
  mutate(N = n()) %>%
  ungroup() %>%
  group_by(exp) %>%
  mutate(maxN = max(N),
         Nnorm = N/max(N))%>% 
  distinct(exp, rankbin, Nnorm) %>%
  pivot_wider(names_from = exp, values_from = Nnorm)%>%
  mutate(log2FC = log2(`JvB-133-DHX36KO-mESC-15E-60WO-15E-F`/  `JvB-133-WT-mESC-15E-60WO-15E-F`)) %>%
  ggplot()+geom_point(aes(x=rankbin, y=log2FC), color = "#F9A11B")+
  geom_smooth(aes(x=rankbin, y=log2FC), fill = "#F9A11B", color = "#F9A11B", span =0.25, se = 0.8)+
  theme_bw()+ggtitle("DHX36")+ylim(-1.5,0.3)+
  scale_fill_manual(values = "#E92724")+
  scale_color_manual(values = "#E92724")



bbb3/bbb1/bbb2



rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
              function(x){fread(paste0(base_folder, x)
              )[, exp := gsub(".*/|_.*", "", x)
              ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
              )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
              ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 1.5| is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4 
][str_detect(exp,"133-WT-mESC|133-FANCJKO-mESC|133-DHX36KO-mESC"), rank := abs(rank - 1)  
][rank < 1.1
][, s := predict(S_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
][] %>% group_by(exp) %>%
  mutate(total =  n()) %>%
  mutate(rankbin = (round(rank*25))/25,
         plate = ceiling(cell/384)) %>%
  dplyr::filter(str_detect(exp, "133-WT-mESC|133-DHX36-FANCJKO-mESC")) %>%
  group_by(rankbin, exp, plate)%>%
  mutate(N = n()) %>%
  ungroup() %>%
  group_by(exp, plate) %>%
  mutate(maxN = max(N),
         Nnorm = N/max(N)) %>%
  mutate(DJ_G4 =  if_else(rankbin < 0.2, true = "DJ-G4 enriched", false = "G4 enriched")) %>%
  group_by(DJ_G4, exp, plate) %>%
  summarize(Nnorm = median(Nnorm))%>%
  ggplot()+
  geom_boxplot(aes(x=paste(DJ_G4, exp), y=Nnorm, fill = exp), alpha=0.1, outlier.shape = NA, width=0.2)+
  geom_jitter(aes(x=paste(DJ_G4, exp), y=Nnorm, color = exp),width = 0.05, height=0.02)+
  theme_bw()+
  scale_fill_manual(values = c( "#C0BDBA",# "#6D90CA", "#F9A11B", 
                                "#E92724"))+
  scale_color_manual(values = c("#C0BDBA", #"#6D90CA", "#F9A11B", 
                                "#E92724" ))




#library(tidyverse)
poisson <- bind_rows(
  (data.table::fread('~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/JvB-133-WT-mESC-15E-60WO-15E-F/JvB-133-WT-mESC-15E-60WO-15E-F_poissontest.tsv.gz')%>% arrange(cell)%>%
  mutate(cellid = 1:n(), exp = "JvB-133-WT-mESC-15E-60WO-15E-F")),
  (data.table::fread('~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/JvB-133-DHX36KO-mESC-15E-60WO-15E-F/JvB-133-DHX36KO-mESC-15E-60WO-15E-F_poissontest.tsv.gz')%>% arrange(cell)%>%
    mutate(cellid = 1:n(), exp= "JvB-133-DHX36KO-mESC-15E-60WO-15E-F")),
  (data.table::fread('~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/JvB-133-FANCJKO-mESC-15E-60WO-15E-F/JvB-133-FANCJKO-mESC-15E-60WO-15E-F_poissontest.tsv.gz')%>% arrange(cell)%>%
    mutate(cellid = 1:n(), exp = "JvB-133-FANCJKO-mESC-15E-60WO-15E-F")),
  (data.table::fread('~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/JvB-133-DHX36-FANCJKO-mESC-15E-60WO-15E-F/JvB-133-DHX36-FANCJKO-mESC-15E-60WO-15E-F_poissontest.tsv.gz')%>% arrange(cell)%>%
    mutate(cellid = 1:n(), exp = "JvB-133-DHX36-FANCJKO-mESC-15E-60WO-15E-F"))
)%>% as.data.table()


roll_median <- function(x, windowsize){
  
  
  frollapply(c(rep(NA, windowsize), x, rep(NA, windowsize)), 
             n = windowsize * 2 + 1, 
             align = "center",
             FUN = function(y){
               xx <- y[!is.na(y)]
               
               approx((1:length(xx)) / length(xx), sort(xx), 0.5)$y}
  )[(1:length(x)) + windowsize]
}


rollfun            <- function(x, w, FUN = "se", window = 0.25){
  
  was <- ceiling(length(x) * window / 2) 
  
  rowe  <- frollsum(c(rep(NA, was), w, rep(NA, was)),  
                    n = was*2+1, 
                    align = "center",
                    na.rm = T, hasNA = T)[1:length(x) + was]
  
  roweme  <- frollsum(c(rep(NA, was), x*w, rep(NA, was)),  
                      n = was*2+1, 
                      align = "center",
                      na.rm = T, hasNA = T)[1:length(x) + was] / rowe
  
  if(FUN == "mean"){return(roweme)}
  
  if(FUN == "se"){
    return(
      sqrt(frollsum(c(rep(NA, was), (x - roweme)^2 * w, rep(NA, was)),  
                    n = was*2+1, 
                    align = "center",
                    na.rm = T)[1:length(x) + was] / rowe))}
}


cellnumbers = seq(1:381)
CN <- cellnumbers[!cellnumbers %in% 358:360]

library(data.table)
idx <- bind_rows(
  
 (fread("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/index/20230328WT ES WT_001_index.csv") %>%
    mutate(condition = "15E-60WO-15E-F", 
    treatment = "WT",
    DAPI = `*[405] 460/50`,
    exp ="JvB-133-WT-mESC-15E-60WO-15E-F",
    cell= paste("JvB-133-WT-mESC-15E-60WO-15E-F-pl1", CN, sep = "_"))),
 
 
 (fread("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/index/20230328WT ES WT_002_index.csv") %>%
    mutate(condition = "15E-60WO-15E-F", 
           treatment = "WT",
           DAPI = `*[405] 460/50`,
           exp ="JvB-133-WT-mESC-15E-60WO-15E-F",
           cell= paste("JvB-133-WT-mESC-15E-60WO-15E-F-pl2", CN, sep = "_"))),
 
 
 (fread("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/index/20230328WT ES WT_003_index.csv") %>%
    mutate(condition = "15E-60WO-15E-F", 
           treatment = "WT",
           DAPI = `*[405] 460/50`,
           exp ="JvB-133-WT-mESC-15E-60WO-15E-F",
           cell= paste("JvB-133-WT-mESC-15E-60WO-15E-F-pl3", CN, sep = "_"))),
 
 
 (fread("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/index/20230328WT ES WT_004_index.csv") %>%
    mutate(condition = "15E-60WO-15E-F", 
           treatment = "WT",
           DAPI = `*[405] 460/50`,
           exp ="JvB-133-WT-mESC-15E-60WO-15E-F",
           cell= paste("JvB-133-WT-mESC-15E-60WO-15E-F-pl4", CN, sep = "_"))),
 
 
 
 
 
 (fread("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/index/20230328WT ES FANCJ_001_index.csv") %>%
    mutate(condition = "15E-60WO-15E-F", 
           treatment = "FANCJ-KO",
           DAPI = `*[405] 460/50`,
           exp = "JvB-133-FANCJKO-mESC-15E-60WO-15E-F",
           cell= paste("JvB-133-FANCJKO-mESC-15E-60WO-15E-F-pl1", CN, sep = "_"))),
 
 (fread("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/index/20230328WT ES FANCJ_002_index.csv") %>%
    mutate(condition = "15E-60WO-15E-F", 
           treatment = "FANCJ-KO",
           DAPI = `*[405] 460/50`,
           exp = "JvB-133-FANCJKO-mESC-15E-60WO-15E-F",
           cell= paste("JvB-133-FANCJKO-mESC-15E-60WO-15E-F-pl2", CN, sep = "_"))),
 
 (fread("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/index/20230328WT ES FANCJ_003_index.csv") %>%
    mutate(condition = "15E-60WO-15E-F", 
           treatment = "FANCJ-KO",
           DAPI = `*[405] 460/50`,
           exp = "JvB-133-FANCJKO-mESC-15E-60WO-15E-F",
           cell= paste("JvB-133-FANCJKO-mESC-15E-60WO-15E-F-pl3", CN, sep = "_"))),
 
 (fread("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/index/20230328WT ES FANCJ_004_index.csv") %>%
    mutate(condition = "15E-60WO-15E-F", 
           treatment = "FANCJ-KO",
           DAPI = `*[405] 460/50`,
           exp = "JvB-133-FANCJKO-mESC-15E-60WO-15E-F",
           cell= paste("JvB-133-FANCJKO-mESC-15E-60WO-15E-F-pl4", CN, sep = "_"))),
 
 
 (fread("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/index/20230328WT ES DHX36_001_index.csv") %>%
    mutate(condition = "15E-60WO-15E-F", 
           treatment = "DHX36-KO",
           DAPI = `*[405] 460/50`,
           exp= "JvB-133-DHX36KO-mESC-15E-60WO-15E-F",
           cell= paste("JvB-133-DHX36KO-mESC-15E-60WO-15E-F-pl1", CN, sep = "_"))),
 
 (fread("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/index/20230328WT ES DHX36_002_index.csv") %>%
    mutate(condition = "15E-60WO-15E-F", 
           treatment = "DHX36-KO",
           DAPI = `*[405] 460/50`,
           exp= "JvB-133-DHX36KO-mESC-15E-60WO-15E-F",
           cell= paste("JvB-133-DHX36KO-mESC-15E-60WO-15E-F-pl2", CN, sep = "_"))),
 (fread("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/index/20230328WT ES DHX36_003_index.csv") %>%
    mutate(condition = "15E-60WO-15E-F", 
           treatment = "DHX36-KO",
           DAPI = `*[405] 460/50`,
           exp= "JvB-133-DHX36KO-mESC-15E-60WO-15E-F",
           cell= paste("JvB-133-DHX36KO-mESC-15E-60WO-15E-F-pl3", CN, sep = "_"))),
 (fread("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/index/20230328WT ES DHX36_004_index.csv") %>%
    mutate(condition = "15E-60WO-15E-F", 
           treatment = "DHX36-KO",
           DAPI = `*[405] 460/50`,
           exp= "JvB-133-DHX36KO-mESC-15E-60WO-15E-F",
           cell= paste("JvB-133-DHX36KO-mESC-15E-60WO-15E-F-pl4", CN, sep = "_"))),
 
 
 
 (fread("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/index/20230328WT ES DHX36 FANCJ_001_index.csv") %>%
    mutate(condition = "15E-60WO-15E-F", 
           treatment = "FANCJ-DHX36-KO",
           DAPI = `*[405] 460/50`,
           exp= "JvB-133-DHX36-FANCJKO-mESC-15E-60WO-15E-F",
           cell= paste("JvB-133-DHX36-FANCJKO-mESC-15E-60WO-15E-F-pl1", CN, sep = "_"))),
 
 (fread("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/index/20230328WT ES DHX36 FANCJ_002_index.csv") %>%
    mutate(condition = "15E-60WO-15E-F", 
           treatment = "FANCJ-DHX36-KO",
           DAPI = `*[405] 460/50`,
           exp= "JvB-133-DHX36-FANCJKO-mESC-15E-60WO-15E-F",
           cell= paste("JvB-133-DHX36-FANCJKO-mESC-15E-60WO-15E-F-pl2", CN, sep = "_"))),
 (fread("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/index/20230328WT ES DHX36 FANCJ_003_index.csv") %>%
    mutate(condition = "15E-60WO-15E-F", 
           treatment = "FANCJ-DHX36-KO",
           DAPI = `*[405] 460/50`,
           exp= "JvB-133-DHX36-FANCJKO-mESC-15E-60WO-15E-F",
           cell= paste("JvB-133-DHX36-FANCJKO-mESC-15E-60WO-15E-F-pl3", CN, sep = "_"))),
 (fread("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/index/20230328WT ES DHX36 FANCJ_004_index.csv") %>%
    mutate(condition = "15E-60WO-15E-F", 
           treatment = "FANCJ-DHX36-KO",
           DAPI = `*[405] 460/50`,
           exp= "JvB-133-DHX36-FANCJKO-mESC-15E-60WO-15E-F",
           cell= paste("JvB-133-DHX36-FANCJKO-mESC-15E-60WO-15E-F-pl4", CN, sep = "_")))
)%>% group_by(exp, treatment)%>% as.data.table()#mutate()


a_idx <- merge(idx, poisson, by.x = "id", by.y = "cell" )%>% mutate(EdUpos = (poisson_score > 0.1 ) & log(mean) < 2 & log(mean) > - 2.5 ) %>%
  ggplot() +
  geom_point(aes(x = log(mean), y = log(sqrt(variance) / mean), col= EdUpos),
             size = 0.3, alpha=0.4) + theme_bw()+
  # facet_wrap(vars(exp), scales = "free") +
  theme(legend.position = "none") + 
  geom_abline(slope = -0.5, size =0.25, linetype = "dashed",  color ="black", intercept = 0.1)+
  #geom_hline(yintercept = 0.1) + 
  geom_vline(xintercept = c(-2.5, 2), size =0.25, linetype = "dashed",  color ="black")+
  scale_color_manual(values = c("grey50", "blue1"))


background <- bind_rows(
  (flowCore::read.FCS(filename = "~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/index/20230328WT ES WT_001.fcs",transformation="linearize", truncate_max_range = FALSE)%>%
     {.@exprs} %>% as.tibble %>%
     mutate(condition = "15E-60WO-15E-F", 
            treatment = "WT",
            DAPI = `*[405] 460/50`,
            exp ="JvB-133-WT-mESC-15E-60WO-15E-F")),
  
  (flowCore::read.FCS(filename = "~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/index/20230328WT ES DHX36 FANCJ_001.fcs",transformation="linearize", truncate_max_range = FALSE)%>%
     {.@exprs} %>% as.tibble %>%
     mutate(condition = "15E-60WO-15E-F", 
            treatment = "FANCJ-DHX36-KO",
            DAPI = `*[405] 460/50`,
            exp= "JvB-133-DHX36-FANCJKO-mESC-15E-60WO-15E-F")),
  (flowCore::read.FCS(filename = "~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/index/20230328WT ES DHX36_001.fcs",transformation="linearize", truncate_max_range = FALSE)%>%
     {.@exprs} %>% as.tibble %>%
     mutate(condition = "15E-60WO-15E-F", 
            treatment = "DHX36-KO",
            DAPI = `*[405] 460/50`,
            exp= "JvB-133-DHX36-KO-mESC-15E-60WO-15E-F")),
  (flowCore::read.FCS(filename = "~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/index/20230328WT ES FANCJ_001.fcs",transformation="linearize", truncate_max_range = FALSE)%>%
     {.@exprs} %>% as.tibble %>%
     mutate(condition = "15E-60WO-15E-F", 
            treatment = "FANCJ-KO",
            DAPI = `*[405] 460/50`,
            exp= "JvB-133-FANCJKO-mESC-15E-60WO-15E-F")))%>% as.data.table()%>%
  mutate(ratio_up =`[405] 460/50 Area`/`*[405] 460/50`)%>% dplyr::filter(ratio_up<1, DAPI >1e4, DAPI <4e4) %>%
  group_by(exp)%>%
  slice_sample(n=1000)

 merge.data.table(idx, poisson, by.x = "cell", by.y = "cell" )%>% 
  mutate(EdUpos = (poisson_score > 0.1 ) & log(mean) < 2 & log(mean) > - 2.5 ) %>% 
  #dplyr::filter(str_detect(exp.x, "133-WT-mESC|133-DHX36-FANCJKO-mESC"))%>%
  ggplot()+
  geom_histogram(aes(x=DAPI/1e4, fill = EdUpos), position= position_identity(), alpha = 0.6, bins=100)+
  geom_histogram(data = background, aes(x=DAPI/1e4), fill = "grey40", position= position_identity(), alpha = 0.2, bins=100)+
  facet_grid(cols=vars(fct_rev(exp.x)))+
  geom_vline(xintercept = 1.75, col = "red", linetype=2, size = 0.5)+
  theme_bw()+ ylab("Number of cells") + xlab("DAPI (A.U., x10e3)")+theme(legend.position = "none") + 
  scale_fill_manual(values = c("orange3", "green3"))#+xlim(1.75,4.25)+coord_fixed(ratio =0.02500)#+  facet_grid(cols=vars(exp.x))

a_idx|b_idx


merge.data.table(idx, poisson, by.x = "cell", by.y = "cell" )%>% 
  mutate(EdUpos = (poisson_score > 0.1 ) & log(mean) < 2 & log(mean) > - 2.5 ) %>% 
  #dplyr::filter(str_detect(exp.x, "133-WT-mESC|133-DHX36-FANCJKO-mESC"))%>%
  ggplot()+
  geom_histogram(aes(x=DAPI/1e4, fill = EdUpos), position= position_identity(), alpha = 0.6, bins=100)+
  geom_histogram(data = background, aes(x=DAPI/1e4), fill = "grey40", position= position_identity(), alpha = 0.2, bins=100)+
  facet_grid(cols=vars(fct_rev(exp.x)))+
  geom_vline(xintercept = 1.75, col = "red", linetype=2, size = 0.5)+
  theme_bw()+ ylab("Number of cells") + xlab("DAPI (A.U., x10e3)")+theme(legend.position = "none") + 
  scale_fill_manual(values = c("orange3", "green3"))#+xlim(1.75,4.25)+coord_fixed(ratio =0.02500)#+  facet_grid(cols=vars(exp.x))

merge.data.table(idx, poisson, by.x = "cell", by.y = "cell" )%>% 
  mutate(EdUpos = (poisson_score > 0.1 ) & log(mean) < 2 & log(mean) > - 2.5 ) %>% 
  #dplyr::filter(str_detect(exp.x, "133-WT-mESC|133-DHX36-FANCJKO-mESC"))%>%
  ggplot()+
  geom_histogram(aes(x=DAPI/1e4, fill = EdUpos), position= position_identity(), alpha = 0.6, bins=100)+
  geom_histogram(data = background, aes(x=DAPI/1e4), fill = "grey40", position= position_identity(), alpha = 0.2, bins=100)+
  facet_grid(cols=vars(fct_rev(exp.x)))+
  geom_vline(xintercept = 1.75, col = "red", linetype=2, size = 0.5)+
  theme_bw()+ ylab("Number of cells") + xlab("DAPI (A.U., x10e3)")+theme(legend.position = "none") + 
  scale_fill_manual(values = c("orange3", "green3"))#+xlim(1.75,4.25)+coord_fixed(ratio =0.02500)#+  facet_grid(cols=vars(exp.x))



pois_idx <- merge.data.table(idx, poisson, by.x = "cell", by.y = "cell" )%>% 
  mutate(EdUpos = (poisson_score > 0.1 ) & log(mean) < 2 & log(mean) > - 2.5 )%>%
  mutate(cellid2 = paste0(exp.x ,"_", cellid))%>% as.data.table()



RFS <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
              function(x){fread(paste0(base_folder, x)
              )[, exp := gsub(".*/|_.*", "", x)
              ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
              )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
              ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2| is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4 
][str_detect(exp,"133-WT-mESC|133-FANCJKO-mESC|133-DHX36KO-mESC"), rank := abs(rank - 1)  
][
][, s := predict(S_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
][] %>% distinct(rank, cell, exp, u) %>% mutate(cellid2 = paste0(exp,"_", cell)) %>% as.data.table()

background2 <- bind_rows(
  (flowCore::read.FCS(filename = "~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/index/20230328WT ES WT_001.fcs",transformation="linearize", truncate_max_range = FALSE)%>%
     {.@exprs} %>% as.tibble %>%
     mutate(condition = "15E-60WO-15E-F", 
            treatment = "WT",
            DAPI = `*[405] 460/50`,
            exp ="JvB-133-WT-mESC-15E-60WO-15E-F")),
  
  (flowCore::read.FCS(filename = "~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/index/20230328WT ES DHX36 FANCJ_001.fcs",transformation="linearize", truncate_max_range = FALSE)%>%
     {.@exprs} %>% as.tibble %>%
     mutate(condition = "15E-60WO-15E-F", 
            treatment = "FANCJ-DHX36-KO",
            DAPI = `*[405] 460/50`,
            exp= "JvB-133-DHX36-FANCJKO-mESC-15E-60WO-15E-F")),
  (flowCore::read.FCS(filename = "~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/index/20230328WT ES DHX36_001.fcs",transformation="linearize", truncate_max_range = FALSE)%>%
     {.@exprs} %>% as.tibble %>%
     mutate(condition = "15E-60WO-15E-F", 
            treatment = "DHX36-KO",
            DAPI = `*[405] 460/50`,
            exp= "JvB-133-DHX36-KO-mESC-15E-60WO-15E-F")),
  (flowCore::read.FCS(filename = "~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/index/20230328WT ES FANCJ_001.fcs",transformation="linearize", truncate_max_range = FALSE)%>%
     {.@exprs} %>% as.tibble %>%
     mutate(condition = "15E-60WO-15E-F", 
            treatment = "FANCJ-KO",
            DAPI = `*[405] 460/50`,
            exp= "JvB-133-FANCJKO-mESC-15E-60WO-15E-F")))%>% as.data.table()%>%
  mutate(ratio_up =`[405] 460/50 Area`/`*[405] 460/50`)%>% dplyr::filter(ratio_up<1, DAPI >1e4, DAPI <4e4) %>%
  group_by(exp)%>%
  slice_sample(n=500)

left_join(pois_idx, RFS, by = "cellid2") %>%
  mutate(RFSpos = u > 0) %>%
  mutate_at(vars(RFSpos), ~replace_na(., FALSE)) %>%
  
  ggplot()+
  geom_histogram(aes(x=DAPI/1e4, fill = RFSpos), position= position_identity(), alpha = 0.6, bins=100)+
  geom_histogram(data = background2, aes(x=DAPI/1e4), fill = "grey40", position= position_identity(), alpha = 0.2, bins=100)+
  facet_grid(cols=vars(fct_rev(exp.x)))+
  geom_vline(xintercept = 1.77, col = "red", linetype=2, size = 0.5)+
  theme_bw()+ ylab("Number of cells") + xlab("DAPI (A.U., x10e3)")+theme(legend.position = "top") + 
  scale_fill_manual(values = c("orange3", "green3"))

#dt <- left_join(pois_idx, RFS, by = "cellid2")

#133-WT-mESC|133-FANCJKO-mESC|133-DHX36KO-mESC

forks_1 <-  fread("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/JvB-133-WT-mESC-15E-60WO-15E-F/JvB-133-WT-mESC-15E-60WO-15E-F_singlefork.tsv.gz"
)[state == 2
][, width_est := width + width / N
][order(center)
][, rank := as.data.table(readRDS("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/JvB-133-WT-mESC-15E-60WO-15E-F/JvB-133-WT-mESC-15E-60WO-15E-F_sphaseorder.RDS"
)$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]][, treatment := 'WT'
]%>% mutate(rank = abs(rank-1))%>% as.data.table()

forks_1[, bin := round(center / 2e5) * 2e5]
forks_1[order(rank), cumfire := cumsum(rep(1, .N)), .(chr, bin)]
forks_1[, cumfire := cumfire / uniqueN(forks_10$cell)]
forks_1[order(rank), cumreads := cumsum(N), .(chr, bin)]


forks_2 <-  fread("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/JvB-133-DHX36KO-mESC-15E-60WO-15E-F/JvB-133-DHX36KO-mESC-15E-60WO-15E-F_singlefork.tsv.gz"
)[state == 2
][, width_est := width + width / N
][order(center)
][, rank := as.data.table(readRDS("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/JvB-133-DHX36KO-mESC-15E-60WO-15E-F/JvB-133-DHX36KO-mESC-15E-60WO-15E-F_sphaseorder.RDS"
)$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]][, treatment := 'DHX36'
]%>% mutate(rank = abs(rank-1)) %>% as.data.table()
          

forks_2[, bin := round(center / 2e5) * 2e5]
forks_2[order(rank), cumfire := cumsum(rep(1, .N)), .(chr, bin)]
forks_2[, cumfire := cumfire / uniqueN(forks_10$cell)]
forks_2[order(rank), cumreads := cumsum(N), .(chr, bin)]




forks_3 <-   fread("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/JvB-133-FANCJKO-mESC-15E-60WO-15E-F/JvB-133-FANCJKO-mESC-15E-60WO-15E-F_singlefork.tsv.gz"
)[state == 2
][, width_est := width + width / N
][order(center)
][, rank := as.data.table(readRDS("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/JvB-133-FANCJKO-mESC-15E-60WO-15E-F/JvB-133-FANCJKO-mESC-15E-60WO-15E-F_sphaseorder.RDS"
)$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]][, treatment := 'FANCJ'
]%>% mutate(rank = abs(rank-1)) %>% as.data.table()

forks_3[, bin := round(center / 2e5) * 2e5]
forks_3[order(rank), cumfire := cumsum(rep(1, .N)), .(chr, bin)]
forks_3[, cumfire := cumfire / uniqueN(forks_10$cell)]
forks_3[order(rank), cumreads := cumsum(N), .(chr, bin)]


forks_4 <-   fread("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/JvB-133-DHX36-FANCJKO-mESC-15E-60WO-15E-F/JvB-133-DHX36-FANCJKO-mESC-15E-60WO-15E-F_singlefork.tsv.gz"
)[state == 2
][, width_est := width + width / N
][order(center)
][, rank := as.data.table(readRDS("~/Desktop/post-doc/Experiments/exp133-mESC-FANCJ-DHX36KO-15E-60WO-15E-F/JvB-133-DHX36-FANCJKO-mESC-15E-60WO-15E-F/JvB-133-DHX36-FANCJKO-mESC-15E-60WO-15E-F_sphaseorder.RDS"
)$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]][, treatment := 'DKO'
]%>% #mutate(rank = abs(rank-1)) %>% 
  as.data.table()


forks_4[, bin := round(center / 2e5) * 2e5]
forks_4[order(rank), cumfire := cumsum(rep(1, .N)), .(chr, bin)]
forks_4[, cumfire := cumfire / uniqueN(forks_10$cell)]
forks_4[order(rank), cumreads := cumsum(N), .(chr, bin)]


forks_plus <- bind_rows(forks_1, forks_2, forks_3, forks_4)

forks_plus[chr == 2 & bin > 4e7 & bin < 7e7 & width_est < 1e6] %>%
  mutate(rank = abs(rank-1))%>%
  ggplot() +
  geom_segment(aes(x = (center - width_est / 2) / 1e6, xend = (center + width_est / 2) / 1e6, y = rank, yend = rank), color ="grey40", size =0.33) +
  ylab("S-phase progression") +
  xlab("Genomic location [Mb]") +
  coord_cartesian(expand = F) + scale_fill_manual("grey40")+
  theme_bw()+facet_grid(rows = vars(factor(treatment, c("WT",
                                                        "DHX36",
                                                        "FANCJ",
                                                        "DKO"))))

forks_plus %>% distinct(treatment, rank) %>% group_by(treatment) %>% summarise(n=n())

forksa <- forks_plus %>% #mutate(rank = rank/max(rank))%>% 
  #dplyr::filter(rank< 0.8)%>%
  group_by(treatment, rank)%>%
  summarize(n= n())%>%
  ggplot()+geom_violin(aes(y = factor(treatment, c("WT",
                                                   "DHX36",
                                                   "FANCJ",
                                                   "DKO")),
                           x=n,
                           fill= factor(treatment, c("WT",
                                                     "DHX36",
                                                     "FANCJ",
                                                     "DKO"))), alpha =0.3)+
  geom_jitter(aes(y = factor(treatment, c("WT",
                                          "DHX36",
                                          "FANCJ",
                                          "DKO")),
                  x=n, 
                  col= factor(treatment, c("WT",
                                           "DHX36",
                                           "FANCJ",
                                           "DKO"))), alpha=0.5, size=0.75)+
  #xlim(0,8000)+ 
  xlab("DNA replication forks per cell")+ylab("")+
  theme_bw()+  
  scale_fill_manual(values = c("grey60", "#FF7B45", "magenta", "magenta4", "blue", "yellow4" ))+
  scale_color_manual(values = c( "grey60", "#FF7B45", "magenta", "magenta4", "blue", "yellow4" ))+ theme(legend.position = "none")



forksb <- forks_plus %>% #mutate(rank = rank/max(rank))%>% 
  dplyr::filter(rank< 0.9)%>%
  group_by(treatment, rank)%>%
  summarize(n= n())%>%
  ggplot()+
  geom_smooth(aes(x = rank, y=n, col= factor(treatment, c("WT",
                                                          "DHX36",
                                                          "FANCJ",
                                                          "DKO")), fill=factor(treatment, c("WT",
                                                                                            "DHX36",
                                                                                            "FANCJ",
                                                                                            "DKO"))), method = 'gam')+
  geom_point(aes(x = rank, y=n, col=factor(treatment, c("WT",
                                                        "DHX36",
                                                        "FANCJ",
                                                        "DKO"))), size=0.6, alpha = 0.6)+
  #ylim(0,8000)+ 
  xlab("S-phase Progression")+ ylab("DNA replication forks per cell")+
  theme_bw()+xlab("")+  
  theme_bw()+  
  scale_fill_manual(values = c("grey60", "#FF7B45", "magenta", "magenta4", "blue", "yellow4" ))+
  facet_wrap(vars(treatment),ncol = 2)+
  scale_color_manual(values = c( "grey60", "#FF7B45", "magenta", "magenta4", "blue", "yellow4" ))+ theme(legend.position = "none")


(forksa/forksb)+plot_layout(heights =  c(6,6))


Early_Sphase_corr <- forks_plus %>% dplyr::filter(rank < 0.25) %>%#dplyr::filter(rank < 0.6)%>% 
  group_by(chr, treatment) %>% 
  mutate(bin  = round(center / 2e5) * 2e5)%>% 
  group_by(chr, bin, treatment) %>% 
  summarize(n= n()) %>%
  group_by(treatment)%>%
  mutate(loc = paste0(chr, "-", bin))%>%
  pivot_wider(names_from = treatment, values_from = n)%>%
  mutate(WT = WT,
         DHX36 = DHX36,
         FANCJ = FANCJ,
         DKO = DKO)%>% ungroup() %>% select(-chr, -bin, -loc)%>% drop_na()%>%
  cor(method="pearson") %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname) %>% 
  #mutate(rowname = factor(rowname,levels = c("ref_early","DOT_early","ref_late","DOT_late")),
  #       name = factor(name,levels = rev(c("ref_early","DOT_early","ref_late","DOT_late")))) %>%
  ggplot(aes(x=rowname, y=name))+ 
  geom_tile(aes(fill=value),color="black", size=0.35)+
  theme_bw()+
  ylab("")+
  xlab("")+
  scale_fill_distiller(palette = "Spectral", direction = 1) +
  geom_text(aes(label=round(value,2)), size=5)+coord_fixed()+
  #guides(shape = guide_legend(override.aes = list(size = 1.5)),
  #       color = guide_legend(override.aes = list(size = 1.5))) +
  theme(legend.position = 'bottom')+ggtitle('Early S-phase RT Pearson correlation')  



Late_Sphase_corr <- forks_plus %>% dplyr::filter(rank > 0.75) %>%#dplyr::filter(rank < 0.6)%>% 
  group_by(chr, treatment) %>% 
  mutate(bin  = round(center / 2e5) * 2e5)%>% 
  group_by(chr, bin, treatment) %>% 
  summarize(n= n()) %>%
  group_by(treatment)%>%
  mutate(loc = paste0(chr, "-", bin))%>%
  pivot_wider(names_from = treatment, values_from = n)%>%
  mutate(WT = WT,
         DHX36 = DHX36,
         FANCJ = FANCJ,
         DKO = DKO)%>% ungroup() %>% select(-chr, -bin, -loc)%>% drop_na()%>%
  cor(method="pearson") %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname) %>% 
  #mutate(rowname = factor(rowname,levels = c("ref_early","DOT_early","ref_late","DOT_late")),
  #       name = factor(name,levels = rev(c("ref_early","DOT_early","ref_late","DOT_late")))) %>%
  ggplot(aes(x=rowname, y=name))+ 
  geom_tile(aes(fill=value),color="black", size=0.35)+
  theme_bw()+
  ylab("")+
  xlab("")+
  scale_fill_distiller(palette = "Spectral", direction = 1) +
  geom_text(aes(label=round(value,2)), size=5)+coord_fixed()+
  #guides(shape = guide_legend(override.aes = list(size = 1.5)),
  #       color = guide_legend(override.aes = list(size = 1.5))) +
  theme(legend.position = 'bottom')+ggtitle('Late S-phase RT pearson correlation')  



Mid_Sphase_corr <- forks_plus %>% dplyr::filter(rank >0.375 & rank < 0.625) %>%#dplyr::filter(rank < 0.6)%>% 
  group_by(chr, treatment) %>% 
  mutate(bin  = round(center / 2e5) * 2e5)%>% 
  group_by(chr, bin, treatment) %>% 
  summarize(n= n()) %>%
  group_by(treatment)%>%
  mutate(loc = paste0(chr, "-", bin))%>%
  pivot_wider(names_from = treatment, values_from = n)%>%
  mutate(WT = WT,
         DHX36 = DHX36,
         FANCJ = FANCJ,
         DKO = DKO)%>% ungroup() %>% select(-chr, -bin, -loc)%>% drop_na()%>%
  cor(method="pearson") %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname) %>% 
  #mutate(rowname = factor(rowname,levels = c("ref_early","DOT_early","ref_late","DOT_late")),
  #       name = factor(name,levels = rev(c("ref_early","DOT_early","ref_late","DOT_late")))) %>%
  ggplot(aes(x=rowname, y=name))+ 
  geom_tile(aes(fill=value),color="black", size=0.35)+
  theme_bw()+
  ylab("")+
  xlab("")+
  scale_fill_distiller(palette = "Spectral", direction = 1) +
  geom_text(aes(label=round(value,2)), size=5)+coord_fixed()+
  #guides(shape = guide_legend(override.aes = list(size = 1.5)),
  #       color = guide_legend(override.aes = list(size = 1.5))) +
  theme(legend.position = 'bottom')+ggtitle('Mid S-phase RT Pearson correlation')  

forks_plus %>%dplyr::filter(rank >0.66) %>% 
  distinct(cell, .keep_all = T)%>%
  group_by(treatment)%>%
  summarise(n=n())

Late <- forks_plus %>% dplyr::filter(rank >0.66) %>% 
  group_by(chr, treatment) %>% 
  mutate(bin  = round(center / 1.5e5) * 1.5e5)%>% 
  group_by(chr, bin, treatment) %>% 
  summarize(n= n()) %>%
  group_by(treatment)%>%
  mutate(loc = paste0(chr, "-", bin))%>% 
  pivot_wider(names_from = treatment, values_from = n)%>%
  mutate(WT = WT/331,
         DHX36 = DHX36/261,
         FANCJ = FANCJ/178,
         DKO = DKO/114)%>% select(loc, DHX36, DKO, FANCJ, WT)%>%
  setNames(c("loc",  "Late_DHX36", "Late_DKO", "Late_FANCJ", "Late_WT"))%>% as.data.table()

forks_plus %>% dplyr::filter(rank < 0.66 & rank > 0.33 ) %>% 
  distinct(cell, .keep_all = T)%>%
  group_by(treatment)%>%
  summarise(n=n())

Mid <- forks_plus %>% dplyr::filter(rank < 0.66 & rank > 0.33 ) %>% 
  group_by(chr, treatment) %>% 
  mutate(bin  = round(center / 1.5e5) * 1.5e5)%>% 
  group_by(chr, bin, treatment) %>% 
  summarize(n= n()) %>%
  group_by(treatment)%>%
  mutate(loc = paste0(chr, "-", bin))%>% 
  pivot_wider(names_from = treatment, values_from = n)%>%
  mutate(WT = WT/322,
         DHX36 = DHX36/254,
         FANCJ = FANCJ/178,
         DKO = DKO/104)%>% select(loc, DHX36, DKO, FANCJ, WT)%>%
  setNames(c("loc",  "Mid_DHX36", "Mid_DKO", "Mid_FANCJ", "Mid_WT"))%>% as.data.table()

forks_plus %>% dplyr::filter(rank < 0.333 )%>%
  distinct(cell, .keep_all = T)%>%
  group_by(treatment)%>%
  summarise(n=n())


Early <- forks_plus %>% dplyr::filter(rank < 0.333 ) %>% 
  group_by(chr, treatment) %>% 
  mutate(bin  = round(center / 1.5e5) * 1.5e5)%>% 
  group_by(chr, bin, treatment) %>% 
  summarize(n= n()) %>%
  group_by(treatment)%>%
  mutate(loc = paste0(chr, "-", bin))%>% 
  pivot_wider(names_from = treatment, values_from = n)%>%
  mutate(WT = WT/326,
         DHX36 = DHX36/251,
         FANCJ = FANCJ/179,
         DKO = DKO/107)%>% select(loc, DHX36, DKO, FANCJ, WT)%>%
  setNames(c("loc",  "Early_DHX36", "Early_DKO", "Early_FANCJ", "Early_WT"))%>% as.data.table()


EM <- merge(Early, Mid, by = "loc" )

All <- merge(EM, Late, by = "loc"  )

corrplot <- All %>% ungroup() %>% select(-loc)%>% drop_na()%>%
  cor(method="pearson") %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname) %>% 
  #mutate(rowname = factor(rowname,levels = c("ref_early","DOT_early","ref_late","DOT_late")),
  #       name = factor(name,levels = rev(c("ref_early","DOT_early","ref_late","DOT_late")))) %>%
  ggplot(aes(x=factor(rowname, c("Early_WT",
                                 "Early_DHX36",
                                 "Early_FANCJ",
                                 "Early_DKO",
                                 "Mid_WT",
                                 "Mid_DHX36",
                                 "Mid_FANCJ",
                                 "Mid_DKO",
                                 "Late_WT",
                                 "Late_DHX36",
                                 "Late_FANCJ",
                                 "Late_DKO"
                                 
                                 )),
             y= factor(name, c(
               "Late_DKO", "Late_FANCJ",   "Late_DHX36","Late_WT",
                            
               "Mid_DKO", "Mid_FANCJ",     "Mid_DHX36", "Mid_WT",        
                             
               "Early_DKO","Early_FANCJ", "Early_DHX36","Early_WT"        
                              ))))+ 
  geom_tile(aes(fill=value),color="black", size=0.35)+
  theme_bw()+
  ylab("")+
  xlab("")+
  scale_fill_distiller(palette = "RdYlBu", direction = -1, name = "Pearson \nCorrelation") +
  geom_text(aes(label=round(value,2)), size=3
            )+coord_fixed()+
  #guides(shape = guide_legend(override.aes = list(size = 1.5)),
  #       color = guide_legend(override.aes = list(size = 1.5))) +
  theme(legend.position = 'right')+ggtitle('Replication Timing') 

corrplot

#"#C0BDBA", "#6D90CA", "#F9A11B", "#E92724"


# quadruple 
# 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
# em nm nv nm nv um pe pn pn pu

rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
              function(x){fread(paste0(base_folder, x)
              )[, exp := gsub(".*/|_.*", "", x)
              ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
              )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
              ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][#log2(est_sd / stde_mean) > 2| is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 6e5 & est_sd < 10e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][str_detect(exp,"133-WT-mESC|133-FANCJKO-mESC|133-DHX36KO-mESC"), rank := abs(rank - 1)  
] %>% dplyr::filter(parameter %in% c(7,9)) %>% select(exp, cell,  rank, parameter, Estimate) %>%
  pivot_wider(names_from =parameter, values_from = Estimate) %>%
  setNames(c("exp", "cell", "rank","pe", "pn")) %>%
  mutate(ratio = pn/pe) %>%
  ggplot()+
  geom_jitter(aes(x=rank, y=ratio, col= fct_rev(exp)), alpha=0.8, size=1)+
  geom_smooth(aes(x=rank, y=ratio, col = fct_rev(exp), fill= fct_rev(exp)),  alpha = 0.33)+
  theme_bw()+theme(legend.position= 'none')+
  scale_fill_manual(values = c( "#C0BDBA", "#6D90CA", "#F9A11B", 
                                "#E92724"))+
  scale_color_manual(values = c("#C0BDBA", "#6D90CA", "#F9A11B", 
                                "#E92724" ))+
  xlab('Replication Timing')+
  facet_grid(cols=vars(fct_rev(exp)))+
  ylab('PN/PE')




rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
              function(x){fread(paste0(base_folder, x)
              )[, exp := gsub(".*/|_.*", "", x)
              ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
              )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
              ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][str_detect(exp,"133-WT-mESC|133-FANCJKO-mESC|133-DHX36KO-mESC"), rank := abs(rank - 1)  
][] %>% dplyr::filter(parameter %in% c(7,9)) %>% select(exp, cell,  rank, parameter, Estimate) %>%
  pivot_wider(names_from =parameter, values_from = Estimate) %>%
  setNames(c("exp", "cell", "rank","pe", "pn")) %>%
  mutate(ratio = pn/pe) %>%
  ggplot()+
  geom_jitter(aes(x=rank, y=ratio, col= fct_rev(exp)), alpha=0.8, size=1)+
  geom_smooth(aes(x=rank, y=ratio, col = fct_rev(exp), fill= fct_rev(exp)),  alpha = 0.33)+
  theme_bw()+theme(legend.position= 'none')+
  scale_fill_manual(values = c( "#C0BDBA", "#6D90CA", "#F9A11B", 
                                "#E92724"))+
  scale_color_manual(values = c("#C0BDBA", "#6D90CA", "#F9A11B", 
                                "#E92724" ))+
  xlab('Replication Timing')+
  facet_grid(cols=vars(fct_rev(exp)))+
  ylab('PN/PE')




