Summary Plots
================

This markdown document goes over how summary plots are constructed from
a table containing aggregate data.

## Load File

The code below will load the aggregate data file and create an ordered
list of factors from the constructs present.It also creates a list of
names associated with each channel of interest that will later be
mappepd onto these factors for plot aesthetics.

``` r
library(tidyverse)

dat <- read.csv('agg_bychannel.csv', header=TRUE)
dat$construct <- factor(dat$construct,levels=  c('mp1','mp2','mtraak','mtrek1','mscl','osca12',
                                                 'tmem63a','tmem63b','trpa1','trpv4','pkd2l1','yfp'),
                        ordered = T)

names = c('mPiezo1','mPiezo2','mTraak','mTrek1','Mscl','Osca1.2','Tmem63a','Tmem63b','TrpA1','TrpV4','Pkd2L1','YFP')

dat_triage <- subset(dat, construct %in% c('mtraak','mtrek1') | seal >= 0.3)
dat_signal <- subset(dat_triage, peaki >= 150)
```

### Subset of Table:

| uniqueID           |     date | construct | cell | protocol | sweep | velocity | kcant | dkcant | osm |
| :----------------- | -------: | :-------- | ---: | :------- | ----: | -------: | ----: | -----: | --: |
| 20190917-1-scan-80 | 20190917 | mtraak    |    1 | scan-80  |     5 |       40 |  0.82 |   0.05 | 332 |
| 20190918-1-scan-80 | 20190918 | mtraak    |    1 | scan-80  |     7 |       40 |  0.86 |   0.05 |   0 |
| 20190918-2-scan-80 | 20190918 | mtraak    |    2 | scan-80  |     9 |       40 |  0.86 |   0.05 | 318 |
| 20190918-3-scan-80 | 20190918 | mtraak    |    3 | scan-80  |    11 |       40 |  0.75 |   0.04 | 321 |
| 20191001-1-scan-80 | 20191001 | mtraak    |    1 | scan-80  |     6 |       40 |  1.03 |   0.06 | 316 |
| 20191001-2-scan-80 | 20191001 | mtraak    |    2 | scan-80  |     6 |       40 |  1.03 |   0.06 | 318 |

## Prepare plot themes and aesthetics

The following code will load the plotting libraries and font libraries
to be used in creating the final plot. It also creates a custom theme
function with starting plot parameters for publication quality figures.

``` r
library(ggplot2)
library(ggbeeswarm)
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
library(extrafont)
loadfonts(device = "win")

theme_paper <- function(base_size=10,base_family="Arial") {
  library(grid)
  (theme_bw(base_size = base_size, base_family = base_family)+
      theme(text=element_text(color="black"),
            axis.title=element_text(size = 10),
            axis.text=element_text(size = 8, color = "black"),
            legend.position = "none",
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid=element_blank()
      ))
}



counts <- data.frame(aggregate(wthresh ~ construct,dat_signal,length),
                     aggregate(wthresh ~ construct, dat_signal, max)[2])
colnames(counts)<- c('construct','n','max')

noSig <- data.frame(c(aggregate(wthresh ~ construct, dat_triage, 'length')))
colnames(noSig) <- c('construct','nnosig')

fullcounts <- merge(counts, noSig, by = 'construct', all = TRUE)
fullcounts$n[is.na(fullcounts$n)] <- 0
fullcounts$max[is.na(fullcounts$max)] <- 0
```

## Plot Data

The following code will create the final plot of, in this case,
Threshold as a function of the expressed construct. The plot is a
boxplot with overlayed individual datapoints.

``` r
library(pals)

p1 <- ggplot(dat_signal, aes(x=construct, y=wthresh,fill=construct)) +
  geom_boxplot(alpha=0.4, lwd=0.25)+
  geom_beeswarm(shape=21, size = 1.5, alpha = 0.75,lwd=0.25)+
  scale_x_discrete(breaks=c('mp1','mp2','mtraak','mtrek1','mscl','osca12','tmem63a','tmem63b',
                            'trpa1','trpv4','pkd2l1','yfp'), labels = names)+
  geom_text(data = fullcounts, aes(x=construct, y= max + 0.12*1200, 
                                   label = paste('(',n,'/',nnosig,')', sep = "")),size=2)+
  scale_fill_manual(values = glasbey(12))+
  scale_y_continuous(breaks=c(0,200,400,600,800,1000,1200))+
  xlab("Construct") +
  ylab("Work Threshold (fJ)")+
  theme_paper() +
  theme(axis.text.x = element_text(angle=45,hjust=1))

p1
```

![](C:/Users/HAL/afm-ephys/docs/summary_plots_files/figure-gfm/pressure-1.png)<!-- -->

``` r
ggsave('agg_wthresh.pdf', plot=last_plot(), dpi=300, width=5, height=3, units= "in", dev=cairo_pdf)
```
