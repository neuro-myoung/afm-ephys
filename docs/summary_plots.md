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
dat$construct <- factor(dat$construct,levels=  c('mp1','mp2','mtraak','mtrek1','mscl','osca12','tmem63a','tmem63b','trpa1','trpv4','pkd2l1','yfp'))

names = c('mPiezo1','mPiezo2','mTRAAK','mTREK1','Mscl','Osca1.2','TMEM63a','TMEM63b','TRPA1','TRPV4','PKD2L1','YFP')
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
```

## Plot Data

The following code will create the final plot of, in this case,
Threshold as a function of the expressed construct. The plot is a
boxplot with overlayed individual datapoints.

``` r
library(pals)

ggplot(dat, aes(x=construct, y=thresh,fill=construct)) +
  geom_boxplot(alpha=0.4, lwd=0.25)+
  geom_beeswarm(shape=21, size = 1.5, alpha = 0.75,lwd=0.25)+
  scale_x_discrete(labels = names)+
  scale_fill_manual(values = glasbey(12))+
  xlab("Construct") +
  ylab("Threshold (pA)")+
  ylim(c(0,175))+
  theme_paper() +
  theme(axis.text.x = element_text(angle=45,hjust=1))
```

![](C:/Users/HAL/afm-ephys/docs/summary_plots_files/figure-gfm/pressure-1.png)<!-- -->

``` r
ggsave('agg_threshold.pdf', plot=last_plot(), dpi=300, width=5, height=3, units= "in", dev=cairo_pdf)
```

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
