This markdown document goes over how summary plots are constructed from
a table containing aggregate data.

Load File
---------

The code below will load the aggregate data file and create an ordered
list of factors from the constructs present.It also creates a list of
names associated with each channel of interest that will later be
mappepd onto these factors for plot aesthetics.

    library(tidyverse)

    dat <- read.csv('agg_bychannel.csv', header=TRUE)
    dat$construct <- factor(dat$construct,levels=  c('mp1','mp2','mtraak','mtrek1','mscl','osca12','tmem63a','tmem63b','trpa1','trpv4','pkd2l1','yfp'))

    names = c('mPiezo1','mPiezo2','mTRAAK','mTREK1','Mscl','Osca1.2','TMEM63a','TMEM63b','TRPA1','TRPV4','PKD2L1','YFP')

### Subset of Table:

<table>
<thead>
<tr class="header">
<th style="text-align: left;">uniqueID</th>
<th style="text-align: right;">date</th>
<th style="text-align: left;">construct</th>
<th style="text-align: right;">cell</th>
<th style="text-align: left;">protocol</th>
<th style="text-align: right;">sweep</th>
<th style="text-align: right;">velocity</th>
<th style="text-align: right;">kcant</th>
<th style="text-align: right;">dkcant</th>
<th style="text-align: right;">osm</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">20190917-1-scan-80</td>
<td style="text-align: right;">20190917</td>
<td style="text-align: left;">mtraak</td>
<td style="text-align: right;">1</td>
<td style="text-align: left;">scan-80</td>
<td style="text-align: right;">5</td>
<td style="text-align: right;">40</td>
<td style="text-align: right;">0.82</td>
<td style="text-align: right;">0.05</td>
<td style="text-align: right;">332</td>
</tr>
<tr class="even">
<td style="text-align: left;">20190918-1-scan-80</td>
<td style="text-align: right;">20190918</td>
<td style="text-align: left;">mtraak</td>
<td style="text-align: right;">1</td>
<td style="text-align: left;">scan-80</td>
<td style="text-align: right;">7</td>
<td style="text-align: right;">40</td>
<td style="text-align: right;">0.86</td>
<td style="text-align: right;">0.05</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">20190918-2-scan-80</td>
<td style="text-align: right;">20190918</td>
<td style="text-align: left;">mtraak</td>
<td style="text-align: right;">2</td>
<td style="text-align: left;">scan-80</td>
<td style="text-align: right;">9</td>
<td style="text-align: right;">40</td>
<td style="text-align: right;">0.86</td>
<td style="text-align: right;">0.05</td>
<td style="text-align: right;">318</td>
</tr>
<tr class="even">
<td style="text-align: left;">20190918-3-scan-80</td>
<td style="text-align: right;">20190918</td>
<td style="text-align: left;">mtraak</td>
<td style="text-align: right;">3</td>
<td style="text-align: left;">scan-80</td>
<td style="text-align: right;">11</td>
<td style="text-align: right;">40</td>
<td style="text-align: right;">0.75</td>
<td style="text-align: right;">0.04</td>
<td style="text-align: right;">321</td>
</tr>
<tr class="odd">
<td style="text-align: left;">20191001-1-scan-80</td>
<td style="text-align: right;">20191001</td>
<td style="text-align: left;">mtraak</td>
<td style="text-align: right;">1</td>
<td style="text-align: left;">scan-80</td>
<td style="text-align: right;">6</td>
<td style="text-align: right;">40</td>
<td style="text-align: right;">1.03</td>
<td style="text-align: right;">0.06</td>
<td style="text-align: right;">316</td>
</tr>
<tr class="even">
<td style="text-align: left;">20191001-2-scan-80</td>
<td style="text-align: right;">20191001</td>
<td style="text-align: left;">mtraak</td>
<td style="text-align: right;">2</td>
<td style="text-align: left;">scan-80</td>
<td style="text-align: right;">6</td>
<td style="text-align: right;">40</td>
<td style="text-align: right;">1.03</td>
<td style="text-align: right;">0.06</td>
<td style="text-align: right;">318</td>
</tr>
</tbody>
</table>

Prepare plot themes and aesthetics
----------------------------------

The following code will load the plotting libraries and font libraries
to be used in creating the final plot. It also creates a custom theme
function with starting plot parameters for publication quality figures.

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

Plot Data
---------

The following code will create the final plot of, in this case,
Threshold as a function of the expressed construct. The plot is a
boxplot with overlayed individual datapoints.

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

![](C:/Users/HAL/afm-ephys/docs/summary_plots_files/figure-markdown_strict/pressure-1.png)

    ggsave('agg_threshold.pdf', plot=last_plot(), dpi=300, width=5, height=3, units= "in", dev=cairo_pdf)

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
