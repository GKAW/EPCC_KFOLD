library(data.table)
library(plyr)
library(reshape2)
library(ggplot2)

FONT_SIZE <- 18

dataFile <- "../data/cirrus.csv" 

frame <- read.csv(dataFile,header=FALSE,col.names=c("sample","basepairs","simulation time","memory","runtime"),sep=",",stringsAsFactors = FALSE)

# bytes -> gigabytes
frame$memory <- frame$memory/1000000000

# microseconds -> milliseconds
frame$simulation.time <- frame$simulation.time/1000

# prune the unfinished runs
frame <- frame[ which(frame$simulation.time < 5), ]
frame <-frame[!(frame$sample=="stmv"),]

# split timing and memory data
time_frame <- subset(frame, select=-c(memory))
memory_frame <- subset(frame, select=-c(runtime))

timing_plot <- ggplot(data=time_frame, mapping= aes(x=as.factor(simulation.time), y=runtime,  color=as.factor(basepairs))) + 
    theme_bw() +
    geom_line(aes(group=basepairs))+ 
    geom_point(size=4,shape=19) +
    scale_color_manual(values=c("#ca0020","#0571b0","#00FF00"),name="Number of Basepairs")+  #,"#00FF00"
    #scale_x_discrete(name="Simulation Time (ms)", limits=c("1.0","2.5","3.0","4.5","7.0","8.5","10.0"),expand=c(.05,0))+
    scale_y_continuous(breaks = round(seq(0, max(frame$runtime), by = 5),1))+
    ylab("Runtime (s)") +
    xlab("Simulation Time (ms)")+
    theme(
        strip.background=element_rect(fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(size=FONT_SIZE+5),
        axis.title.x=element_text(size=FONT_SIZE+5, vjust=1),
        axis.text.y=element_text(size=FONT_SIZE+5, hjust=1),
        axis.title.y=element_text(size=FONT_SIZE+5),
        strip.text=element_text(size=FONT_SIZE+5),
        strip.text.x = element_text(margin = margin(2,2,2,2, "mm")),
        legend.position="top",
        legend.text=element_text(size=FONT_SIZE+5),
        legend.title=element_text(size=FONT_SIZE+5)
        )

ggsave("../figures/kfold_runtime_by_basepair_cirrus.pdf", plot = timing_plot,width = 400, height = 300, units="mm", pointsize=25, useDingbats=FALSE)

memory_plot <- ggplot(data=memory_frame, mapping= aes(x=as.factor(simulation.time), y=memory,  color=as.factor(basepairs))) + 
    theme_bw() +
    geom_line(aes(group=basepairs))+ 
    geom_point(size=4,shape=19) +
    scale_color_manual(values=c("#ca0020","#0571b0","#00FF00"),name="Number of Basepairs")+  #,"#00FF00"
    #scale_x_discrete(name="Simulation Time (ms)", limits=c("1.0","2.5","3.0","4.5","7.0","8.5","10.0"),expand=c(.05,0))+
    scale_y_continuous(breaks = round(seq(min(frame$memory), max(frame$memory), by = .5),1))+
    ylab("Memory (GB)") +
    xlab("Simulation Time (ms)")+
    theme(
        strip.background=element_rect(fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(size=FONT_SIZE+5),
        axis.title.x=element_text(size=FONT_SIZE+5, vjust=1),
        axis.text.y=element_text(size=FONT_SIZE+5, hjust=1),
        axis.title.y=element_text(size=FONT_SIZE+5),
        strip.text=element_text(size=FONT_SIZE+5),
        strip.text.x = element_text(margin = margin(2,2,2,2, "mm")),
        legend.position="top",
        legend.text=element_text(size=FONT_SIZE+5),
        legend.title=element_text(size=FONT_SIZE+5)
        )

ggsave("../figures/kfold_memory_by_basepair_cirrus.pdf", plot = memory_plot,width = 400, height = 300, units="mm", pointsize=25, useDingbats=FALSE)
