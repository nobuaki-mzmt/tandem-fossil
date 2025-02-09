## Data analysis for C.formosanus tandem movements on sticky trap
## N. Mizumoto

#------------------------------------------------------------------------------#
# This script displays the results 
# The processed data will be used, so run Preprocess.R before this.
# All the results will be stored at img/
#------------------------------------------------------------------------------#

rm(list = ls())
outputAll()

#---------------------------------------------------------------------------#
{
  library(ggplot2)
  library(viridis)
  
  library(MASS)
  
  library(exactRankTests)
  library("survminer")
  library(survival)
  library(coxme)
  library(car)
  
  library(extrafont)
  #font_import(pattern="PT") # for running for the first time
  loadfonts()

  ## constants
  pdfHeight  <- 7 
  pdfWidth   <- 10
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
outputAll <- function(){
  plot.trajectories()
  show.escapeEvents()
  compare.speed()
  plot.relative.pos()
  plot.interindividual.distance()
  plot.angle.diff()
  Compare.traveled.dis()
  plot.PCA()
  head.tip.dis()
}
#---------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
plot.trajectories <- function(){
  load("data/df.pos.rda")
  df <- subset(df.pos, treat=="sticky")
  ggplot(df) +
    geom_path(aes(x=fx, y=fy, col=viridis(2)[1])) +
    geom_path(aes(x=mx, y=my, col=viridis(2)[2])) +
    ylab("y (body length)") +
    xlab("x (body length)") +
    theme_bw()+
    theme(aspect.ratio = 1, legend.position = "none") +
    theme(strip.text = element_text(size = 8, margin = margin()))+
    theme(strip.background = element_rect(colour="#00000000", fill="#00000000"))+
    facet_wrap(.~name, ncol=4) 
  ggsave(file.path("img/", paste0("All_trajectories", ".png")),
       width=7, height=10)  
  
  df <- subset(df.pos, treat=="sticky")
  print("frames of each video clip")
  print(tapply(df$time, df$name, max))
  
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
show.escapeEvents<-function(){

  d.trap.obs <- data.frame(fread("data/SticklyTrapObservationNote.csv", header=T))
  d.trap.obs <- d.trap.obs[!is.na(d.trap.obs$Trap.10min),]
  
  # Survival plot
  {
    # reorganize data for plot
    {
      l = dim(d.trap.obs)[1]
      f = f.cens = m = m.cens = rep(0, l)
      for(i in 1:l){
        d.temp <- d.trap.obs[i,]
        
        # female
        if       (d.temp$Trap.10min == 0 || d.temp$Trap.10min == "male"){
          f[i] = 0 ; f.cens[i] = 1;
        } else if(d.temp$Trap.20min == 0 || d.temp$Trap.20min == "male"){
          f[i] = 10; f.cens[i] = 1;
        } else if(is.na(d.temp$Trap.30min)){
          f[i] = 20; f.cens[i] = 0;
        } else if(d.temp$Trap.30min == 0 || d.temp$Trap.10min == "male"){
          f[i] = 20; f.cens[i] = 1;
        } else {
          f[i] = 30; f.cens[i] = 0;
        }
        
        # male
        if       (d.temp$Trap.10min == 0 || d.temp$Trap.10min == "female"){
          m[i] = 0 ; m.cens[i] = 1;
        } else if(d.temp$Trap.20min == 0 || d.temp$Trap.20min == "female"){
          m[i] = 10; m.cens[i] = 1;
        } else if(is.na(d.temp$Trap.30min)){
          m[i] = 20; m.cens[i] = 0;
        } else if(d.temp$Trap.30min == 0 || d.temp$Trap.10min == "female"){
          m[i] = 20; m.cens[i] = 1;
        } else {
          m[i] = 30; m.cens[i] = 0;
        }
      }
      d.surv <- rbind(
        data.frame(
          d.trap.obs[,1:4],
          sex = "Female",
          duration = f,
          cens = f.cens
        ),
        data.frame(
          d.trap.obs[,1:4],
          sex = "male",
          duration = m,
          cens = m.cens
        )
      )
    }
    
    df<-survfit(Surv(duration,cens)~sex, type = "kaplan-meier", data=d.surv)
    ggsurvplot(fit = df, data = d.surv,
               pval = F, pval.method = TRUE,
               risk.table = F, conf.int = FALSE,
               ncensor.plot = FALSE, size = 1, linetype = 1:3,
               xlab="Time (min)", ggtheme = theme_bw()  + theme(aspect.ratio = 0.75))
  }
  
  # Output results as text file
  {
    fname <- paste("img/countEscapeEvents.txt", sep = "")
    sink(fname)
    {
      cat("\n----10 min----")
      print(table(d.trap.obs$Trap.10min))
      
      cat("\n----20 min----")
      print(table(d.trap.obs$Trap.20min))
      
      cat("\n----30 min----")
      print(table(d.trap.obs$Trap.30min))
    }
    
    
    cat("\n----Survival analysis----\n")
    cat("Cox mixed effect model\n")
    cat("> m <- coxme(Surv(duration, cens) ~ sex + (1|ID), data = d.surv)\n")
    m <- coxme(Surv(duration, cens) ~ sex + (1|ID), data = d.surv)
    cat("> summary(m)\n")
    summary(m)
    cat("> Anova(m)\n")
    res = Anova(m)
    cat("Chisq =", res$Chisq, "Df =", res$Df, "P =", res$`Pr(>Chisq)`, "\n")
    sink()
  }
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
compare.speed <- function(){
  load("data/dfspeed.rda")
  fname <- paste("img/speedComparison.txt", sep = "")
  sink(fname)
  
  cat("--Speed comparison between treatments--\n")
  cat("*pool sexes\n")
  cat("Mean speed\n")
  print(tapply(df.speed$speed, df.speed[,c("treat")], mean))
  cat("SD speed\n")
  print(tapply(df.speed$speed, df.speed[,c("treat")], sd))
  res <- wilcox.test(df.speed[df.speed$treat=="sticky", "speed"],
                     df.speed[df.speed$treat=="control", "speed"],
                     paired = F, alternative = "two.sided")
  print(res)
  
  
  cat("\n-- Comparison between sexes on Sticky trap\n")
  res <- wilcox.test(df.speed[df.speed$treat=="sticky" & df.speed$sex=="female", "speed"],
              df.speed[df.speed$treat=="sticky" & df.speed$sex=="male", "speed"],
              paired = TRUE, alternative = "two.sided")
  print(res)
  sink()
  
  ggpaired(df.speed[df.speed$treat=="sticky",], x = "sex", y = "speed",
           color = "sex",
           line.color = "gray", line.size = 0.4,
           palette = "viridis",
           ylab="Speed (bodylength/sec)")+
    stat_compare_means(paired = TRUE)
  fname <- paste("img/speedComparison.pdf", sep = "")
  ggsave(fname, height = pdfHeight, width = pdfWidth)
  
  
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
plot.relative.pos <- function(){
  plot.range = 1.5
  load("data/dplot_relativepos.rda")
  load("data/df_fossil.rda")
  df.relative.fossil <- data.frame(
    direction = c("FtoM", "MtoF"),
    rx = c(df.center.fossil[1,2] - df.center.fossil[2,2],
           df.center.fossil[2,2] - df.center.fossil[1,2]),
    ry = c(df.center.fossil[1,3] - df.center.fossil[2,3],
           df.center.fossil[2,3] - df.center.fossil[1,3])
  )
  
  # Sticky surface
  ggplot(df.plot.relative.pos[df.plot.relative.pos$treat=="sticky",],
         aes(rpy,rpx)) + 
    stat_density_2d(aes(fill = stat(level)), geom="polygon", contour = TRUE, bins=7, col=1)+
    scale_fill_viridis() +
    scale_x_continuous(expand = c(0, 0.05), limits = c(-plot.range, plot.range)) +
    scale_y_continuous(expand = c(0, 0.05), limits = c(-plot.range, plot.range)) +
    theme_bw() +
    coord_fixed() +
    #scale_color_viridis() +
    xlab("Distance left-right (body length)") +
    ylab("Distance back-front (body length)") +
    facet_grid(~direction)+
    geom_point(data = df.relative.fossil, aes(x = rx, y = ry, col="red"))
  
  ggsave(filename = paste0("img/RelativeDensity-Sticky.pdf"),
         width=6, height = 4, family="PT Sans")
  
  
  # Control
  ggplot(df.plot.relative.pos[df.plot.relative.pos$treat=="control",], aes(rpy,rpx)) + 
    #stat_bin_hex(bins=50) +
    #geom_pointdensity(adjust=0.05, size=0.5)+
    #stat_density_2d_filled()+
    stat_density_2d(aes(fill = stat(level)), geom="polygon", contour = TRUE, bins=7, col=1)+    scale_fill_viridis() +
    scale_x_continuous(expand = c(0, 0.05), limits = c(-plot.range, plot.range)) +
    scale_y_continuous(expand = c(0, 0.05), limits = c(-plot.range, plot.range)) +
    theme_bw() +
    coord_fixed() +
    #scale_color_viridis() +
    xlab("Distance left-right (body length)") +
    ylab("Distance back-front (body length)") +
    facet_grid(~direction)
  ggsave(filename = paste0("img/RelativeDensity-Control.pdf"),
         width=6, height = 4, family="PT Sans")
  
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
plot.interindividual.distance <- function(){
  load("data/dplot_relativepos.rda")
  load("data/df_fossil.rda")
  
  fossil_dis = sqrt(diff(df.center.fossil[,2])^2+diff(df.center.fossil[,3])^2)
  
  ## Distance
  ggplot(df.plot.relative.pos, aes(x=dis,fill=treat)) + 
    geom_density(alpha=0.3, draw_baseline = T) +
    scale_fill_viridis(discrete = T, option = "A") +
    xlim(0,2) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    xlab("Distance (body length)")+
    ylab("Density")+
    geom_vline(xintercept = fossil_dis)
  
  ggsave(filename = paste0("img/RelativeDistance.pdf"),
         width=4, height = 4, family="PT Sans")
  
  ## Compare the distance between control and sticky
  ## use the data when dis < 2 bl to remove any separation events
  df.temp <- subset(df.plot.relative.pos, dis < 2 & direction == "FtoM")
  df <- na.omit(
    rbind(data.frame(treat = "control", dis = tapply(df.temp$dis, df.temp[,c("id", "treat")], mean)[,1]),
    data.frame(treat = "sticky", dis = tapply(df.temp$dis, df.temp[,c("id", "treat")], mean)[,2])))
  r <- wilcox.exact(dis~treat, data=df, paired=F)
  cat("---Compare the distance between control and sticky---\n")
  cat("---Exact Wilcoxon rank sum test---\n")
  cat(paste("W=",r$statistic, "P=", r$p.value, "\n"))
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
plot.angle.diff <- function(){
  load("data/dplot_relativepos.rda")
  df.temp <- subset(df.plot.relative.pos, direction=="FtoM")
  anglediff <- subset(df.plot.relative.pos, direction=="FtoM")$angle -
    subset(df.plot.relative.pos, direction=="MtoF")$angle
  anglediff[anglediff<0] = anglediff[anglediff<0] + 2*pi
  df.temp$anglediff = anglediff

  ## Angle diff
  ggplot(df.temp, aes(x=anglediff,fill=treat)) + 
    geom_histogram(alpha=0.3, draw_baseline = T, position = "identity", bins=30) +
    scale_fill_viridis(discrete = T, option = "A") +
    xlim(0,2*pi)  +
    theme_bw() +
    theme(aspect.ratio = 0.75) +
    xlab("Angle difference (rad)")+
    ylab("Density")
  ggsave(filename = paste0("img/AngleDiff.pdf"),
         width=4, height = 4, family="PT Sans")
  
  ## Compare the distance between control and sticky
  ## use the data when dis < 2 bl to remove any separation events
  df.temp <- subset(df.temp, dis < 2)
  df.temp$anglediff[df.temp$anglediff > pi] = df.temp$anglediff[df.temp$anglediff > pi] - (2*pi)
  df.temp$anglediff = abs(df.temp$anglediff)

  ks.test(df.temp[df.temp$treat=="sticky","anglediff"], 
          df.temp[df.temp$treat=="control","anglediff"])
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
Compare.traveled.dis <- function(){
  load("data/df.pos.rda")
  
  sink("img/TraveledDistance.txt")
  
  # comparison of traveled distance
  print("comparison of traveled distance")
  fdis <- c(NA, sqrt(diff(df.pos$fx)^2 + diff(df.pos$fy)^2))
  mdis <- c(NA, sqrt(diff(df.pos$mx)^2 + diff(df.pos$my)^2))
  fdis[df.pos$time==0] <- NA
  mdis[df.pos$time==0] <- NA
  fdis.ind <- tapply(fdis, df.pos$id, mean, na.rm=T)
  mdis.ind <- tapply(mdis, df.pos$id, mean, na.rm=T)
  print(t.test(fdis.ind, mdis.ind, paird=T))
  
  # comparison of traveled distance when two are within 2 body length
  print("comparison of traveled distance when two are within 2 body length")
  rx <- df.pos$mx - df.pos$fx
  ry <- df.pos$my - df.pos$fy
  dis = sqrt(rx^2+ry^2)
  
  df <- rbind(
    data.frame(id=df.pos$id, sex="female", dis,  traveldis=fdis),
    data.frame(id=df.pos$id, sex = "male", dis,  traveldis=mdis)
  )
  
  df2 <- df[df$dis < 2,]
  dfsum <- tapply(df2$traveldis, df2[,1:2], sum, na.rm=T)
  print(t.test(dfsum[,1], dfsum[,2], paired=T))
  
  sink()
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
plot.PCA <- function(){
  load("data/df_pca.rda")

  ggplot(df.pca.res[3:dim(df.pca.res)[1]-2,],aes(x=PC1, y=PC2)) + 
    stat_density_2d(geom = "polygon",  aes(alpha = ..level.., 
                                           fill = as.factor(relative)))+
    scale_fill_viridis(discrete = T, direction = 1, end=0.5)+
    geom_point(data=df.pca.res[(dim(df.pca.res)[1]-1):(dim(df.pca.res)[1]),], 
               aes(x=PC1, y=PC2))+
    coord_fixed(ylim = c(-1, 1), xlim = c(-1.2,1.), expand = F) +
    theme(aspect.ratio = 1)+
    theme_classic()
  ggsave(file.path("img/", paste0("PCA.pdf")),
         width=5, height=4)  
  
  train.data = df.pca.res[, c("head", "pron", "tip", "relative")]
  test.data = df.fossil[,c("head", "pron", "tip", "relative")]
  
  train.data$relative = train.data$relative == "relative_leader"
  test.data$relative = test.data$relative == "relative_leader"
  
  lr.model <- glm(relative ~ ., data=train.data, family = binomial)
  summary(lr.model)
  probs <- predict(lr.model, test.data, type = 'response')
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
head.tip.dis <- function(){
  load("data/df_pca.rda")
  load("data/Cf-control-sleap/rda/df_sleap_dis.rda")
  
  ## trapped tandem
  ggplot(df.pca.res[df.pca.res$data_set == "Cf",], 
         aes(x=tip, fill=as.factor(relative)))+
    geom_density(alpha=0.4)+
    scale_fill_viridis(discrete=T, end=0.5)+
    theme_classic()+
    geom_vline(xintercept = df.pca.res[df.pca.res$data_set == "fossil", "tip"] )+
    xlim(c(0,2.5))+
    theme(legend.position = "none")
  ggsave(file.path("img/", paste0("Head-Tip_dis_trap.pdf")),
         width=4, height=3)  
  
  dftemp = df.pca.res[df.pca.res$data_set == "Cf",]
  dftemp = data.frame(
    dftemp[dftemp$relative=="relative_leader",1:2],
    fhead_mtip = dftemp[dftemp$relative=="relative_leader","tip"],
    ftip_mhead = dftemp[dftemp$relative=="relative_follower","tip"])
  t.test(dftemp$fhead_mtip, dftemp$ftip_mhead, paired=T)
  
  sum(dftemp$fhead_mtip > dftemp$ftip_mhead) / dim(dftemp)[1]
  
  ## regular tandem
  # focus on the time when the distance between female head and male head < 1.5 body length
  ggplot(df_sleap_dis[df_sleap_dis$fhead_mhead < 1.5,])+
    geom_density(aes(x=ftip_mhead), fill=viridis(3)[1], alpha=0.4)+
    geom_density(aes(x=fhead_mtip), fill=viridis(3)[2], alpha=0.4)+
    theme_classic()+
    xlim(c(0,2.5))+
    theme(legend.position = "none")
  ggsave(file.path("img/", paste0("Head-Tip_dis_regular.pdf")),
         width=4, height=3)  
  
  dftemp = df_sleap_dis[df_sleap_dis$fhead_mhead < 1.5,]
  sum(dftemp$fhead_mtip > dftemp$ftip_mhead) / dim(dftemp)[1]
}
#------------------------------------------------------------------------------#
  


