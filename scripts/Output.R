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
  regular.tandem.posture()
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
  d.trap.obs$Tandem.enter
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
  odir <- "img/"
  
  fname <- paste(odir, "speedComparison.txt", sep = "")
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
  fname <- paste(odir, "speedComparison.pdf", sep = "")
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
  
  # Surrogate
  ggplot(df.plot.relative.pos[df.plot.relative.pos$treat=="surrogate",], aes(rpy,rpx)) + 
    #stat_bin_hex(bins=100) +
    stat_density_2d(aes(fill = stat(level)), geom="polygon", contour = TRUE, bins=7, col=1)+    scale_fill_viridis() +
    scale_x_continuous(expand = c(0, 0.05), limits = c(-plot.range, plot.range)) +
    scale_y_continuous(expand = c(0, 0.05), limits = c(-plot.range, plot.range)) +
    theme_bw() +
    coord_fixed() +
    xlab("Distance left-right (body length)") +
    ylab("Distance back-front (body length)") +
    facet_grid(~direction)
  ggsave(filename = paste0("img/RelativeDensity-Surrogate.pdf"),
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
  df.temp <- subset(df.plot.relative.pos, dis < 2 & treat != "surrogate" & direction == "FtoM")
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
  #df.temp <- subset(df.temp, treat!="surrogate")
  
  
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
  
  ks.test(df.temp[df.temp$treat=="sticky","anglediff"], 
               df.temp[df.temp$treat=="surrogate","anglediff"])
  
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
  # NOTE: the results relating to PCA was removed during revision.
  load("data/df_pca.rda")
  ggplot(df.pca.res[df.pca.res$swap<2,],aes(x=PC1, y=PC2)) + 
    stat_density_2d(geom = "polygon",  aes(alpha = ..level.., 
                    fill = as.factor(swap)))+
    scale_fill_viridis(discrete = T, direction = 1, end=0.5)+
    geom_point(data=df.pca.res[df.pca.res$swap==2,], aes(x=PC1, y=PC2))+
    coord_fixed(ylim = c(-1,1), xlim = c(-2,2), expand = F) +
    theme(aspect.ratio = 1)+
    theme_classic()
  ggsave(file.path("img/", paste0("PCA.pdf")),
         width=7, height=10)  
  
  g1 <- ggplot(df.pca.res[df.pca.res$swap<2,],aes(x=PC1, fill=as.factor(swap))) + 
    geom_histogram(position  = "identity", alpha=0.4)+
    scale_fill_viridis(discrete = T, direction = 1, end=0.5)+
    theme_classic()
  
  g2 <- ggplot(df.pca.res[df.pca.res$swap<2,],aes(x=PC2, fill=as.factor(swap))) + 
    geom_histogram(position  = "identity", alpha=0.4)+
    scale_fill_viridis(discrete = T, direction = 1, end=0.5)+
    theme_classic()
  
  g3 <- ggplot(df.pca.res[df.pca.res$swap<2,],aes(x=PC3, fill=as.factor(swap))) + 
    geom_histogram(position  = "identity", alpha=0.4)+
    scale_fill_viridis(discrete = T, direction = 1, end=0.5)+
    theme_classic()
  
  g4 <- ggplot(df.pca.res[df.pca.res$swap<2,],aes(x=PC4, fill=as.factor(swap))) + 
    geom_histogram(position  = "identity", alpha=0.4)+
    scale_fill_viridis(discrete = T, direction = 1, end=0.5)+
    theme_classic()
  g <- gridExtra::grid.arrange(g1,g2,g3,g4, ncol=1)
  ggsave(file.path("img/", paste0("PCA_each.pdf")),
         width=5, height=6, plot = g)  
  
  # representative posture
  if(F){
    plotall <- function(dftemp){
      #par(mfrow=c(5,5), pin=c(2,2))
      expand = 10
      for(i in 1:dim(dftemp)[1]){
        dftempplot = subset(df.all, name == dftemp[i,1] & time == dftemp[i,2])
        xrange = mean(as.numeric(dftempplot[,c(3,5,7,9,11,13)]))
        ytange = mean(as.numeric(dftempplot[,c(4,6,8,10,12,14)]))
        plot(0, pch="none", xlim=c(xrange-expand, xrange+expand), 
             ylim=c(ytange-expand, ytange+expand), col="#58c8acff",)
        title(paste(round(c(dftemp[i,"PC1"], dftemp[i,"PC2"]),2), 
                    collapse = ", "), line = -2)
        draw.arrow = function(x1,y1, x2, y2){ 
          arrows(dftempplot[,x1], dftempplot[,y1],
                 dftempplot[,x2], dftempplot[,y2], length=0)}
        draw.arrow(3,4,5,6)
        draw.arrow(5,6,7,8)
        draw.arrow(9,10,11,12)
        draw.arrow(11,12,13,14)
        
        #draw.arrow(3,4,9,10)
        #draw.arrow(3,4,11,12)
        #draw.arrow(3,4,13,14)

        #draw.arrow(5,6,9,10)
        #draw.arrow(5,6,11,12)
        #draw.arrow(5,6,13,14)
        
        #draw.arrow(7,8,9,10)
        #draw.arrow(7,8,11,12)
        #draw.arrow(7,8,13,14)
        
        #draw.arrow(5,6,7,8)
        draw.arrow(9,10,11,12)
        draw.arrow(11,12,13,14)
        
        points(dftempplot[,5:6], col="#2c7dd3ff", pch=19)
        points(dftempplot[,3:4], col="#7009d3ff", pch=19)
        points(dftempplot[,7:8], col="#58c8acff", pch=19)
        points(dftempplot[,9:10], col="#8fb972ff", pch=19)
        points(dftempplot[,11:12], col="#b86d3bff", pch=19)
        points(dftempplot[,13:14], col="#d51315ff", pch=19)
        
      }
    }
    
    set.seed(2)
    dftemp = df.pca.res[df.pca.res$PC1 < -1 &
                          df.pca.res$PC2 < -0.3 &
                          df.pca.res$swap == 0,]
    dftemp = dftemp[sample(1:dim(dftemp)[1],25),]
    plotall(dftemp)
     
    dftemp = df.pca.res[df.pca.res$PC1 < -0.8 &
                          df.pca.res$PC2 > -0.5 &
                          df.pca.res$PC1 > -1.2 &
                          df.pca.res$PC2 < 0.3 &
                          df.pca.res$swap == 0,]
    dftemp = dftemp[sample(1:dim(dftemp)[1],25),]
    plotall(dftemp)
    
    dftemp = df.pca.res[df.pca.res$PC1 < -0.4 &
                          df.pca.res$PC1 > -0.8 &
                          df.pca.res$swap == 0,]
    dftemp = dftemp[sample(1:dim(dftemp)[1],25),]
    plotall(dftemp)
    
    dftemp = df.pca.res[df.pca.res$PC1 > -0.02 &
                          df.pca.res$PC1 < 0.02 &
                          df.pca.res$swap == 0,]
    dftemp = dftemp[sample(1:dim(dftemp)[1],25),]
    plotall(dftemp)
    
    # a
    dftemp = df.pca.res[round(df.pca.res$PC1,2) == -1.60 & 
                          round(df.pca.res$PC2,2) == -0.58,]
    plotall(dftemp)
    
    # b
    dftemp = df.pca.res[round(df.pca.res$PC1,2) == -1.14 &
                        round(df.pca.res$PC2,2) == -0.32,]
    plotall(dftemp)
    
    # c (-0.61, -0.01)
    dftemp = df.pca.res[round(df.pca.res$PC1,2) == -0.61 &
                          round(df.pca.res$PC2,2) == -0.01,]
    plotall(dftemp)
    
    # d
    dftemp = df.pca.res[round(df.pca.res$PC1,2) == 0 &
                          round(df.pca.res$PC2,2) == 0.52,]
    plotall(dftemp)
    
    # e
    dftemp = df.pca.res[round(df.pca.res$PC1,2) == -0.01 &
                          round(df.pca.res$PC2,2) == 0.13,]
    plotall(dftemp)
    
    d.example = data.frame(PC1 = c(-1.60, -1.14, -0.61, 0, -0.01),
               PC2 = c(-0.58, -0.32, -0.01, 0.52, 0.13))
    ggplot(df.pca.res[df.pca.res$swap<2,],aes(x=PC1, y=PC2)) + 
      stat_density_2d(geom = "polygon",  aes(alpha = ..level.., 
                                             fill = as.factor(swap)))+
      scale_fill_viridis(discrete = T, direction = 1, end=0.5)+
      geom_point(data=df.pca.res[df.pca.res$swap==2,], aes(x=PC1, y=PC2))+
      coord_fixed(ylim = c(-1,1), xlim = c(-2,2), expand = F) +
      theme(aspect.ratio = 1)+
      theme_classic()+
      geom_point(data = d.example, aes(x = PC1, y = PC2), pch=2)
    
    
  }
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
plot.LDA <- function(){
  load("data/df_pca.rda")
  
  # Linear Discriminant Analysis 
  train.data = df.pca.combined[df.pca.combined$swap < 2, 
                               c("fhead_mtip", "ftip_mhead", 
                                 "fpron_mtip", "ftip_mpron",
                                 "fhead_mpron", "fpron_mhead", "swap")]
  (Z<- lda(swap~ .,data=train.data))
  
  p <- predict(Z, train.data)
  mean(p$class==train.data$swap)
  
  df.lda <- data.frame(
    train.data[,1:6],
    datasets = train.data$swap,
    lda_class = p$class,
    p$x,
    correct = (p$class == train.data$swap),
    female_prob = p$posterior[,1]
  )
  
  max(subset(df.lda, female_prob > 0.8)$LD1)
  max(subset(df.lda, female_prob > 0.6)$LD1)
  max(subset(df.lda, female_prob > 0.4)$LD1)
  max(subset(df.lda, female_prob > 0.2)$LD1)
  
  df.lda.original = df.lda[df.lda$datasets == 0,]
  sum(df.lda.original$female_prob > 0.4 & df.lda.original$female_prob < 0.6) / dim(df.lda.original)
  sum(df.lda.original$female_prob < 0.4) / dim(df.lda.original)
  sum(df.lda.original$female_prob > 0.6) / dim(df.lda.original)
  
  sum(df.lda$female_prob < 0.4) / dim(df.lda)
  sum(df.lda$female_prob > 0.6) / dim(df.lda)
  
  test.data = df.pca.combined[df.pca.combined$swap == 2,
                              c("fhead_mtip", "ftip_mhead",
                                "fpron_mtip", "ftip_mpron", 
                                "fhead_mpron", "fpron_mhead", "swap")]
  
  Y<-predict(Z,test.data)
  table(test.data[,1],Y$class)
  print("Linear Discriminant Analysis ")
  print(Y)
  
  ggplot(df.lda, aes(x=female_prob, col=as.factor(datasets), fill=as.factor(datasets)))+
    geom_step(stat="bin", binwidth=0.01, alpha=0.8, position = "identity")+
    scale_color_viridis(discrete = T, end = 0.5)+
    scale_fill_viridis(discrete = T, end=0.5) +
    coord_cartesian(xlim = c(0,1), expand = T) +
    theme_classic()+
    #facet_grid(datasets ~ .)+
    geom_vline(data=data.frame(Y$posterior), aes(xintercept = X1), color = "red")+
    theme(legend.position = "bottom")
  ggsave(file.path("img/", paste0("LDA6.pdf")), width=7, height=4)  
  
  # representative posture
  #a
  df.lda["14634",]
  df.all["14634",]
  
  df.lda["15134",]
  df.all["15134",]
  plotall(df.pca.res["15134",])
  
  #b
  df.lda["251113",]
  
  # c
  df.lda["26869",]
  plotall(df.pca.res["26869",])
  
  # d
  df.lda["17011",]
  plotall(df.pca.res["17011",])
  
  # e
  df.lda["3816101",]
  plotall(df.pca.res["3816101",])
  
  # f
  df.lda["2716141",]
  plotall(df.pca.res["2716141",])
  
  # LDA with only head-tip data
  train.data = df.pca.combined[df.pca.combined$swap < 2, 
                               c("fhead_mtip", "ftip_mhead", "swap")]
  (Z<- lda(swap~ .,data=train.data))
  sum.table = table(train.data[,1],predict(Z)$class)
  (sum.table[1,2]*2) / ((sum.table[1,2]*2)+(sum.table[1,1]*2))
  
  test.data = df.pca.combined[df.pca.combined$swap == 2, 
                              c("fhead_mtip", "ftip_mhead", "swap")]
  Y<-predict(Z,test.data)
  table(test.data[,1],Y$class)
  print("Linear Discriminant Analysis ")
  print(Y)
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
regular.tandem.posture <- function(){
  
  load("data/Cf-control-sleap/rda/df_sleap_dis.rda")
  load("data/df_pca.rda")
  
  df_sleap_dis = df_sleap_dis[df_sleap_dis$fhead_mhead < 1.5, ]
  
  test.data = df_sleap_dis[df_sleap_dis$swap == 3, 
                           c("fhead_mtip", "ftip_mhead",
                             "fpron_mtip", "ftip_mpron", 
                             "fhead_mpron", "fpron_mhead", "swap")]
  Y<-predict(Z,test.data)
  table(test.data[,1],Y$class)
  table(Y$class == 0)
  
  dd <- rbind(df.pca.res[,c(6,7,9,11,12,13,18)],
              df_sleap_dis[,c(6,7,9,11,12,13,18)])
  pl2 = prcomp(dd[,1:6], scale= F)

  ddf <- cbind(dd, pl2$x[,1:4])

  ggplot() + 
    stat_density_2d(data = ddf[ddf$swap<2,],
                    aes(x=PC1, y=PC2, fill = as.factor(swap), 
                        alpha = ..level..), geom = "polygon")+
    scale_fill_viridis(discrete = T, direction = 1, end=1)+
    geom_point(data=df.pca.res[df.pca.res$swap==2,], aes(x=PC1, y=PC2))+
    coord_fixed(xlim=c(-2,2), ylim=c(-1,1))+
    theme(aspect.ratio = 1)+
    theme_classic()+
    stat_density_2d(data = ddf[ddf$swap>2,],
                    aes(x=PC1, y=PC2, fill = as.factor(swap),
                        alpha = ..level..),
                    geom = "polygon")
  ggsave(file.path("img/PCA_w_regular.pdf"), width=5, height=10) 
  
  ## LDA only with regular tandem data
  train.data = df_sleap_dis[c("fhead_mtip", "ftip_mhead", 
                              "fpron_mtip", "ftip_mpron",
                              "fhead_mpron", "fpron_mhead", "swap")]
  
  (Z<- lda(swap~ .,data=train.data))
  
  p <- predict(Z, train.data)
  mean(p$class==train.data$swap)
  
  df.lda <- data.frame(
    train.data[,1:6],
    datasets = train.data$swap,
    lda_class = p$class,
    p$x,
    correct = (p$class == train.data$swap),
    female_prob = p$posterior[,1]
  )
  
  
  test.data = df.pca.combined[df.pca.combined$swap == 2,
                              c("fhead_mtip", "ftip_mhead",
                                "fpron_mtip", "ftip_mpron", 
                                "fhead_mpron", "fpron_mhead", "swap")]
  Y<-predict(Z,test.data)
  table(test.data[,1],Y$class)
  print("Linear Discriminant Analysis ")
  print(Y)
  
  
  ggplot(df.lda, aes(x=LD1, fill=as.factor(correct), col=as.factor(datasets)))+
    geom_histogram(binwidth=0.05, alpha=0.5)+
    scale_color_viridis(discrete=T, end = 0.5)+
    scale_fill_viridis(discrete = T, option = "A", direction = -1) +
    coord_cartesian(xlim = c(-8,8), expand = T) +
    theme_classic()+
    facet_grid(datasets ~ .)+
    theme(legend.position = "bottom")
  ggsave(file.path("img/", paste0("LDA-regular.pdf")), width=7, height=4)  
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
head.tip.dis <- function(){
  load("data/df_pca.rda")
  load("data/Cf-control-sleap/rda/df_sleap_dis.rda")
  
  ## trapped tandem
  ggplot(df.pca.combined[df.pca.combined$swap < 2 & df.pca.combined$fhead_mhead < 1.5,], 
         aes(x=ftip_mhead, fill=as.factor(swap)))+
    geom_density(alpha=0.4)+
    scale_fill_viridis(discrete=T, end=0.5)+
    theme_classic()+
    geom_vline(xintercept = df.pca.combined[df.pca.combined$swap == 2,"ftip_mhead"])+
    geom_vline(xintercept = df.pca.combined[df.pca.combined$swap == 2,"fhead_mtip"])+
    xlim(c(0,2.5))+
    theme(legend.position = "none")
  ggsave(file.path("img/", paste0("Head-Tip_dis_trap.pdf")),
         width=4, height=3)  
  
  dftemp = df.pca.combined[df.pca.combined$swap == 0 & df.pca.combined$fhead_mhead < 1.5,]
  t.test(dftemp$fhead_mtip, dftemp$ftip_mhead, paired=T)
  
  sum(dftemp$fhead_mtip < dftemp$ftip_mhead) / dim(dftemp)[1]
  
  ggplot(df.pca.combined[df.pca.combined$swap < 2 & df.pca.combined$fhead_mhead < 1.5,], 
         aes(x=fhead_mpron, fill=as.factor(swap)))+
    geom_density(alpha=0.4)+
    scale_fill_viridis(discrete=T, end=0.5)+
    theme_classic()+
    geom_vline(xintercept = df.pca.combined[df.pca.combined$swap == 2,"fhead_mpron"])+
    geom_vline(xintercept = df.pca.combined[df.pca.combined$swap == 2,"fpron_mhead"])+
    xlim(c(0,2.5))+
    theme(legend.position = "none")
  
  ggplot(df.pca.combined[df.pca.combined$swap < 2 & df.pca.combined$fhead_mhead < 1.5,], 
         aes(x=fpron_mtip, fill=as.factor(swap)))+
    geom_density(alpha=0.4)+
    scale_fill_viridis(discrete=T, end=0.5)+
    theme_classic()+
    geom_vline(xintercept = df.pca.combined[df.pca.combined$swap == 2,"fpron_mtip"])+
    geom_vline(xintercept = df.pca.combined[df.pca.combined$swap == 2,"ftip_mpron"])+
    xlim(c(0,2.5))+
    theme(legend.position = "none")
  
  ## regular tandem
  # focus on the time when the distance between female head and male head < 1.5 body length
  ggplot(df_sleap_dis[df_sleap_dis$fhead_mhead < 1.5,], 
         aes(x=ftip_mhead, fill=as.factor(swap)))+
    geom_density(alpha=0.4)+
    scale_fill_viridis(discrete=T, end=0.5)+
    theme_classic()+
    xlim(c(0,2.5))+
    theme(legend.position = "none")
  ggsave(file.path("img/", paste0("Head-Tip_dis_regular.pdf")),
         width=4, height=3)  
  
  dftemp = df_sleap_dis[df_sleap_dis$swap == 3 & df_sleap_dis$fhead_mhead < 1.5,]
  sum(dftemp$fhead_mtip > dftemp$ftip_mhead) / dim(dftemp)[1]
  
  ggplot(df_sleap_dis[df_sleap_dis$fhead_mhead < 1.5,], 
         aes(x=ftip_mpron, fill=as.factor(swap)))+
    geom_density(alpha=0.4)+
    scale_fill_viridis(discrete=T, end=0.5)+
    theme_classic()+
    xlim(c(0,2.5))+
    theme(legend.position = "none")
  
  ggplot(df_sleap_dis[df_sleap_dis$fhead_mhead < 1.5,], 
         aes(x=fpron_mhead, fill=as.factor(swap)))+
    geom_density(alpha=0.4)+
    scale_fill_viridis(discrete=T, end=0.5)+
    theme_classic()+
    xlim(c(0,2.5))+
    theme(legend.position = "none")
}
#------------------------------------------------------------------------------#
  
  