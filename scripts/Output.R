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
  
  # Sticky surface
  ggplot(df.plot.relative.pos[df.plot.relative.pos$treat=="sticky",],
         aes(rpy,rpx)) + 
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
  
  ## Distance
  ggplot(df.plot.relative.pos, aes(x=dis,fill=treat)) + 
    geom_density(alpha=0.3, draw_baseline = T) +
    scale_fill_viridis(discrete = T, option = "A") +
    xlim(0,2) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    xlab("Distance (body length)")+
    ylab("Density")
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
  df <- na.omit(
    rbind(data.frame(treat = "control", anglediff = tapply(df.temp$anglediff, df.temp[,c("id", "treat")], mean)[,1]),
          data.frame(treat = "sticky", anglediff = tapply(df.temp$anglediff, df.temp[,c("id", "treat")], mean)[,2])))
  r <- wilcox.exact(anglediff~treat, data=df, paired=F)
  cat("---Compare the distance between control and sticky---\n")
  cat("---Exact Wilcoxon rank sum test---\n")
  cat(paste("W=",r$statistic, "P=", r$p.value, "\n"))
  
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
  load("data/df_pca.rda")
  ggplot(df.pca.res[df.pca.res$swap<2,],aes(x=PC1, y=PC2)) + 
    stat_density_2d(geom = "polygon",  aes(alpha = ..level.., 
                    fill = as.factor(swap)))+
    scale_fill_viridis(discrete = T, direction = 1, end=0.5)+
    geom_point(data=df.pca.res[df.pca.res$swap==2,], aes(x=PC1, y=PC2))+
    scale_y_continuous(limits = c(-2,2))+
    scale_x_continuous(limits = c(-2,2))+
    coord_fixed()+
    theme(aspect.ratio = 1)+
    theme_bw()
  
  ggsave(file.path("img/", paste0("PCA.pdf")),
         width=5, height=10)  
  
  g1 <- ggplot(df.pca.res[df.pca.res$swap<2,],aes(x=PC1, fill=as.factor(swap))) + 
    geom_histogram(position  = "identity", alpha=0.4)+
    scale_fill_viridis(discrete = T, direction = 1, end=0.5)+
    theme_bw()
  
  g2 <- ggplot(df.pca.res[df.pca.res$swap<2,],aes(x=PC2, fill=as.factor(swap))) + 
    geom_histogram(position  = "identity", alpha=0.4)+
    scale_fill_viridis(discrete = T, direction = 1, end=0.5)+
    theme_bw()
  
  g3 <- ggplot(df.pca.res[df.pca.res$swap<2,],aes(x=PC3, fill=as.factor(swap))) + 
    geom_histogram(position  = "identity", alpha=0.4)+
    scale_fill_viridis(discrete = T, direction = 1, end=0.5)+
    theme_bw()
  
  g4 <- ggplot(df.pca.res[df.pca.res$swap<2,],aes(x=PC4, fill=as.factor(swap))) + 
    geom_histogram(position  = "identity", alpha=0.4)+
    scale_fill_viridis(discrete = T, direction = 1, end=0.5)+
    theme_bw()
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
        
        draw.arrow(3,4,9,10)
        draw.arrow(3,4,11,12)
        draw.arrow(3,4,13,14)

        draw.arrow(5,6,9,10)
        draw.arrow(5,6,11,12)
        draw.arrow(5,6,13,14)
        
        draw.arrow(7,8,9,10)
        draw.arrow(7,8,11,12)
        draw.arrow(7,8,13,14)
        
        draw.arrow(5,6,7,8)
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
    dftemp = df.pca.res[df.pca.res$PC1 < -1.2 &
                          df.pca.res$PC2 < -0.5 &
                          df.pca.res$swap == 0,]
    dftemp = dftemp[sample(1:dim(dftemp)[1],25),]
    plotall(dftemp)
    
    dftemp = df.pca.res[df.pca.res$PC1 < -0.4 &
                          df.pca.res$PC2 > -0.5 &
                          df.pca.res$PC1 > -1.2 &
                          df.pca.res$PC2 < 0.3 &
                          df.pca.res$swap == 0,]
    dftemp = dftemp[sample(1:dim(dftemp)[1],25),]
    plotall(dftemp)
    
    dftemp = df.pca.res[df.pca.res$PC1 > -0.4 &
                          df.pca.res$PC2 > 0.3 &
                          df.pca.res$PC1 < 0.4 &
                          df.pca.res$swap == 0,]
    dftemp = dftemp[sample(1:dim(dftemp)[1],25),]
    plotall(dftemp)
  }
    
  
  # Linear Discriminant Analysis 
  train.data = df.pca.res[df.pca.res$swap < 2, 18:22]
  (Z<- lda(swap~ .,data=train.data))
  sum.table = table(train.data[,1],predict(Z)$class)
  (sum.table[1,2]*2) / ((sum.table[1,2]*2)+(sum.table[1,1]*2))
  
  test.data = df.pca.res[df.pca.res$swap == 2, 18:22]
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
  
  
  dd <- rbind(df.pca.res[,3:18],df_sleap_dis[,3:18])
  pl2 = prcomp(dd[,1:15], scale= F)

  ddf <- cbind(dd, pl2$x[,1:4])

  ggplot() + 
    stat_density_2d(data = ddf[ddf$swap<2,],aes(x=PC1, y=PC2, fill = "black", alpha = ..level..), geom = "polygon")+
    scale_fill_viridis(discrete = T, direction = 1, end=1)+
    geom_point(data=df.pca.res[df.pca.res$swap==2,], aes(x=PC1, y=PC2))+
    scale_y_continuous(limits = c(-2,2))+
    scale_x_continuous(limits = c(-2,2))+
    coord_fixed()+
    theme(aspect.ratio = 1)+
    theme_bw()+
    stat_density_2d(data = ddf[ddf$swap>2,],aes(x=PC1, y=PC2, fill = as.factor(swap-3), alpha = ..level..),
                    geom = "polygon")
  ggsave(file.path("img/PCA_w_regular.pdf"), width=5, height=10) 
  
  # plot posture
  if(F){
    i = 800
    {
      expand = 1
      df_temp_f <- subset(df_all, id == id_list[1] & sex == "f")[i,]
      df_temp_m <- subset(df_all, id == id_list[1] & sex == "m")[i,]
      xrange = mean(as.numeric(c(df_temp_f[,c(2,4,6)],df_temp_m[,c(2,4,6)])))
      yrange = mean(as.numeric(c(df_temp_f[,c(3,5,7)],df_temp_m[,c(3,5,7)])))
      plot(0, pch="none", xlim=c(xrange-expand, xrange+expand), 
           ylim=c(yrange-expand, yrange+expand), col="#58c8acff",)
      dftempplot = as.numeric(c(df_temp_f[,c(2:7)], df_temp_m[,c(2:7)]))
      draw.arrow = function(x1,y1, x2, y2){ 
        arrows(dftempplot[x1-2], dftempplot[y1-2],
               dftempplot[x2-2], dftempplot[y2-2], length=0)}
      draw.arrow(3,4,5,6)
      draw.arrow(5,6,7,8)
      draw.arrow(9,10,11,12)
      draw.arrow(11,12,13,14)
      
      draw.arrow(3,4,9,10)
      draw.arrow(3,4,11,12)
      draw.arrow(3,4,13,14)
      
      draw.arrow(5,6,9,10)
      draw.arrow(5,6,11,12)
      draw.arrow(5,6,13,14)
      
      draw.arrow(7,8,9,10)
      draw.arrow(7,8,11,12)
      draw.arrow(7,8,13,14)
      
      draw.arrow(5,6,7,8)
      draw.arrow(9,10,11,12)
      draw.arrow(11,12,13,14)
      
      points(dftempplot[3],dftempplot[4], col="#2c7dd3ff", pch=19)
      points(dftempplot[1],dftempplot[2], col="#7009d3ff", pch=19)
      points(dftempplot[5],dftempplot[6], col="#58c8acff", pch=19)
      points(dftempplot[7],dftempplot[8], col="#8fb972ff", pch=19)
      points(dftempplot[9],dftempplot[10], col="#b86d3bff", pch=19)
      points(dftempplot[11],dftempplot[12], col="#d51315ff", pch=19)
    }
  }
}
#------------------------------------------------------------------------------#
