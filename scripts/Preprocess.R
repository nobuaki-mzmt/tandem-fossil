## Data analysis for C.formosanus tandem movements on sticky trap
## N. Mizumoto

#------------------------------------------------------------------------------#
# This file is for preprocess all data for statistical analysis and plot
#------------------------------------------------------------------------------#

rm(list = ls())
preprocessAll()

#------------------------------------------------------------------------------#
{
  library(data.table)
  library(stringr)
  library(dlcpr)
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
preprocessAll <- function(){
  data.convert()
  data.convert_sleap()
  get.spatial.organization()
  get.movement.speed()
  fossil.dataset()
  posture.organization()
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# This function read csv files for sticky trap and control, and then
# 1) scale from pixel to bodylength
# 2) save as rda
#------------------------------------------------------------------------------#
data.convert <- function(Plot=T, Dataframe = T){
  
  Max.sec <- 1800 # only focus on first 30 minutes or less
  
  ## file location
  raw.place <- file.path("data/", c("Cf-sticky", "Cf-control"))

  ## body size data in pixel
  ## scale all coordinate data in the unit of bodylength
  d.scale = data.frame(fread("data/bodylength-in-pixel.csv", header=T))

  ## calculation
  df.pos <- data.frame()
  df.arrow <- data.frame()
  for(i in 1:2){ # sticky / control
    ## data
    {
      rawdata <- list.files(raw.place[i], full.names = TRUE, pattern = ".csv")
      dataname <- list.files(raw.place[i], full.names = F, pattern = ".csv")
      
      position.data <- rawdata[str_detect(dataname, "position")]
      arrow.data <- rawdata[str_detect(dataname, "arrow")]
      position.data.name <- dataname[str_detect(dataname, "position")]
      arrow.data.name <- dataname[str_detect(dataname, "arrow")]
    }
    
    for(v in 1:length(position.data)){
      # file info
      if(i == 1){ # sticky
        date = str_sub(position.data.name[v], start=1, end=6)
        rep = substr(position.data.name[v], 14, 15)
        event = substr(position.data.name[v], 17, 17) 
        if(event == "p"){ event = 1}
        d.scale.temp = subset(d.scale, Treat=="Sticky" & Date == date &
                                Rep == as.numeric(rep) & Event == as.numeric(event))
        treat = "sticky"
        name = paste(treat,date,rep,event, sep="-")
      } else { # control
        id = str_sub(position.data.name[v], start=13, end=14)
        d.scale.temp = subset(d.scale, Treat=="Control" & ID == as.numeric(id))
        treat = "control"
        name = paste(treat,id, sep="-")
      }
      
      bodylength = (d.scale.temp$Ind1 + d.scale.temp$Ind2)/2
      fps = d.scale.temp$FPS
        
      # data read and scale
      d.pos <- data.frame(fread(position.data[v], header=T))
      d.arrow <- data.frame(fread(arrow.data[v], header=T))
      d.pos <- d.pos[d.pos$position%%fps == 0,]
      d.arrow <- d.arrow[d.arrow$arrow%%fps == 0,]
      
      video.length     <- dim(d.pos)
      d.pos[,2:5]      <- d.pos[,2:5]   / bodylength
      d.arrow[,2:5]    <- d.arrow[,2:5] / bodylength
      d.arrow[,c(2,4)] <- d.arrow[,c(2,4)] - d.pos[1,2]
      d.arrow[,c(3,5)] <- d.arrow[,c(3,5)] - d.pos[1,3]
      d.pos[,c(2,4)]   <- d.pos[,c(2,4)]   - d.pos[1,2]
      d.pos[,c(3,5)]   <- d.pos[,c(3,5)]   - d.pos[1,3]
      
      print(paste(v, "/", length(position.data), "->", name, "video.length", video.length[1]))
      
      # dataframe
      d.pos[,1]   <- d.pos[,1]/fps
      d.arrow[,1] <- d.arrow[,1]/fps
      d.pos       <- d.pos[d.pos[,1] < min(Max.sec, max(d.pos[,1])),]
      d.arrow     <- d.arrow[d.arrow[,1] < min(Max.sec, max(d.arrow[,1])),]
      colnames(d.pos) = colnames(d.arrow) = c("time", "fx", "fy", "mx", "my")
      d.pos <- data.frame(
        d.pos, name, id = v, treat = treat
      )
      d.arrow <- data.frame(
        d.arrow, name, id = v, treat = treat
      )
      
      df.pos <- rbind(df.pos, d.pos)
      df.arrow <- rbind(df.arrow, d.arrow)
    }
  }
  save(df.pos, file=file.path("data/df.pos.rda"))
  save(df.arrow, file=file.path("data/df.arrow.rda"))
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# This function read csv files for sleap results, and then
# 1) scale from pixel to bodylength
# 2) save as rda
#------------------------------------------------------------------------------#
data.convert_sleap <- function(){
  csv_data <- list.files("data/Cf-control-sleap/csv", full.names = TRUE, pattern = ".csv")
  f_names <- list.files("data/Cf-control-sleap/csv", full.names = FALSE, pattern = ".csv")
  d.scale = data.frame(fread("data/bodylength-in-pixel.csv", header=T))
  
  df_all <- NULL
  for(v in 1:length(csv_data)){
    id = str_split(f_names[v], "_")[[1]][6]
    sex = str_split(str_split(f_names[v], "_")[[1]][8], ".csv")[[1]][1]
    d <- data.frame(fread(csv_data[v], header=T))
    d$headtip_x = (d$headtip_x+d$pronotumfront_x)/2
    d$headtip_y = (d$headtip_y+d$pronotumfront_y)/2
    
    d.scale.temp = subset(d.scale, Treat=="Control-SLEAP" & ID == as.numeric(id))
    bodylength = (d.scale.temp$Ind1 + d.scale.temp$Ind2)/2
    
    d = d[seq(1,60*60*15,60),]
    d$V1 = d$V1/60
    colnames(d)[1] = "time"
    
    d[,2:7] = d[,2:7] / (bodylength)
    
    d$sex = sex
    d$id = id
    
    df_all = rbind(df_all, d)
  }
  row.names(df_all) <- NULL
  
  # get pairwise distances
  id_list = unique(df_all$id)
  df = NULL
  
  for( i in id_list){
    df_temp_f <- subset(df_all, id == i & sex == "f")
    df_temp_m <- subset(df_all, id == i & sex == "m")
    
    df_temp = data.frame(i,
                         df_temp_f[,c("time", "headtip_x", "headtip_y", "pronotumfront_x", "pronotumfront_y",
                                      "abdomentip_x", "abdomentip_y")],
                         df_temp_m[,c("headtip_x", "headtip_y", "pronotumfront_x", "pronotumfront_y",
                                      "abdomentip_x", "abdomentip_y")])
    
    colnames(df_temp) = c("name", "time", 
                          "fhead_x", "fhead_y",
                          "fpron_x", "fpron_y",
                          "ftip_x", "ftip_y",
                          "mhead_x", "mhead_y",
                          "mpron_x", "mpron_y", 
                          "mtip_x", "mtip_y")
    df = rbind(df, df_temp)
  }
  
  for(i in seq(3,13,2)){
    for(j in seq(3,13,2)){
      if( i < j){
        dis = sqrt( (df[,i] - df[,j])^2 +
                      (df[,i+1] - df[,j+1])^2)
        cname = paste( str_remove(colnames(df)[i], "_x"), 
                       str_remove(colnames(df)[j], "_x"),
                       sep="_")
        df.temp = data.frame(dis)
        colnames(df.temp) = cname
        df = cbind(df, df.temp)
      }
    }
  }
  
  df_sleap_dis <- na.omit(df)[c(1:2, 15:29)]
  
  # create surrogate data
  swapname = str_replace_all(colnames(df_sleap_dis)[3:17], c("f"="z", "m" = "f", "z" = "m"))
  revfm = str_detect(swapname, "m.*_f.*")
  swapname[revfm] = paste(str_remove(str_extract(swapname[revfm], "_.*"), "_"),
                          str_remove(str_extract(swapname[revfm], ".*_"), "_"),
                          sep="_")
  df_sleap_dis_swap = df_sleap_dis[,c("name","time", swapname)]
  colnames(df_sleap_dis_swap) = colnames(df_sleap_dis)
  
  df_sleap_dis$swap = 3
  df_sleap_dis_swap$swap = 4
  df_sleap_dis = rbind(df_sleap_dis,df_sleap_dis_swap)
  
  f.name = "data/Cf-control-sleap/rda/df_sleap_dis.rda"
  save(df_sleap_dis, file = f.name)
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# This function read df.pos.rda, and df.arrow.rda and then
# 1)  calculate relative position (distance and relative position)
#     and output as dplot_relativepos.rda
# 2)  calculate relative orientation and output as dplot_relativeorient.rda
#------------------------------------------------------------------------------#
get.spatial.organization <- function(){
  
  get.relative.pos <- function(ax, ay, rx, ry, treat, direction, time=NA, id=NA){
    angle <- atan2(ay, ax)
    rp <- atan2(ry, rx)
    
    rpa <- rp-angle
    
    rpa[rpa > pi] <- rpa[rpa > pi] - 2*pi
    rpa[rpa < -pi] <- rpa[rpa < -pi] + 2*pi
    
    dis <- sqrt(rx^2 + ry^2)
    rpx <- cos(rpa)*dis
    rpy <- sin(rpa)*dis
    
    return(data.frame(dis, angle, rpx, rpy, treat, direction, time, id))
  }
  
  ## Original data (sticky trap)
  {
    load("data/df.pos.rda")
    load("data/df.arrow.rda")
    df.pos <- df.pos[df.pos$treat=="sticky",]
    df.arrow <- df.arrow[df.arrow$treat=="sticky",]
    
    # From female to male
    ax <- df.arrow$fx - df.pos$fx
    ay <- df.arrow$fy - df.pos$fy
    rx <- df.pos$mx - df.pos$fx
    ry <- df.pos$my - df.pos$fy
    df.temp1 <- get.relative.pos(ax, ay, rx, ry,
                                 treat="sticky", direction = "FtoM",
                                 time=df.pos$time, id=df.pos$id)
    
    # From male to female
    ax <- df.arrow$mx - df.pos$mx
    ay <- df.arrow$my - df.pos$my
    rx <- df.pos$fx - df.pos$mx
    ry <- df.pos$fy - df.pos$my
    df.temp2 <- get.relative.pos(ax, ay, rx, ry, 
                                 treat="sticky", direction = "MtoF",
                                 time=df.pos$time, id=df.pos$id)
    
    # Orientation
    f.orient = atan2(df.arrow$fy - df.pos$fy,  df.arrow$fx - df.pos$fx)
    m.orient = atan2(df.arrow$my - df.pos$my,  df.arrow$mx - df.pos$mx)
    orient.diff = m.orient - f.orient
    orient.diff[orient.diff > pi] = orient.diff[orient.diff > pi] - 2*pi
    orient.diff[orient.diff < -pi] = orient.diff[orient.diff < -pi] + 2*pi
    
    # others
    dis = sqrt(rx^2+ry^2)
    
    # output
    df.temp3 <- data.frame(
      relative.direction = abs(atan2(df.temp1$rpy, df.temp1$rpx)), 
      relative.orientation = abs(orient.diff),
      dis, treat="sticky", direction = "FtoM", time=df.pos$time, id=df.pos$id
      )
    
    df.temp4 <- data.frame(
      relative.direction = abs(atan2(df.temp2$rpy, df.temp2$rpx)), 
      relative.orientation = abs(orient.diff),
      dis, treat="sticky", direction = "MtoF", time=df.pos$time, id=df.pos$id
      )
    
    df.plot.sticky <- rbind(df.temp1, df.temp2)
    df.orient.sticky <- rbind(df.temp3, df.temp4)
  }
  
  ## Surrogate data (sticky trap)
  {
    load("data/df.pos.rda")
    load("data/df.arrow.rda")
    df.pos <- df.pos[df.pos$treat=="sticky",]
    df.arrow <- df.arrow[df.arrow$treat=="sticky",]
    
    iteration = 1000
    set.seed(1) 
    fem.random <- sample(1:24, iteration, replace=T)
    set.seed(2) 
    mal.random <- sample(1:24, iteration, replace=T)
    set.seed(3) 
    direction <- runif(iteration, 0, 1)*pi*2
    
    videos <- unique(df.pos$id)
    df.plot = df.orient <- data.frame()
    for(i in 1:iteration){
      print(paste(i, "/", iteration))
      f.df.pos <- df.pos[df.pos$id == fem.random[i],]
      m.df.pos <- df.pos[df.pos$id == mal.random[i],]
      f.df.arrow <- df.arrow[df.arrow$id == fem.random[i],]
      m.df.arrow <- df.arrow[df.arrow$id == mal.random[i],]
      
      cut.dim <- min(dim(f.df.pos)[1],   dim(m.df.pos)[1],
                     dim(f.df.arrow)[1], dim(m.df.arrow)[1])
      f.df.pos <- f.df.pos[1:cut.dim,]
      m.df.pos <- m.df.pos[1:cut.dim,]
      f.df.arrow <- f.df.arrow[1:cut.dim,]
      m.df.arrow <- m.df.arrow[1:cut.dim,]
      
      fx <- f.df.pos$fx
      fy <- f.df.pos$fy
      
      dis <- sqrt(f.df.pos$mx[1]^2+f.df.pos$my[1]^2)
      
      mx <- m.df.pos$mx - m.df.pos$mx[1] + cos(direction[i])*dis
      my <- m.df.pos$my - m.df.pos$my[1] + sin(direction[i])*dis
      
      fax <- f.df.arrow$fx
      fay <- f.df.arrow$fy
      max <- m.df.arrow$mx - m.df.pos$mx[1] + cos(direction[i])*dis
      may <- m.df.arrow$my - m.df.pos$my[1] + sin(direction[i])*dis
      
      ax <- fax - fx
      ay <- fay - fy
      rx <- mx - fx
      ry <- my - fy
      df.temp1 <- get.relative.pos(ax, ay, rx, ry, 
                                   treat="surrogate", direction = "FtoM", 
                                   time=NA, id=NA)
      
      ax <- max - mx
      ay <- may - my
      rx <- fx - mx
      ry <- fy - my
      df.temp2 <- get.relative.pos(ax, ay, rx, ry, 
                                   treat="surrogate", direction = "MtoF", 
                                   time=NA, id=NA)
      
      df.plot <- rbind(df.plot, df.temp1, df.temp2)
      
      
      # Orientation
      f.orient = atan2(fay - fy,  fax - fx)
      m.orient = atan2(may - my,  max - mx)
      orient.diff = m.orient - f.orient
      orient.diff[orient.diff > pi] = orient.diff[orient.diff > pi] - 2*pi
      orient.diff[orient.diff < -pi] = orient.diff[orient.diff < -pi] + 2*pi
      
      # others
      dis = sqrt(rx^2+ry^2)
      
      # output
      df.temp3 <- data.frame(
        relative.direction = abs(atan2(df.temp1$rpy, df.temp1$rpx)), 
        relative.orientation = abs(orient.diff),
        dis, treat="surrogate", direction = "FtoM", time=NA, id=NA
        )
      
      df.temp4 <- data.frame(
        relative.direction = abs(atan2(df.temp2$rpy, df.temp2$rpx)), 
        relative.orientation = abs(orient.diff),
        dis = dis, treat="surrogate", direction = "MtoF", time=NA, id=NA)
      df.orient <- rbind(df.orient, df.temp3, df.temp4)
    }
    
    df.plot.surrogate <- df.plot
    df.orient.surrogate <- df.orient
  }
  
  ## Original data (control tandem)
  {
    load("data/df.pos.rda")
    load("data/df.arrow.rda")
    df.pos <- df.pos[df.pos$treat=="control",]
    df.arrow <- df.arrow[df.arrow$treat=="control",]
    
    # female to male
    ax <- df.arrow$fx - df.pos$fx
    ay <- df.arrow$fy - df.pos$fy
    rx <- df.pos$mx - df.pos$fx
    ry <- df.pos$my - df.pos$fy
    
    df.temp1 <- get.relative.pos(ax, ay, rx, ry,
                                 treat="control", direction="FtoM",
                                 time=df.pos$time, id=df.pos$id)
    
    # male to female
    ax <- df.arrow$mx - df.pos$mx
    ay <- df.arrow$my - df.pos$my
    rx <- df.pos$fx - df.pos$mx
    ry <- df.pos$fy - df.pos$my
    
    df.temp2 <- get.relative.pos(ax, ay, rx, ry, 
                                 treat="control", direction="MtoF",
                                 time=df.pos$time, id=df.pos$id)
    
    df.plot.control <- rbind(df.temp1, df.temp2)
    
    # Orientation
    f.orient = atan2(df.arrow$fy - df.pos$fy,  df.arrow$fx - df.pos$fx)
    m.orient = atan2(df.arrow$my - df.pos$my,  df.arrow$mx - df.pos$mx)
    orient.diff = m.orient - f.orient
    orient.diff[orient.diff > pi] = orient.diff[orient.diff > pi] - 2*pi
    orient.diff[orient.diff < -pi] = orient.diff[orient.diff < -pi] + 2*pi
    dis = sqrt(rx^2+ry^2)
    df.temp3 <- data.frame(
      relative.direction = abs(atan2(df.temp1$rpy, df.temp1$rpx)),
      relative.orientation = abs(orient.diff),
      dis = dis, treat="control", direction = "FtoM", time=df.pos$time, id=df.pos$id)
    df.temp4 <- data.frame(
      relative.direction = abs(atan2(df.temp2$rpy, df.temp2$rpx)), 
      relative.orientation = abs(orient.diff),
      dis = dis, treat="control", direction = "MtoF", time=df.pos$time, id=df.pos$id)
    
    df.orient.control <- rbind(df.temp3, df.temp4)
    
  }
  
  
  df.plot.relative.pos <- 
    rbind(df.plot.sticky, df.plot.control, df.plot.surrogate)
  df.plot.relative.orient <- 
    rbind(df.orient.sticky, df.orient.control, df.orient.surrogate)
  save(df.plot.relative.pos, file="data/dplot_relativepos.rda")
  save(df.plot.relative.orient, file="data/dplot_relativeorient.rda")
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# This function read df.pos.rda, and then
# 1) calculate speed
#------------------------------------------------------------------------------#
get.movement.speed <- function(){
  load("data/df.pos.rda")
  fsl <- c(NA, sqrt( diff(df.pos$fx)^2 + diff(df.pos$fy)^2 ))
  msl <- c(NA, sqrt( diff(df.pos$mx)^2 + diff(df.pos$my)^2 ))
  
  rx  <- df.pos$fx - df.pos$mx
  ry  <- df.pos$fy - df.pos$my
  dis <- sqrt(rx^2+ry^2)
  
  fsl[df.pos$time==0] = NA
  msl[df.pos$time==0] = NA
  
  df.temp <- data.frame(
    fspeed = fsl,
    mspeed = msl,
    dis    = dis, 
    treat  = df.pos$treat,
    time   = df.pos$time, 
    id     = df.pos$id
  )
  dff <- tapply(df.temp$fspeed, df.temp[,c("id", "treat")], mean, na.rm=T)
  dfm <- tapply(df.temp$mspeed, df.temp[,c("id", "treat")], mean, na.rm=T) 
  df.speed = 
    rbind(
      data.frame(
        sex   = "female",
        treat = "control",
        speed = dff[,1],
        id    = row.names(dff)
      ),
      data.frame(
        sex   = "female",
        treat = "sticky",
        speed = dff[,2],
        id    = row.names(dff)
      ),
      data.frame(
        sex   = "male",
        treat = "control",
        speed = dfm[,1],
        id    = row.names(dff)
      ),
      data.frame(
        sex   = "male",
        treat = "sticky",
        speed = dfm[,2],
        id    = row.names(dff)
      )
    )
  df.speed <- na.omit(df.speed)
  save(df.speed, file="data/dfspeed.rda")
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# This function process data from the fossil
#------------------------------------------------------------------------------#
fossil.dataset <- function(){
  BodyLength <- c(469.3591375, 570.0614226)  ## male, female

  mbd <- mean(BodyLength) 
  Position.x <- c(1069.665, 799.07) 
  Position.y <- c(489.159, 435.622) 
  Angle <- c(-78.148, 55.486)/180*pi
  
  df.center.fossil <- data.frame(
    sex = c("male", "female"),
    x = c(1069.665, 799.07)/mbd,
    y = c(489.159, 435.622)/mbd,
    angle = c(-78.148, 55.486)/180*pi
  )
  
  d.fossil <- data.frame(parts = c("fhead", "fpron", "ftip", "mhead", "mpron", "mtip"),
                         x = c(950, 902, 723, 1119, 1109, 1016),
                         y = c(274, 297, 683, 670, 622, 258) )
  df.all.fossil = NULL
  for(i in 1:6){
    for(j in 1:6){
      if(j>i){
        dis = sqrt(    (d.fossil[i,2]-d.fossil[j,2])^2 +
                         (d.fossil[i,3]-d.fossil[j,3])^2 )
        names(dis)  = paste( d.fossil[i,1], d.fossil[j,1],
                             sep="_")
        df.all.fossil <- c(df.all.fossil, dis)
      }
    }
  }
  
  df.all.fossil = df.all.fossil / mean(BodyLength)
  save(df.center.fossil, df.all.fossil, file="data/df_fossil.rda")
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# This function read all csv files for DLC, and then
# 1. measure the distance between all body parts.
# 2. create surrogate datasets and integrate with fossil data
# 3. perform PCA
#------------------------------------------------------------------------------#
posture.organization <- function(){
  idir = "data/Cf-sticky-DLC/"
  rawdataDLC = list.files(idir, full.names = T)
  data.name  = list.files(idir, full.names = F)
  
  ## body size data in pixel
  ## scale all coordinate data in the unit of bodylength
  d.scale = data.frame(fread("data/bodylength-in-pixel.csv", header=T))
  
  # data from observations
  {
    df.all = NULL
    
    for(j in 1:length(rawdataDLC)){
      
      d = read.dlc(rawdataDLC[j], fps = 5, spline=T)
      
      date      = str_sub(data.name[j], start=1, end=6)
      species   = "Cf"
      rep       = substr (data.name[j], 14, 15)
      event     = substr (data.name[j], 17, 17) 
      if(event == "c"){
        event = 1
      }
      
      d.scale.temp = subset(d.scale, Treat=="Sticky-DLC" & Date == date &
                              Rep == as.numeric(rep) & Event == as.numeric(event))
      
      treat = "sticky"
      name = paste(species,treat,date,rep,event, sep="-")
      print(name)
      
      df = data.frame(
        name, 
        d[,c("time", "f1_x", "f1_y", "f2_x", "f2_y", "f3_x", "f3_y", 
             "m1_x", "m1_y",  "m2_x", "m2_y", "m3_x", "m3_y")])
      
      colnames(df) = c("name", "time", 
                       "fhead_x", "fhead_y",
                       "fpron_x", "fpron_y",
                       "ftip_x", "ftip_y",
                       "mhead_x", "mhead_y",
                       "mpron_x", "mpron_y", 
                       "mtip_x", "mtip_y")
      df[,3:14] <- df[,3:14] / ((d.scale.temp$Ind1 + d.scale.temp$Ind2)/2)
      df = df[df$time%%1 == 0,]
      
      df.all <- rbind(df.all, df)
    }
    
    for(i in seq(3,13,2)){
      for(j in seq(3,13,2)){
        if( i < j){
          dis = sqrt( (df.all[,i] - df.all[,j])^2 +
                        (df.all[,i+1] - df.all[,j+1])^2)
          cname = paste( str_remove(colnames(df)[i], "_x"), 
                         str_remove(colnames(df)[j], "_x"),
                         sep="_")
          df.temp = data.frame(dis)
          colnames(df.temp) = cname
          df.all = cbind(df.all, df.temp)
        }
      }
    }
    
    #df.all = df.all[df.all$time%%1 == 0,]
    
    df.pca <- na.omit(df.all)[c(1:2, 15:29)]
    
    # focus on the time when the distance between female head and male head < 1.5 body length
    df.pca = df.pca[ df.pca[,5] < 1.5, ]
    
    # create surrogate data
    swapname = str_replace_all(colnames(df.pca)[3:17], c("f"="z", "m" = "f", "z" = "m"))
    revfm = str_detect(swapname, "m.*_f.*")
    swapname[revfm] = paste(str_remove(str_extract(swapname[revfm], "_.*"), "_"),
                            str_remove(str_extract(swapname[revfm], ".*_"), "_"),
                            sep="_")
    df.pca.swap = df.pca[,c("name","time", swapname)]
    colnames(df.pca.swap) = colnames(df.pca)
  
    df.pca.swap$swap = 1
    df.pca$swap = 0
    df.pca.combined = rbind(df.pca,df.pca.swap)
  }
  
  # data from fossil
  {
    load("data/df_fossil.rda")
    df.fossil = data.frame(name="fossil", time = 0, 
                           data.frame(matrix(df.all.fossil, nrow=1)),
                           swap=2)
    colnames(df.fossil) = c("name", "time", names(df.all.fossil),"swap")
  }
  
  # PCA
  {
    df.pca.combined = rbind(
      df.pca.combined,
      df.fossil
    )
    pl = prcomp(df.pca.combined[,3:17], scale= F)
    summary(pl)
    
    df.pca.res <- cbind(df.pca.combined, pl$x[,1:4])
  }
  
  save(pl, df.all, df.pca.res, file="data/df_pca.rda")
}
#------------------------------------------------------------------------------#
