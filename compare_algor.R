# libraries
library(factoextra) # the clustring algorithms (pam, kmeans, hac ...etc)
library(clusterCrit) # the package with numerous clustering validation crits

# The call to the function
linkage_methods = 'complete' # for the HAC clustering algo
plot_criteria_comparison(R.vect = R.vect, krange = c(2:10), methods = c("pam", "kmeans", "hclust"), linkage_methods = linkage_methods)



# function to compare diffrent cluster algoritms
plot_criteria_comparison<- function(R.vect, krange, methods = c("pam", "kmeans", "hclust", "fanny"), linkage_methods){
  
  if (missing(linkage_methods)) {
    stop("provide hclust linkage method")
  }
  distance = dist(as.matrix(R.vect), method = 'manhattan')
  crit_summary <- data.frame(array(NA, c(max(krange),5))) # to store the crits
  cluster_crits = c("silhouette","dunn","xie_beni","s_dbw","davies_bouldin") # the critesrias
  names(crit_summary) <- cluster_crits
  
  # looping on the algorithms
  for (method in methods) {
    if (method == "pam") {
      for (k in krange) {
        mres<- eclust(as.matrix(R.vect), "pam", k=k ,  graph = FALSE) # fucntion in factorextra package
        result=intCriteria(traj = as.matrix(R.vect), part = mres$cluster , cluster_crits) # function in Cluster crit package
        crit_summary$silhouette[k] = result$silhouette
        crit_summary$dunn[k] = result$dunn
        crit_summary$xie_beni[k] = result$xie_beni
        crit_summary$s_dbw[k] = result$s_dbw
        crit_summary$davies_bouldin[k] = result$davies_bouldin
        
        pam_summary = crit_summary[-1,] #drop the first row, containing NA
      }
      
    } else if(method == "kmeans") {
      for (k in krange) {
        mres<- eclust(as.matrix(R.vect), "kmeans", k=k ,
                      graph = FALSE)
        result=intCriteria(traj = as.matrix(R.vect), part = mres$cluster ,cluster_crits)
        crit_summary$silhouette[k] = result$silhouette
        crit_summary$dunn[k] = result$dunn
        crit_summary$xie_beni[k] = result$xie_beni
        crit_summary$s_dbw[k] = result$s_dbw
        crit_summary$davies_bouldin[k] = result$davies_bouldin
        kmeans_summary = crit_summary[-1,]
      }
      
    }else if(method == "fanny") {
      for (k in krange) {
        mres<- eclust(as.matrix(R.vect), "fanny", k=k ,
                      graph = FALSE)
        result=intCriteria(traj = as.matrix(R.vect), part = mres$cluster , cluster_crits)
        crit_summary$silhouette[k] = result$silhouette
        crit_summary$dunn[k] = result$dunn
        crit_summary$xie_beni[k] = result$xie_beni
        crit_summary$s_dbw[k] = result$s_dbw
        crit_summary$davies_bouldin[k] = result$davies_bouldin
        fanny_summary = crit_summary[-1,]
      }
      
    }else if(method == "hclust") {
      for (k in krange) {
        mres<- eclust(as.matrix(R.vect), "hclust", k=k , hc_metric = "manhattan",hc_method = linkage_methods, graph = FALSE)
        result=intCriteria(traj = as.matrix(R.vect), part = mres$cluster , cluster_crits)
        results=cluster.stats(d = distance, clustering = mres$cluster)
        crit_summary$silhouette[k] = results$avg.silwidth
        crit_summary$dunn[k] = result$dunn
        crit_summary$xie_beni[k] = result$xie_beni
        crit_summary$s_dbw[k] = result$s_dbw
        crit_summary$davies_bouldin[k] = result$davies_bouldin
        hac_summary = crit_summary[-1,]
      }
      
    }
  }
  # summary of the crits
  best_crit = data.frame(array(NA, dim=c(5,2)))
  colnames(best_crit) = c('method',"K")
  rownames(best_crit) = cluster_crits
  crits = list()
  for (i in 1:length(cluster_crits)) {
    crits[[i]] =data.frame(pam_summary[,i],kmeans_summary[,i],hac_summary[,i])
    if (i %in% c(1,2)) { # (1,2)  correspond to silhoutte and Dunn that are maximized
      best_algo =which.max(sapply(crits[[i]], FUN = max, na.rm=T))
      best_crit[i,1] = methods[best_algo]
      best_crit[i,2] =which.max(crits[[i]][,best_algo])+1
    } else { # for the other three cluster crits that are minimized
      best_algo =which.min(sapply(crits[[i]], FUN = min, na.rm=T))
      best_crit[i,1] = methods[best_algo]
      best_crit[i,2] =which.min(crits[[i]][,best_algo])+1
    }
  }
  
  ## here to generate the plots
  max_d=max(c(pam_summary$dunn, hac_summary$dunn, kmeans_summary$dunn))#, fanny_summary$dunn))
  dunns = ggplot() + 
    geom_point(data = pam_summary, aes(x = as.factor(krange),  y=dunn, color = "pam"))+
    #geom_point(data = fanny_summary, aes(x = as.factor(krange),  y=dunn, color = "fanny"))+
    geom_point(data = kmeans_summary, aes(x = as.factor(krange),  y=dunn, color = "kmeans"))+
    geom_point(data = hac_summary, aes(x = as.factor(krange),  y=dunn, color = "hac"))+ scale_y_continuous(limits =c(0,max_d))+
    labs(x = "Number of clusters") +theme_bw() +theme(legend.title = element_blank() )
  max_x=max(c(pam_summary$xie_beni, hac_summary$xie_beni, kmeans_summary$xie_beni))#, fanny_summary$xie_beni))
  xie =ggplot() + 
    geom_point(data = pam_summary, aes(x = as.factor(krange),  y=xie_beni, color = "pam"))+
    #geom_point(data = fanny_summary, aes(x = as.factor(krange),  y=xie_beni, color = "fanny"))+
    geom_point(data = kmeans_summary, aes(x = as.factor(krange),  y=xie_beni, color = "kmeans"))+
    geom_point(data = hac_summary, aes(x = as.factor(krange),  y=xie_beni, color = "hac"))+ scale_y_continuous(limits =c(0,max_x))+
    labs(x = "Number of clusters")+theme_bw() +theme(legend.title = element_blank() )
  max_c=max(c(pam_summary$s_dbw, hac_summary$s_dbw, kmeans_summary$s_dbw))#, fanny_summary$s_dbw))
  cal =ggplot() + 
    geom_point(data = pam_summary, aes(x = as.factor(krange),  y=s_dbw, color = "pam"))+
    #geom_point(data = fanny_summary, aes(x = as.factor(krange),  y=s_dbw, color = "fanny"))+
    geom_point(data = kmeans_summary, aes(x = as.factor(krange),  y=s_dbw, color = "kmeans"))+
    geom_point(data = hac_summary, aes(x = as.factor(krange),  y=s_dbw, color = "hac"))+ scale_y_continuous(limits =c(0,max_c))+
    labs(x = "Number of clusters")+theme_bw() +theme(legend.title = element_blank() )
  max_db=max(c(pam_summary$davies_bouldin, hac_summary$davies_bouldin, kmeans_summary$davies_bouldin))#, fanny_summary$davies_bouldin))
  db =ggplot() + 
    geom_point(data = pam_summary, aes(x = as.factor(krange),  y=davies_bouldin, color = "pam"))+
    #geom_point(data = fanny_summary, aes(x = as.factor(krange),  y=davies_bouldin, color = "fanny"))+
    geom_point(data = kmeans_summary, aes(x = as.factor(krange),  y=davies_bouldin, color = "kmeans"))+
    geom_point(data = hac_summary, aes(x = as.factor(krange),  y=davies_bouldin, color = "hac"))+ scale_y_continuous(limits =c(0,max_db))+
    labs(x = "Number of clusters")+theme_bw() +theme(legend.title = element_blank() )
  max_s=max(c(pam_summary$silhouette, hac_summary$silhouette, kmeans_summary$silhouette))#, fanny_summary$silhouette))
  sil =ggplot() + 
    geom_point(data = pam_summary, aes(x = as.factor(krange),  y=silhouette, color = "pam"))+
    geom_line(data = pam_summary, aes(x = as.factor(krange),  y=silhouette, color = "pam", group=1))+
    #geom_point(data = fanny_summary, aes(x = as.factor(krange),  y=silhouette, color = "fanny"))+
    #geom_line(data = fanny_summary, aes(x = as.factor(krange),  y=silhouette, color = "fanny", group=1))+
    geom_point(data = kmeans_summary, aes(x = as.factor(krange),  y=silhouette, color = "kmeans"))+
    geom_line(data = kmeans_summary, aes(x = as.factor(krange),  y=silhouette, color = "kmeans", group=1))+
    geom_point(data = hac_summary, aes(x = as.factor(krange),  y=silhouette, color = "hac"))+
    geom_line(data = hac_summary, aes(x = as.factor(krange),  y=silhouette, color = "hac", group=1))+scale_y_continuous(limits =c(0,max_s))+
    labs(x = "Number of clusters", y = "average silhoutte")+theme_bw() +theme(legend.title = element_blank() )
 
  # the elbow plots
  df_elbow  <- data.frame(array(NA, c(max(krange),4)))
  colnames(df_elbow) =c('pam', 'kmeans','hac', 'fanny')
  a=fviz_nbclust(as.matrix(R.vect), pam, method = "wss", k.max = max(krange))
  df_elbow$pam = a$data$y
  a=fviz_nbclust(as.matrix(R.vect), kmeans, method = "wss", k.max = max(krange))
  df_elbow$kmeans = a$data$y
  a=fviz_nbclust(as.matrix(R.vect), hcut, method = "wss", k.max = max(krange))
  df_elbow$hac = a$data$y
  a=fviz_nbclust(as.matrix(R.vect), fanny, method = "wss", k.max = max(krange))
  df_elbow$fanny = a$data$y
  df_elbow$cluster = a$data$clusters
  
  wss= ggplot(data = df_elbow, aes(x=cluster)) + 
    geom_point(aes(y=pam, color = "pam"))+
    geom_line(aes(y=pam, color = "pam", group=1))+ 
    # geom_point(aes(y=fanny, color = "fanny"))+
    # geom_line(aes(y=fanny, color = "fanny", group=1)) + 
    geom_point(aes(y=kmeans, color = "kmeans"))+
    geom_line(aes(y=kmeans, color = "kmeans", group=1)) + 
    geom_point(aes(y=hac, color = "hac"))+
    geom_line(aes(y=hac, color = "hac", group=1))+
    labs(y="Total within sum of squares", x = "Number of clusters") + theme_bw()+ theme(legend.title = element_blank() )
  
  return(list("plots"=ggarrange(wss, sil,dunns, xie, db, cal,  common.legend = T), "summary" = best_crit) )
}

