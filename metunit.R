library(meta)
event.e = c(23,25,43)
n.e = c(102,74,100)
event.c = c(23,25,43)
n.c = c(75,65,80)
study = c("A","B","C")

study1 = data.frame(study,event.e,n.e,event.c,n.c)

 metabin(event.e, n.e, event.c, n.c,
                  data = study1, sm = "OR")
 
 
  e = c(84,25,43)
 ne = c(102,74,100)
 c = c(23,25,43)
 nc = c(75,65,80)
  study = c("A","B","C")
 study1 = data.frame(e,ne,c,nc,study)
 library(metafor)
 rma(measure = "OR" ,ai = e, n1i =ne ,ci =c,n2i = nc , data = study1)