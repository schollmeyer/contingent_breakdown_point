##depth functions

### Tukey depth


Tukey_depth <- function(X){  ### berechnet Levelfunktion fuer begriffliches Quantilkonzept
  n=dim(X)[1]
  
  colmeans <- colMeans(X)
  
  depths <- rep(0,n)
  for(k in (1:n)){
    indexs <- which(X[k,]==0)
    if(length(indexs)>=1){
      depths[k] <- 1-max(colmeans[indexs])
    }
    else{depths[k] <- 1}
  }
  return(list(depths=depths))}


## peeling depth

nu_min <- function(extent,context,timelimit=Inf){
  m <- sum(extent)
  if(m==1){return(which(extent==1))}
  if(m==1 | all(colSums(context[which(extent==1),]) %in% c(0,m))){return(which(extent==1)[1])}
  intent <- oofos:::compute_psi(extent,context=context)
  ans <- oofos:::min_k_obj_generated(extent,intent,context)
  #  TODO: delete the following line
  ans$vtype <- rep("B",ncol(ans$A))
  bns <- gurobi::gurobi(ans,params=list(outputflag=0,timelimit=timelimit))
return(which(bns$x[-(1:dim(context)[2])]==1))}

peeled_nu_min=function(E,context,K=1,VC=Inf){
  i <- NULL
  j <-  NULL
  for(k in (1:K)){
  if(is.null(j) | length(j)==dim(context)[2]){i <- nu_min(E,context)}
  else{i <- nu_min(E,context[,-j])}
  if(length(i)<=VC){return(i)}
  j <- c(j, which(colSums(matrix(context[i,],nrow=length(i)))==length(i)-1))
}
return(i)}

## peeling depth

peeling_depth <- function(context,peeling_operator=nu_min,...){
  m <- dim(context)[1]
  ans <- rep(0,m)
  idx <- (1:m)
  k <- 1
  context2 <- context
  while(TRUE){
 
  i <- peeling_operator(rep(1,dim(context2)[1]),context2,...)
  if(length(i)==0){return(ans)}
  ans[idx[i]] <- k
  idx <- idx[-i]
  if(length(idx)==0){
    number_of_layers <- length(unique(ans))
    return(list(depths=ans,number_of_layers=number_of_layers,cbp_bound = number_of_layers/m))
  }
   k <- k+1  
   context2 <- matrix(context2[-i,],nrow=dim(context2)[1]-length(i))
  
  }
  }
  ######


### enclosing depth

compute_enclosing <- function(context,point_indexs){
   model <- oofos:::compute_extent_vc_dimension(context)
   model$vtype=model$vtypes
   m <- nrow(context)
   n <- ncol(context) 
   #if(class(point_indexs)=="numeric"){point_indexs=matrix(point_indexs,nrow=1)}
   if(length(point_indexs)==1) {zero_indexs=which(context[point_indexs,]==0)}
   else{zero_indexs <- which(colSums(context[point_indexs,]) < length(point_indexs))}
   if(length(zero_indexs)==0){print("arghh");return(NULL)} 
   model$ub <- c(rep(1,m),rep(0,n))
   model$ub[m+ zero_indexs] <- 1
   model$ub[point_indexs] <- 0
   for(k in zero_indexs){
     i <- which(context[,k]==0)
     if(length(i)==0){prinT("GGGG")}
     temp <- rep(0,m+n)
     temp[i] <- 1
     model$A <- rbind(model$A, matrix(temp,nrow=1))
     model$rhs <- c(model$rhs,1)
     model$sense <- c(model$sense,">=")
   }
   model$obj <- NULL
   result <- gurobi::gurobi(model,list(PoolSearchMode=2,PoolSolutions=100000000,outputflag=1))
   print(result$status)
   if(result$status=="OPTIMAL")  {
    m <- nrow(context)
    n <- ncol(context)
    N <- length(result$pool)
    mat2 <- mat <- array(0,c(N,m)) 
    for(k in (1:N)){
      mat[k,] <- round((result$pool[[k]])$poolnx[(1:m)],2)
    }

    for(k in (1:N)){
      mat2[k,] <- oofos:::compute_phi(oofos:::compute_psi(mat[k,],context),context)
    }
     
     rowsums <- rowSums(mat2)
     i <- which.min(rowsums)
     return(list(mat=mat,model=model,result=result,extent=mat[i,]))
   }
   
   else(return(list(result=result)))}

enclosing_depth <- function(context,startindex){
  depths <- rep(0,nrow(context))
  e <- rep(0,nrow(context))
  e[startindex] <- 1
  depths[startindex] <- 0
  t <- 2
  for(k in seq_len(nrow(context))){
    temp <- compute_enclosing(context,which(e==1))
    if(temp$result$status!="OPTIMAL"){return(depths)}
    depths[which(temp$extent==1)] <- t 
    e[which(depths!=0)] <- 1
    t <- t+1 




  }}
  

### Tverberg depth
## Tverberg depth in R^2



Tverberg_depth_basic_model <- function(X){
  n <- nrow(X)
  N <- choose(n,3)
  model <- list()
  obj <- rep(0,N)
  t <- 1
  indexs <- array(0,c(N,3))
  sos <- list()
  print(Sys.time())
  for(k in seq_len(n-2)){
   for(l in seq(k+1,n-1)){
     for(m in seq(l+1,n)){
	   indexs[t,] <- c(k,l,m)
	     t <- t+1
		 }}}
		 
		 
		 
	
	t <- 1
	for(k in seq_len(n)){
	index <- which(indexs[,1]==k | indexs[,2]==k | indexs[,3]==k)
	sos[[k]] <- list(type=1,index=index,weight=rep(1/length(index),length(index)))	
	}
	print( Sys.time())

	

return(list(indexs=indexs,A=Matrix::Matrix(0, nrow = 0, ncol = N, sparse = TRUE),obj=obj,sos=sos,lb=rep(0,N),ub=rep(1,N),vtype=rep("B",N),modelsense="max"))}

update_Tverberg_depth_model <- function(model,X,point_index){
  n <- nrow(X)
  N <- choose(n,3)
  new_model <- model
  new_model$obj <- rep(0,N)
  t <- 1
  for(k in seq_len(n-2)){
   for(l in seq(k+1,n-1)){
     for(m in seq(l+1,n)){
	     if( (point_index %in% c(k,l,m)) | is_in_convex_hull(X[point_index,],X[k,],X[l,],X[m,])     ){new_model$obj[t] <- 1}
		 #if(point_index %in% c(k,l,m)){new_model$ub[t] <- 0}
		 t <- t+1
		 }}}
return(new_model)}		 




## Tverberg depth in R^3
Tverberg_depth_basic_model3d <- function(X){
  n <- nrow(X)
  N <- choose(n,4)
  model <- list()
  obj <- rep(0,N)
  t <- 1
  indexs <- array(0,c(N,4))
  sos <- list()
  print(Sys.time())
  for(k in seq_len(n-3)){
   for(l in seq(k+1,n-2)){
     for(m in seq(l+1,n-1)){
       for(o in seq(m+1,n)){
	     indexs[t,] <- c(k,l,m,o)
	     t <- t+1
       }
		 }}}
		 
		 
		 
	
	t <- 1
	for(k in seq_len(n)){
	index <- which(indexs[,1]==k | indexs[,2]==k | indexs[,3]==k | indexs[,4]==k)
	sos[[k]] <- list(type=1,index=index,weight=rep(1/length(index),length(index)))	
	}
	print( Sys.time())

	

return(list(indexs=indexs,A=Matrix::Matrix(0, nrow = 0, ncol = N, sparse = TRUE),obj=obj,sos=sos,lb=rep(0,N),ub=rep(1,N),vtype=rep("B",N),modelsense="max"))}

update_Tverberg_depth_model3d <- function(model,X,point_index,exclude_point_index=FALSE){
  n <- nrow(X)
  N <- choose(n,4)
  new_model <- model
  new_model$obj <- rep(0,N)
  t <- 1
  for(k in seq_len(n-3)){
   for(l in seq(k+1,n-2)){
     for(m in seq(l+1,n-1)){
      for(o in seq(m+1,n)){
        #chull_indexs <- unique(as.vector(geometry::convhulln(X[c(point_index,k,l,m,o),])))
	      if( (point_index %in% c(k,l,m,o)) | ( is_in_tetrahedron(X[point_index,], X[k,],X[l,],X[m,],X[o,]))){new_model$obj[t] <- 1}
		    if(exclude_point_index & point_index %in% c(k,l,m,o)){new_model$ub[t] <- 0}
		 t <- t+1
      }
		 }}}
return(new_model)}		 


## Tverberg depth in R^4
Tverberg_depth_basic_model4d <- function(X){
  n <- nrow(X)
  N <- choose(n,5)
  model <- list()
  obj <- rep(0,N)
  t <- 1
  indexs <- array(0,c(N,5))
  sos <- list()
  print(Sys.time())
  for(k in seq_len(n-4)){
   for(l in seq(k+1,n-3)){
     for(m in seq(l+1,n-2)){
       for(o in seq(m+1,n-1)){
          for(p in seq(o+1,n)){
	          indexs[t,] <- c(k,l,m,o,p)
	          t <- t+1
          }
       }
		 }}}
		 
		 
		 
	
	t <- 1
	for(k in seq_len(n)){
	index <- which(indexs[,1]==k | indexs[,2]==k | indexs[,3]==k | indexs[,4]==k | indexs[,5]==k)
	sos[[k]] <- list(type=1,index=index,weight=rep(1/length(index),length(index)))	
	}
	print( Sys.time())

	

return(list(indexs=indexs,A=Matrix::Matrix(0, nrow = 0, ncol = N, sparse = TRUE),obj=obj,sos=sos,lb=rep(0,N),ub=rep(1,N),vtype=rep("B",N),modelsense="max"))}

update_Tverberg_depth_model4d <- function(model,X,point_index,exclude_point_index=FALSE){
  n <- nrow(X)
  N <- choose(n,5)
  new_model <- model
  new_model$obj <- rep(0,N)
  t <- 1
  for(k in seq_len(n-4)){
    print(k)
   for(l in seq(k+1,n-3)){
     for(m in seq(l+1,n-2)){
      for(o in seq(m+1,n-1)){
         for(p in seq(o+1,n)){
#print(p)
        #chull_indexs <- unique(as.vector(geometry::convhulln(X[c(point_index,k,l,m,o),])))
	      if( (point_index %in% c(k,l,m,o)) | ( is_in_pentachoron(X[point_index,], X[k,],X[l,],X[m,],X[o,],X[p,]))){new_model$obj[t] <- 1}
        #if( (point_index %in% c(k,l,m,o)) | ( is_point_in_hull_5pt_exact(X[point_index,], rbind(X[k,],X[l,],X[m,],X[o,],X[p,])))){new_model$obj[t] <- 1}
		    if(exlcude_point_index & point_index %in% c(k,l,m,o)){new_model$ub[t] <- 0}
		 t <- t+1
         }}
		 }}}
return(new_model)}		 




## experimental code:

is_in_tetrahedron <- function(x, V1, V2, V3, V4) {
  # Construct the coefficient matrix A (4x4)
  A <- rbind(cbind(V1, V2, V3, V4), c(1, 1, 1, 1))
  
  # Construct the target vector b (4x1)
  b <- c(x, 1)
  
  # Solve for the barycentric coordinates (alpha)
  # Wrap in tryCatch in case the 4 points are coplanar (matrix is singular)
  alpha <- tryCatch({
    solve(A, b)
  }, error = function(e) {
    return(NULL)
  })
  
  # If a valid solution exists, check that all coordinates are >= 0
  # A small tolerance (e.g., -1e-9) handles floating-point precision errors
  #print(alpha)
  #alpha <<- alpha
  if (!is.null(alpha) && all(alpha >= 0)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}




#Pentachoron

is_in_pentachoron <- function(x, V1, V2, V3, V4,V5) {
  # Construct the coefficient matrix A (4x4)
  A <- rbind(cbind(V1, V2, V3, V4, V5), c(1, 1, 1, 1,1))
  
  # Construct the target vector b (4x1)
  b <- c(x, 1)
  
  # Solve for the barycentric coordinates (alpha)
  # Wrap in tryCatch in case the 4 points are coplanar (matrix is singular)
  alpha <- tryCatch({
    solve(A, b)
  }, error = function(e) {
    return(NULL)
  })
  
  # If a valid solution exists, check that all coordinates are >= 0
  # A small tolerance (e.g., -1e-9) handles floating-point precision errors
  if (!is.null(alpha) && all(alpha >= -1e-9)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

##


#neuer versuch:


library(rcdd)
library(MASS)
library(rcdd)

# 1. Hilfsfunktion zur Bereinigung von Strings für GMP (keine Leerzeichen, expliziter Nenner)
safe_double_to_rational <- function(x) {
  # Konvertiert Doubles oder Strings mit Dezimalpunkt exakt in GMP-Rationale Objekte
  bigq_obj <- gmp::as.bigq(x)
  fracs <- as.character(bigq_obj)
  
  # Falls eine Ganzzahl ohne Nenner ausgeworfen wird (z.B. "2"), "/1" anhängen
  has_no_denom <- !grepl("/", fracs)
  fracs[has_no_denom] <- paste0(fracs[has_no_denom], "/1")
  
  return(fracs)
}

library(rcdd)

is_point_in_hull_5pt_exact <- function(query_point, cloud_5_points) {
  # cloud_5_points: Charakter-Matrix (5 Zeilen, 4 Spalten) mit Brüchen (z.B. "1/3")
  # query_point: Charakter-Vektor (Länge 4) mit Brüchen
  
  # Formatierung für GMP sichern (Leerzeichen entfernen, Brüche erzwingen)
  query_clean <- gsub(" ", "", as.character(query_point))
  query_clean[!grepl("/", query_clean)] <- paste0(query_clean[!grepl("/", query_clean)], "/1")
  
  cloud_clean <- matrix(gsub(" ", "", as.character(cloud_5_points)), nrow = 5)
  cloud_clean[!grepl("/", cloud_clean)] <- paste0(cloud_clean[!grepl("/", cloud_clean)], "/1")
  
  # 1. H-Repräsentation der Hyperebenen nach offiziellem rcdd-Standard aufbauen:
  # Zeilenformat: [0 (für <=) | 0 (Code) | 1 (Konstante) | -Koordinaten]
  # Für die 5 Punkte der Wolke:
  hrep <- cbind("0", "0", "1", matrix(paste0("-", cloud_clean), nrow = 5))
  # Negative Vorzeichen bereinigen (Doppel-Minus "--" zu "")
  hrep <- gsub("--", "", hrep)
  
  # 2. Den Testpunkt im Gleichungsformat am Ende anhängen [0 | 1 | 1 | -q]
  # Dadurch verankern wir die Zielfunktion mathematisch fest im System
  query_row <- c("0", "1", "1", paste0("-", query_clean))
  query_row <- gsub("--", "", query_row)
  hrep <- rbind(hrep, query_row)
  
  # 3. Zielfunktion definieren: c(-1, q_1, q_2, q_3, q_4)
  # Die Länge entspricht exakt dem Raum (1 Hilfsvariable + 4 Dimensionen = 5 Elemente)
  obj <- c("-1", query_clean)
  
  # 4. Exakten rationalen LP-Solver aufrufen
  # Minimiert nicht, sondern sucht den Maximalwert der Trennung
  out <- lpcdd(hrep, obj, minimize = FALSE)
  
  # 5. Auswertung laut rcdd-Vignette:
  # Ist der Optimalwert > 0, gibt es eine trennende Ebene -> Punkt ist AUẞERHALB.
  # Ist der Optimalwert == "0", liegt der Punkt INNERHALB oder auf dem Rand.
  if (out$optimal.value == "0") {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# --- BLITZSCHNELLER TEST (4D-SIMPLEX) ---

cloud_5 <- matrix(c(
  "0/1", "0/1", "0/1", "0/1",
  "1/1", "0/1", "0/1", "0/1",
  "0/1", "1/1", "0/1", "0/1",
  "0/1", "0/1", "1/1", "0/1",
  "0/1", "0/1", "0/1", "1/1"
), ncol = 4, byrow = TRUE)

# Test 1: Punkt liegt im Inneren (Ergibt TRUE)
print(is_point_in_hull_5pt_exact(c("1/5", "1/5", "1/5", "1/5"), cloud_5))

# Test 2: Punkt liegt knapp außerhalb (Ergibt FALSE)
print(is_point_in_hull_5pt_exact(c("1/2", "1/2", "1/2", "1/2"), cloud_5))





library(rcdd)



is_point_in_hull_exact <- function(query_point, cloud_points) {
  # 1. V-Repräsentation der Punktwolke bauen
  vrep_cloud <- cbind("0", "1", cloud_points)
  
  # 2. Den Testpunkt als eigene Zeile vorbereiten
  vrep_query <- matrix(c("0", "1", query_point), nrow = 1)
  
  # 3. Zusammenfügen. Der Testpunkt ist genau die LETZTE Zeile (Index: target_index)
  combined_vrep <- rbind(vrep_cloud, vrep_query)
  target_index <- nrow(combined_vrep)
  
  # 4. Redundante Punkte berechnen
  result <- redundant(combined_vrep, representation = "V")
  
  # 5. Prüfen, ob der Index unseres Testpunkts in der Liste der Redundanzen auftaucht
  # result$redundant enthält die Zeilennummern der redundanten Punkte
  if (target_index %in% result$redundant) {
    return(TRUE)   # Punkt ist redundant -> liegt INNERHALB der Hülle
  } else {
    return(FALSE)  # Punkt ist nicht redundant -> liegt AUẞERHALB der Hülle
  }
}

library(MASS)



library(rcdd)
library(MASS)

# 1. Sichere Konvertierung von numerischen Werten in GMP-konforme Strings
safe_double_to_rational <- function(x) {
  # In Brüche umwandeln und als Text formatieren
  fracs <- as.character(MASS::fractions(x))
  
  # ALLE Leerzeichen entfernen (wichtig für GMP)
  fracs <- gsub(" ", "", fracs)
  
  # Wenn kein "/" enthalten ist (Ganzzahl), "/1" anhängen
  has_no_denominator <- !grepl("/", fracs)
  fracs[has_no_denominator] <- paste0(fracs[has_no_denominator], "/1")
  
  return(fracs)
}


# Function to safely convert any numeric double vector to rational strings
double_to_rational <- function(x) {
  # 1. Compute fraction representation
  fracs <- MASS::fractions(x)
  
  # 2. Convert to character strings
  str_fracs <- as.character(fracs)
  
  # 3. Handle whole numbers: "rcdd" requires the denominator explicitly (e.g. "5" -> "5/1")
  has_no_denominator <- !grepl("/", str_fracs)
  str_fracs[has_no_denominator] <- paste0(str_fracs[has_no_denominator], "/1")
  
  return(str_fracs)
}
#experimental code
is_in_convex_hull <- function(x, V1, V2, V3) {
  # Create the matrix of vertices
  M <- cbind(V1, V2, V3)
  
  # Add the constraint rows for the system
  A <- rbind(M, c(1, 1, 1))
  b <- c(x, 1)
  
  # Solve the linear system for alpha
  # (Will throw an error if the 3 vertices are collinear)
  alpha <- tryCatch({
    solve(A, b)
  }, error = function(e) {
    return(NULL)
  })
  
  # Check if a valid solution was found and all barycentric coordinates are >= 0
  if (!is.null(alpha) && all(alpha >= 0)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}





## Tverberg depth for arbitrary contexts

get_minimal_generators <- function(context,point_index,exclude_point_index=FALSE){
   model <- oofos:::compute_extent_vc_dimension(context)
   m <- nrow(context)
   n <- ncol(context) 
   zero_indexs <- which(context[point_index,]==0)
   if(length(zero_indexs)==0){print("arghh");return(NULL)} 
   model$ub <- c(rep(1,m),rep(0,n))
   model$ub[m+ zero_indexs] <- 1
   if(exclude_point_index){model$ub[point_index] <- 0}
   for(k in zero_indexs){
     i <- which(context[,k]==0)
     temp <- rep(0,m+n)
     temp[i] <- 1
     model$A <- rbind(model$A, matrix(temp,nrow=1))
     model$rhs <- c(model$rhs,1)
     model$sense <- c(model$sense,">=")
   }
   model$obj <- NULL
 
 return(model)}

 Tverberg_depth <- function(context,maximal_number_generators=100000){
 depths <- rep(0,nrow(context))
 for(point_index in seq_len(nrow(context))){
 #point_index <- 4
 
 model <- get_minimal_generators(context,point_index)
 result <- gurobi::gurobi(model,list(PoolSearchMode=2,PoolSolutions=maximal_number_generators,NumericFocus=3))
 
 if(result$status=="OPTIMAL")
 
 {m <- nrow(context)
 n <- ncol(context)
 N <- length(result$pool)
 if(N==maximal_number_generators){warning(c("Warning: not all generators computed for data point number",point_index))}
 mat <- array(0,c(N,m)) 
 for(k in (1:N)){
    mat[k,] <- (result$pool[[k]])$poolnx[(1:m)]
  }

sos <- list()
t <- 1
for(k in (1:m)){
  index <- which(mat[,k]>0.5)
  #print(length(index))
  if(length(index)>0){
  sos[[t]] <- list(type=1,index=index,weight=rep(1/length(index),length(index)))
  t <- t+1
  }
}
model <- list(A=Matrix::Matrix(0, nrow = 0, ncol = N, sparse = TRUE),obj=rep(1,N),sos=sos,lb=rep(0,N),ub=rep(1,N),vtype=rep("B",N),modelsense="max")

depths[point_index] <- gurobi::gurobi(model,params=list(NumericFocus=3))$objval
}
}
return(list(depths=depths))}


## parallelized version for Tverberg depth:


library(future.apply)
library(Matrix)

Tverberg_depth_par <- function(context, workers = future::availableCores()) {

  future::plan(future::multisession, workers = workers)

  depths <- future_sapply(
    seq_len(nrow(context)),
    function(point_index) {

      model <- get_minimal_generators(context, point_index)

      result <- gurobi::gurobi(
        model,
        params = list(
          PoolSearchMode = 2,
          PoolSolutions  = 100000000,
          NumericFocus   = 3,
          Threads        = 1      # wichtig!
        )
      )

      if (result$status != "OPTIMAL")
        return(0)

      m <- nrow(context)
      N <- length(result$pool)

      mat <- matrix(0, N, m)

      for (k in seq_len(N))
        mat[k, ] <- result$pool[[k]]$poolnx[1:m]

      sos <- vector("list", m)
      t <- 1

      for (k in seq_len(m)) {
        index <- which(mat[, k] > 0.5)

        if (length(index) > 0) {
          sos[[t]] <- list(
            type   = 1,
            index  = index,
            weight = rep(1 / length(index), length(index))
          )
          t <- t + 1
        }
      }

      sos <- sos[1:(t - 1)]

      model2 <- list(
        A          = Matrix(0, nrow = 0, ncol = N, sparse = TRUE),
        obj        = rep(1, N),
        sos        = sos,
        lb         = rep(0, N),
        ub         = rep(1, N),
        vtype      = rep("B", N),
        modelsense = "max"
      )

      gurobi::gurobi(
        model2,
        params = list(
          NumericFocus = 3,
          Threads = 1
        )
      )$objval

    },
    future.seed = TRUE
  )

  list(depths = depths)
}
