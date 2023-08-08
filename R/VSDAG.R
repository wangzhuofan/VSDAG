emvs = function(seed,n,p,x){
  set.seed(seed)
  v0 = 0.01
  v1 = 1
  beta = matrix(0,p,p)
  #diag(beta) = 0
  #beta = beta-beta
  beta[3,1] = 0
  beta[3,2] = 0
  #beta = rep(1,p)
  sigma2 = rep(1,p)
  test1 = gnlearn::pc(as.data.frame(t(x)),to = "adjacency")
  adj = matrix(as.numeric(test1!=0),p,p)
  theta = sum(adj)/(p*(p-1))
  nu = 1
  lambda = 1
  a = 1
  b = 1
  gamma = matrix(0,p,p)
  #gamma[upper.tri(gamma,diag = TRUE)] <- 0
  c = 0.0001
  eta=0.01
  iters = 10

  # gamma = B0
  # beta= B0
  for (it in 1:iters) {


    loss.f = function(beta){
      temp = t(diag(p)-beta)%*%diag(1/sigma2)
      Theta = temp%*%(diag(p)-beta)
      Theta = Theta+0.01*diag(p)
      # loss = n*(log(det(Theta)))-sum(diag(t(x)%*%Theta%*%x))
      loss = sum(mvtnorm::dmvnorm(t(x),rep(0,p),solve(Theta),log = TRUE))
      d = diag(1/sigma2)%*%((1-gamma)/v0+gamma/v1)
      # loss = loss+sum(beta^2*d)/(2*sigma2)
      loss = loss-sum(beta^2*d)/2
      return(-loss)
    }

    loss.g = function(beta){
      d = diag(1/sigma2)%*%((1-gamma)/v0+gamma/v1)
      beta = matrix(beta,p,p)
      #if(det(t(diag(p)-beta))>1e-10)
      res = n*(diag(p)-beta)%*%solve(t(diag(p)-beta)%*%(diag(p)-beta)+0.01*diag(p))-diag(1/sigma2)%*%(diag(p)-beta)%*%x%*%t(x)+beta*d
      return(res)
    }

    lb = matrix(-100,p,p)
    diag(lb) = -0.0001
    ub = matrix(100,p,p)
    diag(ub) = 0.00001
    beta = optim(beta, loss.f,gr=loss.g, method = "L-BFGS-B",lower = lb,upper = ub)$par
    #beta = optim(beta, loss.f,gr=loss.g, method = "L-BFGS-B",lower = lb,upper = ub,control = list(trace=1,REPORT=1,maxit=100))$par
    #beta = optim(beta, loss.f, method = "L-BFGS-B",lower = lb,upper = ub)$par
    # sigma2 = norm(x-beta%*%x,"F")^2+norm((beta*gamma)^2/v1,"F")^2+nu*lambda
    # sigma2 = sigma2/(p^2+n*p+nu)
    y = (diag(p)-beta)%*%x
    sigma2 = (rowSums(y^2)+rowSums(beta^2*((1-gamma)/v0+gamma/v1)))/(n+p)

    ap = theta*dnorm(beta,0,sqrt(sigma2*v1))
    bp = (1-theta)*dnorm(beta,0,sqrt(sigma2*v0))
    for (i in 1:p) {
      for (j in 1:p) {
        if(j!=i){
          gamma[i,j] = as.numeric(ap[i,j]>bp[i,j])
          if(gamma[i,j]==1){
            temp = igraph::graph_from_adjacency_matrix(gamma)
            if(!igraph::is_dag(temp)) gamma[i,j]=0
          }
        }
      }
    }
    theta = (sum(gamma)+a-1)/(a+b+p^2-p-2)
    if(theta<c)
      theta=c
  }

  return(list(gamma=gamma,beta=beta))
}

