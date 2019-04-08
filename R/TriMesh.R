#' Create Triangles Mesh in 2D Domains
#'
#' This function triangulates the polygonal domain by using Delaunay Triangulation.
#' @importFrom graphics plot
#' @importFrom tripack tri.mesh triangles
#' @importFrom pracma numel isempty meshgrid
#' @importFrom fdaPDE create.MESH.2D
#' @param Pt A two by \code{N} matrix which indicates the outer boundry points of a 2D region.
#' @param n An integer parameter controlling the fineness of the triangulation
#' and subsequent triangulation. As n increases the fineness increases. Usually, \code{n = 8} seems to be a
#' good choice.
#' @param H A list of vertices that are the inner boundary points,
#' default set to '\code{NULL}' if there is no holes.
#'
#' @return
#' \item{V}{an \code{N} by two matrix that lists vertices with the \code{i}th
#' row storing in Cartesian coordinates for the
#' \code{i}th vertex. \code{N} is the number of vertices.}
#' \item{Tr}{a \code{K} by three matrix that each row represents one triangle.
#' All the elements are the integers that stand for the indices of vertices in \code{V}.}
#' @details In the function, we firstly get grid points inside and on the boundary of
#' the polygon with extreme points \code{Pt} and interior holes defined by \code{H}. Then delaunay triangulation
#' is used to generate triangulations by using the grid points.
#' And lastly we delete triangles within the holes or outside the boundary of the region.
#'
#' @examples
#' # rectangular domain
#' bb=rbind(c(0,0),c(1,0),c(1,1),c(0,1))
#' VT=TriMesh(bb,2)
#'
#' # irregular domain
#' data("horseshoe")
#' TriMesh(horseshoe,n=8)
#'
#' # region with holes
#' data("BMP")
#' TriMesh(BMP$bound,n=16,list(as.matrix(BMP$H1),as.matrix(BMP$H2)))
#'
#' data("mymontreal")
#' plot(mymontreal$bound[,1],mymontreal$bound[,2])
#' TriMesh(mymontreal$bound,5)
#' TriMesh(mymontreal$bound,18,list(mymontreal$H1,mymontreal$H2))
#' @export

TriMesh <- function(Pt,n,H=NULL) {
  X <- gridpoly(Pt,n,H)$X
  Y <- gridpoly(Pt,n,H)$Y
  X[abs(X)<1e-12] <- 0
  Y[abs(Y) < 1e-12] <- 0
  tmp <- cbind(X,Y)
  tmp <- unique(tmp)
  X <- tmp[,1]
  Y <- tmp[,2]
  tt <- tri.mesh(X,Y)
  Tr <- triangles(tt)[,1:3]
  Tr <- as.matrix(Tr)
  V <- cbind(X,Y)
  Tr <- del_tri(Pt,V,Tr,1)$NewTr
  Tr <- as.matrix(Tr)
  if (!isempty(H)) {
    for (i in 1:numel(H)) {
      Tr <- del_tri(H[[i]],V,Tr,-1)$NewTr
    }
  }
  tol <- 1e-12
  area <- c()
  Tr <- as.data.frame(Tr)
  for (i in 1:nrow(Tr)) {
    area <- c(area, triarea(V[Tr[i,1],],V[Tr[i,2],],V[Tr[i,3],]))
  }
  dT <- which(area<tol)
  if (length(dT) > 0) {
    Tr<-Tr[-dT,]
  }
  triplot  <-  create.MESH.2D(nodes = V,triangles = Tr)
  plot(triplot)
  return(list(V = V,Tr = Tr))
}
