.First.lib <- function(lib, pkg)
{
  if (!require(survival)) error("'survival' is essential!")
  library.dynam( "eha", pkg, lib )
}
