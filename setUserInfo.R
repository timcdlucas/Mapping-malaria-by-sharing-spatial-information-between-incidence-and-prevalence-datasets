# Detect correct user

Z <- function(path){
  if(Sys.info()["user"] != 'anita'){
    fullpath <- paste0('~/timz/', path)
  } else {
    fullpath <- paste0( '~/Z/', path)
  }
  if(Sys.info()["sysname"] == 'Windows'){
    fullpath <- paste0( 'Z://', path)
  }
  return(fullpath)
}
