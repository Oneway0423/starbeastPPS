seqgenrunner <- function(options, tree, seqgendir) {
    commandline <- c("cd", seqgendir, ";", "./seq-gen", "-q", options)
    temp.file.seqgen <- tempfile("seqgen.")
    newick.tree <- write.tree(tree, digits = 18)
    write.tree(tree, file = temp.file.seqgen)
    commandline <- c(commandline, temp.file.seqgen)
    commandline <- paste(commandline, collapse = " ")
    out <- system(commandline, intern = TRUE)
    unlink(temp.file.seqgen)
    class(out) <- "seqgen"
    out
} 
