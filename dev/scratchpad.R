myfunct <- function(object) {
  object_names <- substitute(object)

  make_named_list(object, object_names)
}

x <- 1
y <- 3
myfunct(list(glm(Y ~ A, data = data), x, y, "c"))
myfunct(glm(Y ~ A, data = data))
