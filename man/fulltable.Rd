\name{fulltable}
\alias{fulltable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
fulltable attaches margins and NA/NaN category to the output of table()
}
\description{
fulltable attaches margins and NA/NaN category to the output of table()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
same as table()
}
\value{
same as returned from table()
}
\author{
Hormuzd A. Katki
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
table
}
\examples{

## The function is currently defined as
function (...)
{
  ## Purpose: Add the margins automatically and don't exclude NA/NaN as its own row/column
  ##          and also add row/column titles.  Works for mixed numeric/factor variables.
  ##          For factors, the exclude option won't include the NAs as columns, that's why
  ##          I need to do more work.
  ## ----------------------------------------------------------------------
  ## Arguments: Same as for table()
  ## ----------------------------------------------------------------------
  ## Author: Hormuzd Katki, Date:  5 May 2006, 19:45

  # This works for purely numeric input, but not for any factors b/c exclude=NULL won't
  # include NAs for them.
  #   return(addmargins(table(...,exclude=NULL)))

  ##
  # Factors are harder.  I have to reconstruct each factor to include NA as a level
  ##
  
  # Put everything into a data frame
  x <- data.frame(...)

  # For each factor (in columns), get the raw levels out, reconstruct to include NAs
  # That is, if there are any NAs -- if none, add it as a level anyway
  for (i in 1:dim(x)[2]) {
    if ( is.factor(x[,i]) )
      if ( any(is.na(x[,i])) )
        x[,i] <- factor(unclass(x[,i]),labels=c(levels(x[,i]),"NA"),exclude=NULL)
      else
        levels(x[,i]) <- c(levels(x[,i]),"NA")
  }

  # Make table with margins.  Since NA is a level in each factor, they'll be included
  return(addmargins(table(x,exclude=NULL)))

  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
