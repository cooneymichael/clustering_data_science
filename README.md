# K-Clustering Algorithm

---



TODO:

- [x] Rewrite clustering algorithm in an object oriented way
- [] Rethink clustering algorithm to prevent occasional stalls (infinite loops)
  - this is rare unless outer_index is set to a large number, which would be preferred
- [] Rethink how radial positions are calculated
  - I think they're calculated naively, and they're always the same
- [] Review algorithm and potentially improve time complexities
- [] Consider writing class for radar charts/graphs in general



---



## This code was written by Michael Cooney (MC509119@ohio.edu)

The data set for this class can be found [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64881) as a text file.  It was [published in the journal nature](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5366070/#!po=5.27638) in 2017.  It is the result of a new genome mapping technique developed by researches, called Genomic Architecture Mapping (GAM).


This code is written primarily in Python.  I thought about using R, and gave a small try, but decided to stay with Python because I know it better.  I've kept the R folder here in case I decide to try my hand at it later.


The Images folder contains screenshots of various charts and graphs which might be interesting to look at.  One of my goals is to no longer have to take screenshots to view multiple graphs together, but instead write the code so that all graphs are displayed at the same time.


The OldCode folder contains my first attempts at this class.  There you will find the code which I'm currently trying to refactor.  It's mostly written as a series of scripts.  Several of the files are almost redundant because I hadn't created a git repository to act a version control system yet.  Hopefully soon, that will all be written in a more object oriented fashion.  You'll also notice that the scripts won't run because I reorganized the file system after making the repo.


CSVs is where you will find all the csv files for this project.  At the time of writing, there is one (important) file missing, but I should have the code to regenerate it and add it to the directory.


Requirements.txt is the list of requirements needed for the python virtual environment.
