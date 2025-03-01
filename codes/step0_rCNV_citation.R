
# Ref
citation("xml2")
capture.output(utils:::print.bibentry(citation("xml2"), style = "Bibtex"),
               file = "xml2.bib")

citation("agricolae")
capture.output(utils:::print.bibentry(citation("agricolae"), style = "Bibtex"),
               file = "agricolae.bib")

citation("rcompanion")
capture.output(utils:::print.bibentry(citation("rcompanion"), style = "Bibtex"),
               file = "rcompanion.bib")


capture.output(utils:::print.bibentry(citation("ggplot2"), style = "Bibtex"),
               file = "ggplot2.bib")
capture.output(utils:::print.bibentry(citation("gghalves"), style = "Bibtex"),
               file = "gghalves.bib")
capture.output(utils:::print.bibentry(citation("ggsci"), style = "Bibtex"),
               file = "ggsci.bib")
capture.output(utils:::print.bibentry(citation("vegan"), style = "Bibtex"),
               file = "vegan.bib")
citation()
capture.output(utils:::print.bibentry(citation(), style = "Bibtex"),
               file = "R.bib")

capture.output(utils:::print.bibentry(citation("lmerTest"), style = "Bibtex"),
               file = "lmerTest.bib")

capture.output(utils:::print.bibentry(citation("MuMIn"), style = "Bibtex"),
               file = "MuMIn.bib")



