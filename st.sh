# 9-12-2019 JHZ

git add README.md
git commit -m "README"
git add BiSeq hearing
git commit -m "Code"
module load gcc/5
R --no-save <<END
  rmarkdown::render("P.Rmd", c("html_document", "pdf_document"))
END
module unload gcc/5
git add P.Rmd P.html P.pdf
git commit -m "P values"
git add st.sh
git commit -m "st.sh"
git push
