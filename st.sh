# 8-7-2020 JHZ

function P()
{
module load gcc/5
R --no-save <<END
  rmarkdown::render("P.Rmd", c("html_document", "pdf_document"))
END
module unload gcc/5
}

git add README.md
git commit -m "README"
git add BiSeq hearing HUA_methylation_analysis/scripts* HUA_methylation_analysis/README.md HUA_methylation_analysis/BiSeq_1.28.1.tar.gz
git commit -m "Code"
git add P.Rmd P.html P.pdf
git commit -m "P values"
git add st.sh
git commit -m "st.sh"
git push
