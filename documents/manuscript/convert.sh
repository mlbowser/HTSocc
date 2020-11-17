pandoc 01_title.md \
02_author_information.md \
03_abstract.md \
04_introduction.md \
05_study_site.md \
06_methods.md \
07_results.md \
08_discussion.md \
10_acknowledgments.md \
11_references.md \
 --filter pandoc-citeproc -o manuscript_01.md -t gfm -s

pandoc 13_figures_and_figure_captions.md -o manuscript_02.md -t gfm -s

pandoc 01_title.md \
manuscript_01.md \
 manuscript_02.md \
 --reference-doc=custom-reference.docx -o manuscript.docx

rm manuscript_01.md
rm manuscript_02.md
