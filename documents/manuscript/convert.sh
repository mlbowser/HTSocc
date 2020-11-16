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
 --filter pandoc-citeproc --reference-doc=custom-reference.docx -o manuscript_01.docx

pandoc 13_figures_and_figure_captions.md --reference-doc=custom-reference.docx -o manuscript_02.docx

pandoc manuscript_01.docx manuscript_02.docx --reference-doc=custom-reference.docx -o manuscript.docx

rm manuscript_01.docx
rm manuscript_02.docx
