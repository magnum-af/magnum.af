# $1 ... abs path to insert_header.html
# $2 ... abs path to insert_footer.html

doxygen -w html header.html footer.html customdoxygen.css
rm customdoxygen.css
sed -i "/<head>/ r $1" header.html
sed -i "/<\/body>/ r $2" footer.html # searches for </body>
