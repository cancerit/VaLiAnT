PKG=src/valiant

pycodestyle --ignore=E501,W504,E266 ${PKG}

stmt=$(pytest --cov=${PKG} 2>/dev/null | grep TOTAL | tr -s ' ' | cut -f2 -d' ')

echo "Statements: ${stmt}"

pytest tests
