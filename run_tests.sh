PKG=src/valiant

pycodestyle --ignore=E501,W504,E266 ${PKG}

pytest --cov=${PKG} tests
