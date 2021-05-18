rm dist/*
pip install wheel setuptools
pip install twine
python setup.py sdist bdist_wheel
twine check dist/*
twine upload --repository-url https://test.pypi.org/legacy/ dist/*