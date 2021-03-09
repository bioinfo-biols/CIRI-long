pip install wheel setuptools
pip install twine
pip install keyring
python setup.py sdist bdist_wheel
twine check dist/*
twine upload dist/*
