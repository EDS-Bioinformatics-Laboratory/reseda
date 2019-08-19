# https://github.com/igraph/python-igraph/issues/88

The trick to apply the following fix for Python 3.7.3
is to download the release first with

pip download python-igraph
Unpack edit the file igraph/drawing/__init__.py and apply the fix, then install

result = io.getvalue()
return result.decode("utf-8")

pip install  python-igraph-0.7.1.post6.tar.gz
