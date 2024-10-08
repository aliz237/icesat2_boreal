
# Table of Contents

1.  [Problems with adding notebooks to a git repository](#orge9d1467)
2.  [How to version control notebooks in a git repository?](#orga6b9d6a)
3.  [Example Workflow](#org3545597)
    1.  [Install `jupytext`](#org70ee133)
    2.  [Convert to `.py`](#org4f2b3c0)
    3.  [The `.py` version updates are automatic](#org0783d28)
    4.  [Convert to `.ipynb`](#org8bf954b)
4.  [More automated workflow](#org4d427fd)
5.  [Reference](#orga7be89d)


<a id="orge9d1467"></a>

# Problems with adding notebooks to a git repository

Notebooks contain code, metadata, and output. We only want to version control the code part of it.
A few reasons why we don't want to directly add a notebook to git:

-   `git diff` of a notebook is not informative or hard to read.
-   We don't really want to track history of outputs/cell executaions.
-   They bloat the git history pretty quickly with information that's not really useful.


<a id="orga6b9d6a"></a>

# How to version control notebooks in a git repository?

We can fairly easily pair a `.ipynb` notebook with a `.py` file such that the original notebook can be fully recovered
from the `.py` equivalent. This way, we can just version control the `.py` file and not worry about adding the `ipynb` in the repository.
One tool that makes this possible is `jupytext`.


<a id="org3545597"></a>

# Example Workflow

Suppose I have been working on a `Review_maps.ipynb` in my local repository, and it's ready to be commited to the repo.
Here are the steps to convert it to the equivalent `.py` representation, and check that in instead.


<a id="org70ee133"></a>

## Install `jupytext`

    pip install jupytext


<a id="org4f2b3c0"></a>

## Convert to `.py`

    jupytext --set-format ipynb,py:light Review_maps.ipynb

This creates a new `Review_maps.py` in the same directory which can then be added to a commit and pushed.


<a id="org0783d28"></a>

## The `.py` version updates are automatic

Suppose I updated the original notebook and added a few markown or code cells. As soon as the notebook is saved,
all the updates, *minus the cell outpus and meta*, will be reflected to the `.py` version. Basically, we never need to manually touch the `.py` files.


<a id="org8bf954b"></a>

## Convert to `.ipynb`

The original notebook was never added to git. Suppose the local working directory was cleaned (e.g with `git clean` or the notebook was accidently removed.
It's possible to recover it from the `.py` version we have in git. The following command recreates the notebook and stores it as `test.ipynb`:

    jupytext --to ipynb Review_maps.py -o test.ipynb

or, if want it stored as `Review_maps.ipynb`:

    jupytext --to ipynb Review_maps.py


<a id="org4d427fd"></a>

# More automated workflow

We can simplify the workflow a bit by creating a `jupytext.toml` config file in the root of the local repository.
Put the following two lines in the `jupytext.toml`:

    # Pair ipynb notebooks to py:light text notebooks
    formats = "ipynb,py:light"

Next time a new notebook is created, it will be automatcally paired with a `.py` file, without running the earlier commands above.
However, this doesn't pair any of the current notebooks in the repository. If a number of notebooks need to paired in one go,
we can loop through the notebooks and run `jupytext --set-format ipynb,py:light path/to/each/notebook.ipynb`.


<a id="orga7be89d"></a>

# Reference

For more information on `jupytext` and formats other than `light` used above see [jupytext](https://jupytext.readthedocs.io/en/latest/index.html).

