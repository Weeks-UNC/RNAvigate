RNAvigate style guide
---------------------

RNAvigate adheres to the Google Python Style Guide found
[here](https://google.github.io/styleguide/pyguide.html). What follows is a
largely abridged version with some of my own additions regarding this project.

In general, keep it as simple as possible, human readable, and consistent.

The biggest difference is the docstring style. I like to keep these very
legible for non-python users, especially the high-level interface in
rnavigate.py and plotting_functions.py so that they are more readable for
non-python users. Another difference is the commit messages, which I'll address
first.

TOC:
- [Commit messages](#commit-messages)
- [Linting](#linting)
- [Imports](#imports)
- [Type Annotated Code](#type-annotated-code)
- [Python Style Rules](#python-style-rules)
    - [long function calls](#long-function-calls)
    - [Docstrings](#docstrings)
    - [Modules](#modules)
    - [Functions and Methods](#functions-and-methods)
    - [Classes](#classes)
- [Naming](#naming)
- [Main](#main)

## Commit messages

Commit messages are generally limited to 50 characters, which is very short.
These short-hand codes help and should be used in all commits. These codes
should be in parentheses, and parentheses should not be used elsewhere.

(+) Means that a feature was added.
(-) Means a feature or name was removed.
(~) Means a behavior was altered.
(fix) Means a bug was fixed

For large changes, a longer message is needed. I don't use a specific format
for this. It should be concise and clear, and describe how and why the changes
were made.

## Linting

Run pylint over your code using the .pylintrc in the top-level directory. Fix
these issues regularly.

## Imports

Use import statements for packages and modules, not for classes or functions.

Use `import x` for importing packages and modules.
Use `from x import y` where x is the package prefix and y is the module name with no prefix.
Use `from x import y` as z in any of the following circumstances:
- Two modules named y are to be imported.
- y conflicts with a top-level name defined in the current module.
- y conflicts with a common parameter name that is part of the public API (e.g., features).
- y is an inconveniently long name.
- y is too generic in the context of your code (e.g., from storage.file_system import options as fs_options).
Use `import y as z` only when z is a standard abbreviation (e.g., import numpy as np).

Import each module using the full pathname location of the module.

Yes:
```python
# Reference absl.flags in code with the complete name (verbose).
import absl.flags
from doctor.who import jodie

_FOO = absl.flags.DEFINE_string(...)
```

No:
```python
# Unclear what module the author wanted and what will be imported.  The actual
# import behavior depends on external factors controlling sys.path.
# Which possible jodie module did the author intend to import?
import jodie
```

## Type Annotated Code

Most of RNAvigate is not type-hinted, but they are fine to use. I think these
are most readable when they are one argument per line with the closing
parenthesis also on it's own line.

# Python Style Rules

Do not terminate your lines with semicolons, and do not use semicolons to put
two statements on the same line.

Maximum line length is 80 characters, except:
- long imports, to keep them explicit
- long URLs, these should get their own line

## long function calls

Function calls occasionally wrap to multiple lines to avoid the 80 character
limit. I prefer this layout with a closing parenthesis indented on it's own
line. Function calls under 80 characters are fine on one line.

```python
short_function(short_arg='')
a_long_function_name_with_a_few_short_arguments(
    x=1, y=2, message="hello there",
    )
a_long_function_name_with_long_arguments(
    renormalization_method="percentile",
    renormalization_bounds=(0.95, 0.99),
    )
```

Not this: (weird indentation)
```python
short_function(
    short_arg='')
a_long_function_name_with_a_few_short_arguments(x=1, y=2,
                                                message="hello there")
a_long_function_name_with_long_arguments(renormalization_method="percentile",
                                         renormalization_bounds=(0.95, 0.99))
```

## Docstrings

This is where RNAvigate diverges a bit from Google style. I decided on these
docstring styles to be as helpful as possible to non-python coders. They are
most important for the high-level API in rnavigate.py and plotting_functions.py
but I rather like them and am transitioning to using them every where.

### Modules

This section pretty much the same as the Google style.

Every file should contain the license boilerplate found in the LICENSE file.

Files start with a docstring describing the contents and usage of the module.

```python
"""A one-line summary of the module or program, terminated by a period.

Leave one blank line.  The rest of this docstring should contain an
overall description of the module or program.  Optionally, it may also
contain a brief description of exported classes and functions and/or usage
examples.

Citations:
    Copy/paste citation for any algorithms taken from published work.
    Copy/paste URLs for unpublished work
    This section is especially important in the analysis module.

Typical usage example:

  foo = ClassFoo()
  bar = foo.FunctionBar()
"""
```

### Functions and Methods
```python
    """a one line description ending in a period.

    a multiline description describing the precise behavior of the function and
    any important details about it's implementation, or a citation and link if
    an implementation of a published algorithm.

    Require arguments:
        first_argument (string)
            a string used for ____

    Optional arguments:
        second argument (integer)
            a number used for ____
            defaults to 0

    Returns a string
        This string represents ___

        For example:
        ...

    Raises:
        IOError if ...
        ValueError if ...
```

This applies to all methods, functions, generators, and properties.

A docstring is mandatory if the function is:
- part of the public API
- nontrivial size
- non-obvious logic

Use imperative style. (`"""Return rows from a Bigtable."""`)

The first word should match the python verb `Return`, `Yield`, `Print`, etc.
if that is the end behavior.

`Arguments:` followed by a list with each parameter, then the type expected in
parentheses, a line break and indent, then a description, the default value
gets its own line.


High-level API that is intended for bench scientist use should be long form and
not use python short-hand: (string) not (str), (True or False) not (bool). The
`Arguments:` section should be split into `Required arguments:` and
`Optional arguments:` separated by a double line break.

`Returns:` or `Yields:` followed by the return type and a description.

`Raises:` followed by a list of exceptions that are relevant to the interface
followed by a description. Exceptions in which an argument input does not match
the expected value do not need to be listed.

### Classes

Classes should have a docstring below the class definition describing the class.
If your class has public attributes, they should be documented here in an
Attributes section and follow the same formatting as a function’s Args section.

```python
class SampleClass:
    """Summary of class here.

    Longer class information...
    Longer class information...

    Attributes:
        likes_spam (True or False)
            whether we like SPAM or not.
        eggs (integer)
            An integer count of the eggs we have laid.
    """

    def __init__(self, likes_spam: bool = False):
        """Initializes the instance based on spam preference.

        Optional arguments:
          likes_spam (True or False)
            whether the SampleClass exhibits this preference.
        """
        self.likes_spam = likes_spam
        self.eggs = 0
```

All class docstrings should start with a one-line summary that describes what
the class instance represents. This implies that subclasses of Exception should
also describe what the exception represents, and not the context in which it
might occur. The class docstring should not repeat unnecessary information,
such as that the class is a class.

## Naming
- module_name
- package_name
- ClassName
- method_name
- ExceptionName
- function_name
- GLOBAL_CONSTANT_NAME
- global_var_name
- instance_var_name
- function_parameter_name
- local_var_name
- query_proper_noun_for_thing
- send_acronym_via_https.

avoid abbreviations. Except for:
- nt for nucleotide
- pos for position (1-indexed inclusive)
- idx for index (0-indexed, exclusive right)
- i and j for interactions or base-pair indices or nested enumerations
- in simple list comprehensions
- e for exception identifiers (except Error as e: Raise e)
- f for open file handlers in with statements
- mathematical notation in small scopes (r for radius, etc.)
  - spell out greek names (delta_t, not dt or d_t, etc.)

Always use a .py filename extension. Never use dashes.

don't include the type in variable names (id_to_name, not id_to_name_dict)

## Main
For now, RNAvigate is explicitly a tool for Jupyter and scripting, not command
line use, no `if __name__ == '__main__':` statements are necessary.
