RNAvigate style guide
---------------------

RNAvigate relies on Black to autoformat code, isort for import sorting, and uses Numpy
style docstrings. Conversion to Numpy style docstrings is work in progress.

Linting is performed with Pylint using the .pylintrc configuration file in the main
directory. This file is taken from Google's pylintrc file
([download](https://google.github.io/styleguide/pylintrc)), with the following
modifications:
- Pylint recommends always using `encoding=` with `open()`. I turned this off.
  - `disable=unspecified-encoding`
- indent string is 4 spaces (Google uses 2). Black enforces 4 spaces anyway.
  - `indent-string='    '`

In general, keep things as simple as possible, human readable, and consistent.

TOC:
- [Commit messages](#commit-messages)
- [Linting](#linting)
- [Type Annotated Code](#type-annotated-code)
- [Long lines](#long-lines)
  - [long imports](#long-imports)
  - [long function calls](#long-function-calls)
- [Naming](#naming)
- [Main](#main)
- [Numpy style docstrings](#numpy-style-docstrings)

## Commit messages

Commit messages are generally limited to 50 characters, which is very short.
These short-hand codes help and should be used in all commits. These codes
should be in parentheses, and parentheses should not be used elsewhere.

(+) means that a feature or name was added or exposed.
(-) means a feature or name was removed.
(~) means a behavior was altered.
(fix) means a bug was fixed
(doc) means only changes are to docstrings and in-line comments
(rtd) means only changes are to documentation website or configuration

For large changes, a longer message is needed. I don't use a specific format
for this. It should be concise and clear, and describe how and why the changes
were made.

## Type Annotated Code

Most of RNAvigate is not type-hinted, but they are fine to use. I think these
are most readable when they are one argument per line with the closing
parenthesis also on it's own line. similar to
[long function calls](#long-function-calls).

# Long lines

Maximum line length is 88 characters. This is enforced by Black, except:
- long imports, to keep them explicit
- long URLs, these should get their own line

## Long imports

In general, only modules and submodules should be imported. This usually means
that imports will fit on one line. In cases where multiple submodules are
imported, use grouping to wrap the line. You'll see this applied often in
RNAvigate's init files, where many names are imported explicitly.

```python
from a_really_long_module_name.with_multiple_submodules import (
    submodule1,
    submodule2,
    submodule3,
    )
```

## Long function calls

Black will automatically format long function calls appropriately.

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
- send_acronym_via_https

avoid abbreviations. Except for:
- nt for nucleotide
- pos for position (1-indexed inclusive)
- idx for index (0-indexed, exclusive right)
- i and j for interactions or base-pair indices or nested enumerations
- anything in simple list comprehensions, as long as it makes sense
- e for exception identifiers (except Error as e: Raise e)
- f for open file handlers in with statements
- mathematical notation in small scopes (r for radius, etc.)
  - spell out greek names (delta_t, not dt or d_t, etc.)

Always use a .py filename extension. Use underscores, never dashes.

don't include the type in variable names (id_to_name, not id_to_name_dict)

## Main

For now, RNAvigate is explicitly a tool for Jupyter and scripting, not command
line use, no `if __name__ == '__main__':` statements are necessary.

## Numpy style docstrings

This is a cheat sheet. For the full version, refer to the [numpydoc style guide][].

[numpydoc style guide]: https://numpydoc.readthedocs.io/en/latest/format.html

### Modules

    """Docstring for the example.py module.
    
    Modules names should have short, all-lowercase names.  The module name may
    have underscores if this improves readability.
    
    Every module should have a docstring at the very top of the file.  The
    module's docstring may extend over multiple lines.  If your docstring does
    extend over multiple lines, the closing three quotation marks must be on
    a line by itself, preferably preceded by a blank line.
    
    """

### Functions

    """Summarize the function in one line.
    
    Several sentences providing an extended description. Refer to
    variables using back-ticks, e.g. `var`.
    
    Parameters
    ----------
    var1 : array_like
        Array_like means all those objects -- lists, nested lists, etc. --
        that can be converted to an array.  We can also refer to
        variables like `var1`.
    var2 : int
        The type above can either refer to an actual Python type
        (e.g. ``int``), or describe the type of the variable in more
        detail, e.g. ``(N,) ndarray`` or ``array_like``.
    *args : iterable
        Other arguments.
    long_var_name : {'hi', 'ho'}, optional
        Choices in brackets, default first when optional.
    
    Returns
    -------
    type
        Explanation of anonymous return value of type ``type``.
    describe : type
        Explanation of return value named `describe`.
    out : type
        Explanation of `out`.
    
    Other Parameters
    ----------------
    only_seldom_used_keyword : int, optional
        Infrequently used parameters can be described under this optional
        section to prevent cluttering the Parameters section.
    **kwargs : dict
        Other infrequently used keyword arguments. Note that all keyword
        arguments appearing after the first parameter specified under the
        Other Parameters section, should also be described under this
        section.
    
    Raises
    ------
    BadException
        An explaination of the bad thing you did.
    
    See Also
    --------
    A list of related functions that are otherwise not easy to find.
    
    Notes
    -----
    Notes about the implementation algorithm (if needed).
    
    References
    ----------
    Citations of the relevant literature, e.g. [1]_.
    
    .. [1] O. McNoleg, "The integration of GIS, remote sensing,
        expert systems and adaptive co-kriging for environmental habitat
        modelling of the Highland Haggis using object-oriented, fuzzy-logic
        and neural-network techniques," Computers & Geosciences, vol. 22,
        pp. 585-588, 1996.
    
    Examples
    --------
    These are written in doctest format, and should illustrate how to use the function.
"""

### Classes

Classes use the same sections as above with the following additions: Attributes and
Methods. Methods are typically not necessary, unless there are many

    """Include the relevent sections mentioned above.
    
    ...
    
    Attributes
    ----------
    attr1 : float
        description of attr1.
    
    Methods
    -------
    do_thing1 (arg1='string')
        do the thing with `arg1`.
    do_thing2 (arg1=1.0)
        do the other thing with `arg1`.
    """
